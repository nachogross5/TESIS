# ============================================================
# CRANIOPLAN - BLOQUE F v5 (UNIFICADO)
#
# Une el generador de cintas de corte y el corte en dominio de voxels
# en un solo archivo y un solo comando.
#
# FLUJO:
#   Bloque A+B v5 -> confirmar_craneo() -> enviar_a_planner()
#   -> dibujar una o mas curvas cerradas sobre Craneo_Final
#      (Shortest distance on surface + Constrain to Model)
#   -> cargar este archivo -> cortar()
#
# COMANDOS:
#   cortar()                    -> genera las cintas y ejecuta el corte completo
#   generar_cintas()            -> solo genera las cintas, para inspeccionarlas antes
#   cortar(regenerar=False)     -> corta usando las cintas que ya estan en la escena
#   exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight
#
# CAMBIOS v5 - CLASIFICACION UNIVERSAL (sin heuristicas de tamano ni forma):
#   Regla 1 (topologia): antes de cortar se etiquetan las islas del hueso.
#     Cada fragmento posterior sabe de que isla venia. Si una isla se partio,
#     TODOS sus pedazos los creo el corte y se conservan siempre, sin filtro
#     de volumen. Las islas que no se partieron son hueso preexistente.
#   Regla 2 (tapon): para saber si un fragmento es interior de una curva se
#     genera la superficie minima que tapa esa curva cerrada, se extruye en
#     prisma y se rasteriza. Un fragmento es interior si mas de la mitad de
#     sus voxels caen dentro. Funciona con curvas concavas, tapas partidas
#     por una sutura, o cualquier cantidad de pedazos.
#   El filtro de volumen queda SOLO para descartar islas preexistentes
#   diminutas (ruido de segmentacion que el Bloque A dejo pasar).
#
# CAMBIOS v4:
#   - Los fragmentos INTERIORES de cada curva (las tapas de hueso que se levantan
#     en la cirugia) se identifican por pertenencia geometrica a la curva y se
#     conservan SIEMPRE, sin importar su volumen. Antes el filtro de tamano se
#     los comia: una tapa de hueso pediatrico fino pesa ~0.1-0.5cm3.
#   - El filtro VOLUMEN_MINIMO_CM3 ahora solo se aplica a lo que queda AFUERA
#     de todas las curvas (ruido real de segmentacion).
#   - Cada cinta lleva el nombre de su curva: HerramientaCorte_CC, etc.
#     Antes todas se llamaban igual y eran indistinguibles en el panel Data.
#   - Reporte final con la clasificacion de cada pieza.
# ============================================================

import numpy as np
import vtk
from scipy import ndimage


# ==================== CONFIGURACION ====================

# --- Geometria del corte ---
ANCHO_CORTE_MM         = 1.5    # kerf: espesor de hueso que se pierde. Piso ~2 voxels.
PROFUNDIDAD_CORTE_MM   = 15.0   # centrada en la superficie: +-7.5mm sobre la normal
RADIO_BUSQUEDA_NORMAL  = 3.0    # radio para promediar normales de la superficie
MUESTRAS_POR_SEGMENTO  = 5      # muestras a lo largo de cada segmento entre puntos de control

# --- Nodos ---
NOMBRE_MODELO_CRANEO   = "Craneo_Final"       # Model node (para calcular normales)
NOMBRE_SEGMENTACION    = "Craneo_Automatico"  # nodo de segmentacion del Bloque A
NOMBRE_SEGMENTO_HUESO  = "Craneo_Final"       # segmento dentro de ese nodo
PREFIJO_HERRAMIENTAS   = "HerramientaCorte"
NOMBRE_SEG_TRABAJO     = "CranioPlan_Corte"
NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"

# --- Procesamiento ---
VOLUMEN_MINIMO_CM3     = 0.10   # solo descarta islas PREEXISTENTES diminutas (ruido)
PISO_ABSOLUTO_CM3      = 0.005  # piso duro: por debajo de esto es ruido numerico
FRACCION_INTERIOR      = 0.5    # fraccion de voxels dentro del tapon para ser interior
MAX_CAVIDAD_MM3        = 2000.0 # cavidades internas mayores NO se rellenan
SUAVIZADO_EXPORT       = "0.3"  # suavizado al convertir a malla
UMBRAL_LIMPIEZA_MALLA  = 0.10


# ============================================================
# PARTE 1 - GENERACION DE LAS CINTAS DE CORTE
# ============================================================

def _normal_promedio_en_punto(p, radio, locator, normalesArray):
    idList = vtk.vtkIdList()
    locator.FindPointsWithinRadius(radio, p, idList)
    if idList.GetNumberOfIds() == 0:
        idList.InsertNextId(locator.FindClosestPoint(p))
    acumulado = np.zeros(3)
    for j in range(idList.GetNumberOfIds()):
        acumulado += np.array(normalesArray.GetTuple3(idList.GetId(j)))
    norma = np.linalg.norm(acumulado)
    if norma < 1e-6:
        return np.array(normalesArray.GetTuple3(locator.FindClosestPoint(p)))
    return acumulado / norma


def _normal_promedio_segmento(p1, p2, muestras, radio, locator, normalesArray):
    acumulado = np.zeros(3)
    for k in range(muestras):
        t = k / float(muestras - 1) if muestras > 1 else 0.5
        acumulado += _normal_promedio_en_punto(
            p1 + t * (p2 - p1), radio, locator, normalesArray)
    norma = np.linalg.norm(acumulado)
    return acumulado / norma if norma > 1e-6 else acumulado


def _construir_cinta(curvaNode, craneoConNormales, locator, normalesArray, nombreSalida):
    """Construye una cinta cerrada perpendicular a la superficie a lo largo de la curva."""
    normalesArrayLocal = normalesArray

    nPuntos = curvaNode.GetNumberOfControlPoints()
    if nPuntos < 3:
        print(f"  '{curvaNode.GetName()}': menos de 3 puntos de control, se omite.")
        return None

    puntosCurva = []
    for i in range(nPuntos):
        p = [0.0, 0.0, 0.0]
        curvaNode.GetNthControlPointPositionWorld(i, p)
        puntosCurva.append(np.array(p))

    # Una normal por segmento (entre punto i y punto i+1), loop cerrado
    normalesSegmento = []
    for i in range(nPuntos):
        p1 = puntosCurva[i]
        p2 = puntosCurva[(i + 1) % nPuntos]
        normalesSegmento.append(_normal_promedio_segmento(
            p1, p2, MUESTRAS_POR_SEGMENTO, RADIO_BUSQUEDA_NORMAL,
            locator, normalesArrayLocal))

    # Normal de union en cada punto: promedio de las 2 planchas que se tocan ahi
    normalesJunta = []
    for i in range(nPuntos):
        suma = normalesSegmento[i - 1] + normalesSegmento[i]
        norma = np.linalg.norm(suma)
        normalesJunta.append(suma / norma if norma > 1e-6 else normalesSegmento[i])

    # Tangente en cada punto (wraparound)
    tangentes = []
    for i in range(nPuntos):
        anterior = puntosCurva[i - 1]
        siguiente = puntosCurva[(i + 1) % nPuntos]
        t1 = puntosCurva[i] - anterior
        t2 = siguiente - puntosCurva[i]
        n1 = np.linalg.norm(t1)
        n2 = np.linalg.norm(t2)
        t = (t1 / n1 if n1 > 1e-6 else t1) + (t2 / n2 if n2 > 1e-6 else t2)
        normaT = np.linalg.norm(t)
        tangentes.append(t / normaT if normaT > 1e-6 else np.array([1.0, 0.0, 0.0]))

    # Direccion de ancho, con correccion de continuidad de signo (evita el twist)
    wCrudos = []
    for i in range(nPuntos):
        w = np.cross(normalesJunta[i], tangentes[i])
        normaW = np.linalg.norm(w)
        wCrudos.append(w / normaW if normaW > 1e-6 else np.array([1.0, 0.0, 0.0]))

    wCorregidos = [wCrudos[0]]
    for i in range(1, nPuntos):
        wActual = wCrudos[i]
        if np.dot(wActual, wCorregidos[i - 1]) < 0:
            wActual = -wActual
        wCorregidos.append(wActual)
    if np.dot(wCorregidos[-1], wCorregidos[0]) < 0:
        print(f"  ADVERTENCIA en '{curvaNode.GetName()}': twist impar en el loop, revisar la curva.")

    # Anillos en cada punto de control
    puntosHerramienta = vtk.vtkPoints()
    celdas = vtk.vtkCellArray()
    anillos = []
    for i in range(nPuntos):
        n = normalesJunta[i]
        w = wCorregidos[i]
        c = puntosCurva[i]
        p1 = c + (ANCHO_CORTE_MM / 2.0) * w + (PROFUNDIDAD_CORTE_MM / 2.0) * n
        p2 = c - (ANCHO_CORTE_MM / 2.0) * w + (PROFUNDIDAD_CORTE_MM / 2.0) * n
        p3 = c - (ANCHO_CORTE_MM / 2.0) * w - (PROFUNDIDAD_CORTE_MM / 2.0) * n
        p4 = c + (ANCHO_CORTE_MM / 2.0) * w - (PROFUNDIDAD_CORTE_MM / 2.0) * n
        anillos.append((
            puntosHerramienta.InsertNextPoint(p1.tolist()),
            puntosHerramienta.InsertNextPoint(p2.tolist()),
            puntosHerramienta.InsertNextPoint(p3.tolist()),
            puntosHerramienta.InsertNextPoint(p4.tolist()),
        ))

    for i in range(nPuntos):
        a = anillos[i]
        b = anillos[(i + 1) % nPuntos]
        for k in range(4):
            k2 = (k + 1) % 4
            quad = vtk.vtkPolygon()
            quad.GetPointIds().SetNumberOfIds(4)
            quad.GetPointIds().SetId(0, a[k])
            quad.GetPointIds().SetId(1, a[k2])
            quad.GetPointIds().SetId(2, b[k2])
            quad.GetPointIds().SetId(3, b[k])
            celdas.InsertNextCell(quad)

    herramientaPoly = vtk.vtkPolyData()
    herramientaPoly.SetPoints(puntosHerramienta)
    herramientaPoly.SetPolys(celdas)

    triangulador = vtk.vtkTriangleFilter()
    triangulador.SetInputData(herramientaPoly)
    triangulador.Update()

    limpiador = vtk.vtkCleanPolyData()
    limpiador.SetInputData(triangulador.GetOutput())
    limpiador.Update()
    herramientaFinal = limpiador.GetOutput()

    fe = vtk.vtkFeatureEdges()
    fe.SetInputData(herramientaFinal)
    fe.BoundaryEdgesOn()
    fe.NonManifoldEdgesOn()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.Update()
    bordesLibres = fe.GetOutput().GetNumberOfCells()

    node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", nombreSalida)
    node.SetAndObservePolyData(herramientaFinal)
    node.CreateDefaultDisplayNodes()
    node.GetDisplayNode().SetColor(1.0, 0.6, 0.0)
    node.GetDisplayNode().SetOpacity(0.5)

    estado = "cerrada" if bordesLibres == 0 else f"ABIERTA ({bordesLibres} bordes libres)"
    print(f"  '{curvaNode.GetName()}' -> {node.GetName()}  "
          f"({nPuntos} puntos, {estado})")
    return node


def generar_cintas():
    """Genera una cinta de corte por cada curva cerrada presente en la escena."""
    try:
        craneoNode = slicer.util.getNode(NOMBRE_MODELO_CRANEO)
    except Exception:
        print(f"ERROR: no encuentro el modelo '{NOMBRE_MODELO_CRANEO}'. Corre enviar_a_planner() primero.")
        return []
    if not craneoNode.IsA("vtkMRMLModelNode"):
        print(f"ERROR: '{NOMBRE_MODELO_CRANEO}' no es un Model node.")
        return []

    curvas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode')]
    abiertas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')]
    if abiertas:
        print(f"AVISO: hay {len(abiertas)} curva(s) ABIERTA(s) en la escena, se ignoran. "
              "El corte necesita curvas cerradas.")
    if not curvas:
        print("ERROR: no hay ninguna curva cerrada en la escena. Dibuja al menos una.")
        return []

    # Borrar cintas de corridas anteriores
    previas = [n for n in slicer.util.getNodesByClass('vtkMRMLModelNode')
               if n.GetName().startswith(PREFIJO_HERRAMIENTAS)]
    for n in previas:
        slicer.mrmlScene.RemoveNode(n)
    if previas:
        print(f"{len(previas)} cinta(s) previa(s) eliminada(s).")

    normalsFilter = vtk.vtkPolyDataNormals()
    normalsFilter.SetInputData(craneoNode.GetPolyData())
    normalsFilter.ComputePointNormalsOn()
    normalsFilter.ComputeCellNormalsOff()
    normalsFilter.SplittingOff()
    normalsFilter.ConsistencyOn()
    normalsFilter.AutoOrientNormalsOn()
    normalsFilter.Update()
    craneoConNormales = normalsFilter.GetOutput()
    normalesArray = craneoConNormales.GetPointData().GetNormals()

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(craneoConNormales)
    locator.BuildLocator()

    print(f"Generando cintas (ancho={ANCHO_CORTE_MM}mm, profundidad={PROFUNDIDAD_CORTE_MM}mm):")
    generadas = []
    for curva in curvas:
        nombreCinta = f"{PREFIJO_HERRAMIENTAS}_{curva.GetName()}"
        node = _construir_cinta(curva, craneoConNormales, locator,
                                normalesArray, nombreCinta)
        if node is not None:
            generadas.append(node)
    return generadas


def _tapon_de_curva(curvaNode, longitud):
    """
    Construye el 'tapon' de una curva cerrada: la superficie minima que la
    tapa, extruida en prisma a lo largo de su normal media.

    Es el reemplazo universal de la esfera envolvente de v4: el prisma sigue
    el contorno REAL de la curva, asi que funciona con curvas concavas y con
    interiores partidos en varios pedazos (por ejemplo, por una sutura).
    """
    cap = vtk.vtkPolyData()
    try:
        slicer.modules.markups.logic().GetClosedCurveSurfaceArea(curvaNode, cap)
    except Exception as e:
        print(f"  no pude generar el tapon de '{curvaNode.GetName()}': {e}")
        return None
    if cap.GetNumberOfPoints() < 3:
        print(f"  el tapon de '{curvaNode.GetName()}' salio vacio.")
        return None

    pts = np.array([cap.GetPoint(i) for i in range(cap.GetNumberOfPoints())])
    centro = pts.mean(axis=0)
    # normal media del parche por descomposicion en valores singulares
    _, _, vt = np.linalg.svd(pts - centro)
    normal = vt[2]
    normal = normal / np.linalg.norm(normal)

    # correr el parche media longitud hacia atras y extruir la longitud entera,
    # para que el prisma quede centrado en la superficie del craneo
    tr = vtk.vtkTransform()
    tr.Translate(*(-normal * longitud / 2.0))
    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetInputData(cap)
    tf.SetTransform(tr)
    tf.Update()

    ext = vtk.vtkLinearExtrusionFilter()
    ext.SetInputData(tf.GetOutput())
    ext.SetExtrusionTypeToVectorExtrusion()
    ext.SetVector(*(normal * longitud))
    ext.CappingOn()
    ext.Update()

    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(ext.GetOutputPort())
    tri.Update()
    limp = vtk.vtkCleanPolyData()
    limp.SetInputConnection(tri.GetOutputPort())
    limp.Update()
    return limp.GetOutput()


def _mascaras_de_tapones(volumeNode, longitud):
    """
    Devuelve {nombreCurva: mascara booleana} rasterizando el tapon de cada
    curva cerrada sobre la grilla del CT. Nodo temporal, se borra al salir.
    """
    curvas = slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode')
    if not curvas:
        return {}

    previo = slicer.mrmlScene.GetFirstNodeByName("CranioPlan_Tapones")
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)

    segTap = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', "CranioPlan_Tapones")
    segTap.CreateDefaultDisplayNodes()
    segTap.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    if segTap.GetDisplayNode():
        segTap.GetDisplayNode().SetVisibility(False)

    modelosTemp = []
    mapaNombre = {}
    for curva in curvas:
        poly = _tapon_de_curva(curva, longitud)
        if poly is None:
            continue
        mNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLModelNode", f"_Tapon_{curva.GetName()}")
        mNode.SetAndObservePolyData(poly)
        mNode.CreateDefaultDisplayNodes()
        mNode.GetDisplayNode().SetVisibility(False)
        modelosTemp.append(mNode)

        antes = set(_ids_de_segmentos(segTap))
        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(mNode, segTap)
        nuevos = [i for i in _ids_de_segmentos(segTap) if i not in antes]
        for i in nuevos:
            mapaNombre[i] = curva.GetName()

    mascaras = {}
    for segId, nombreCurva in mapaNombre.items():
        try:
            arr = slicer.util.arrayFromSegmentBinaryLabelmap(
                segTap, segId, volumeNode).astype(bool)
            mascaras[nombreCurva] = arr
        except Exception as e:
            print(f"  no pude rasterizar el tapon de '{nombreCurva}': {e}")

    for m in modelosTemp:
        slicer.mrmlScene.RemoveNode(m)
    slicer.mrmlScene.RemoveNode(segTap)
    return mascaras


# ============================================================
# PARTE 2 - CORTE EN DOMINIO DE VOXELS
# ============================================================

def _volumen_de_referencia(segNode):
    rol = slicer.vtkMRMLSegmentationNode.GetReferenceImageGeometryReferenceRole()
    volNode = segNode.GetNodeReference(rol)
    if volNode is None:
        vols = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
        volNode = vols[-1] if vols else None
    return volNode


def _ids_de_segmentos(segNode):
    seg = segNode.GetSegmentation()
    return [seg.GetNthSegmentID(i) for i in range(seg.GetNumberOfSegments())]


def _nuevo_editor(segNode, volumeNode):
    w = slicer.qMRMLSegmentEditorWidget()
    w.setMRMLScene(slicer.mrmlScene)
    eNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
    w.setMRMLSegmentEditorNode(eNode)
    w.setSegmentationNode(segNode)
    w.setSourceVolumeNode(volumeNode)
    return w, eNode


def _rellenar_cavidades(arr, voxel_mm3, max_mm3):
    """
    Rellena SOLO cavidades internas por debajo de max_mm3.
    Las que tocan el borde del array (exterior) y las grandes (cavidad craneal)
    quedan intactas.
    """
    fondo = ~arr
    etiquetas, n = ndimage.label(fondo)
    if n == 0:
        return arr, 0

    ids_borde = set()
    for cara in (etiquetas[0, :, :], etiquetas[-1, :, :],
                 etiquetas[:, 0, :], etiquetas[:, -1, :],
                 etiquetas[:, :, 0], etiquetas[:, :, -1]):
        ids_borde.update(np.unique(cara).tolist())
    ids_borde.discard(0)

    conteos = np.bincount(etiquetas.ravel(), minlength=n + 1)
    a_rellenar = [i for i in range(1, n + 1)
                  if i not in ids_borde and conteos[i] * voxel_mm3 <= max_mm3]

    if a_rellenar:
        arr[np.isin(etiquetas, a_rellenar)] = True
    return arr, len(a_rellenar)


def _buscar_o_crear_carpeta(nombre):
    shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    sceneItemId = shNode.GetSceneItemID()
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(sceneItemId, hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if shNode.GetItemName(itemId) == nombre:
            return shNode, itemId
    return shNode, shNode.CreateFolderItem(sceneItemId, nombre)


def _limpiar_malla(modelNode, umbral=UMBRAL_LIMPIEZA_MALLA):
    poly = modelNode.GetPolyData()
    if poly is None or poly.GetNumberOfPoints() == 0:
        return
    conn = vtk.vtkPolyDataConnectivityFilter()
    conn.SetInputData(poly)
    conn.SetExtractionModeToAllRegions()
    conn.ColorRegionsOn()
    conn.Update()
    nRegiones = conn.GetNumberOfExtractedRegions()
    if nRegiones <= 1:
        return

    arrayRegiones = conn.GetOutput().GetPointData().GetArray("RegionId")
    conteo = {}
    for i in range(arrayRegiones.GetNumberOfTuples()):
        rid = int(arrayRegiones.GetTuple1(i))
        conteo[rid] = conteo.get(rid, 0) + 1

    mayor = max(conteo.values())
    conservar = [rid for rid, c in conteo.items() if c >= mayor * umbral]

    conn.SetExtractionModeToSpecifiedRegions()
    conn.InitializeSpecifiedRegionList()
    for rid in conservar:
        conn.AddSpecifiedRegion(rid)
    conn.Update()

    limpiador = vtk.vtkCleanPolyData()
    limpiador.SetInputConnection(conn.GetOutputPort())
    limpiador.Update()
    modelNode.SetAndObservePolyData(limpiador.GetOutput())
    print(f"    limpieza de malla: {nRegiones} regiones -> {len(conservar)}")


def _control_calidad(modelNode):
    poly = modelNode.GetPolyData()
    fe = vtk.vtkFeatureEdges()
    fe.SetInputData(poly)
    fe.BoundaryEdgesOn()
    fe.NonManifoldEdgesOn()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.Update()
    cerrada = (fe.GetOutput().GetNumberOfCells() == 0)

    volumen = 0.0
    if cerrada:
        tri = vtk.vtkTriangleFilter()
        tri.SetInputData(poly)
        tri.Update()
        masa = vtk.vtkMassProperties()
        masa.SetInputConnection(tri.GetOutputPort())
        masa.Update()
        volumen = masa.GetVolume()
    return cerrada, volumen


# ============================================================
# COMANDO PRINCIPAL
# ============================================================

def cortar(regenerar=True):
    """
    Ejecuta el flujo completo: genera las cintas desde las curvas cerradas
    de la escena y corta el hueso en dominio de voxels.
    regenerar=False usa las cintas que ya estan en la escena.
    """
    # ---- 1. Entradas ----
    try:
        segOrigen = slicer.util.getNode(NOMBRE_SEGMENTACION)
    except Exception:
        print(f"ERROR: no encuentro el nodo de segmentacion '{NOMBRE_SEGMENTACION}'.")
        return

    volumeNode = _volumen_de_referencia(segOrigen)
    if volumeNode is None:
        print("ERROR: no encuentro el volumen de referencia del CT.")
        return

    segIdHueso = segOrigen.GetSegmentation().GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_HUESO)
    if not segIdHueso:
        print(f"ERROR: no hay segmento '{NOMBRE_SEGMENTO_HUESO}'. Corre confirmar_craneo() primero.")
        return

    spacing = volumeNode.GetSpacing()
    voxel_mm3 = spacing[0] * spacing[1] * spacing[2]
    piso_kerf = 2.0 * max(spacing)
    print(f"Volumen de referencia: {volumeNode.GetName()}  spacing={spacing}")
    if ANCHO_CORTE_MM < piso_kerf:
        print(f"ADVERTENCIA: ANCHO_CORTE_MM={ANCHO_CORTE_MM} esta por debajo del piso "
              f"de {piso_kerf:.2f}mm (2 voxels). El corte puede no separar las piezas.")

    # ---- 2. Cintas ----
    if regenerar:
        generar_cintas()

    herramientas = [n for n in slicer.util.getNodesByClass('vtkMRMLModelNode')
                    if n.GetName().startswith(PREFIJO_HERRAMIENTAS)]
    if not herramientas:
        print("ERROR: no hay cintas de corte en la escena.")
        return
    print(f"Herramientas de corte: {[h.GetName() for h in herramientas]}")

    # ---- 3. Copia de trabajo ----
    previo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_SEG_TRABAJO)
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)

    segTrabajo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', NOMBRE_SEG_TRABAJO)
    segTrabajo.CreateDefaultDisplayNodes()
    segTrabajo.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segTrabajo.GetSegmentation().CopySegmentFromSegmentation(
        segOrigen.GetSegmentation(), segIdHueso)

    idHueso = _ids_de_segmentos(segTrabajo)[0]
    segTrabajo.GetSegmentation().GetSegment(idHueso).SetName("Hueso")

    # ---- 3b. REGLA 1: islas del hueso ANTES de cortar ----
    mascaraOriginal = slicer.util.arrayFromSegmentBinaryLabelmap(
        segTrabajo, idHueso, volumeNode).astype(bool)
    etiquetasPadre, nPadres = ndimage.label(mascaraOriginal)
    print(f"Islas de hueso antes del corte: {nPadres}")

    # ---- 3c. REGLA 2: tapones de cada curva ----
    taponesMasc = _mascaras_de_tapones(volumeNode, 2.0 * PROFUNDIDAD_CORTE_MM)
    print(f"Tapones generados: {list(taponesMasc.keys()) if taponesMasc else 'ninguno'}")

    # ---- 4. Rasterizar las cintas ----
    idsCortes = []
    for h in herramientas:
        antes = set(_ids_de_segmentos(segTrabajo))
        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(h, segTrabajo)
        idsCortes.extend([i for i in _ids_de_segmentos(segTrabajo) if i not in antes])
    print(f"{len(idsCortes)} cinta(s) rasterizada(s) a voxels.")

    widget, editorNode = _nuevo_editor(segTrabajo, volumeNode)

    idCortes = idsCortes[0]
    if len(idsCortes) > 1:
        widget.setCurrentSegmentID(idCortes)
        widget.setActiveEffectByName("Logical operators")
        efecto = widget.activeEffect()
        for otro in idsCortes[1:]:
            efecto.setParameter("Operation", "UNION")
            efecto.setParameter("ModifierSegmentID", otro)
            efecto.self().onApply()
            segTrabajo.GetSegmentation().RemoveSegment(otro)

    # ---- 5. La resta: hueso AND NOT cortes ----
    widget.setCurrentSegmentID(idHueso)
    widget.setActiveEffectByName("Logical operators")
    efecto = widget.activeEffect()
    efecto.setParameter("Operation", "SUBTRACT")
    efecto.setParameter("ModifierSegmentID", idCortes)
    efecto.self().onApply()
    segTrabajo.GetSegmentation().RemoveSegment(idCortes)
    print("Resta aplicada.")

    # ---- 6. Separar en fragmentos ----
    widget.setCurrentSegmentID(idHueso)
    widget.setActiveEffectByName("Islands")
    efecto = widget.activeEffect()
    efecto.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
    efecto.self().onApply()
    widget.setActiveEffectByName(None)

    ids = _ids_de_segmentos(segTrabajo)
    print(f"Islas generadas por el corte: {len(ids)}")

    # ---- 7. Clasificacion universal ----
    mascaras = {}
    for segId in ids:
        mascaras[segId] = slicer.util.arrayFromSegmentBinaryLabelmap(
            segTrabajo, segId, volumeNode).astype(bool)

    # REGLA 1: de que isla original venia cada fragmento
    padreDe, hijosDe = {}, {}
    for segId, arr in mascaras.items():
        etiquetas = etiquetasPadre[arr]
        padre = int(np.bincount(etiquetas).argmax()) if etiquetas.size else 0
        padreDe[segId] = padre
        hijosDe.setdefault(padre, []).append(segId)

    # REGLA 2: de que curva es interior cada fragmento
    curvaDe = {}
    for segId, arr in mascaras.items():
        nVox = np.count_nonzero(arr)
        mejor, mejorFrac = None, 0.0
        for nombreCurva, masc in taponesMasc.items():
            if masc.shape != arr.shape:
                continue
            frac = np.count_nonzero(arr & masc) / float(nVox) if nVox else 0.0
            if frac > mejorFrac:
                mejor, mejorFrac = nombreCurva, frac
        curvaDe[segId] = mejor if mejorFrac >= FRACCION_INTERIOR else None

    # Decidir que se conserva y con que nombre
    conservados = []          # (segId, nombre, esTapa)
    contadorPorCurva = {}
    contadorResto = 0
    contadorHueso = 0

    for segId in ids:
        arr = mascaras[segId]
        vol = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        creadoPorCorte = len(hijosDe[padreDe[segId]]) > 1
        curva = curvaDe[segId]

        if vol < PISO_ABSOLUTO_CM3:
            print(f"  descartado: ruido numerico, vol={vol:.4f}cm3")
            continue
        if not creadoPorCorte and vol < VOLUMEN_MINIMO_CM3:
            print(f"  descartado: isla preexistente diminuta, vol={vol:.4f}cm3")
            continue

        if curva is not None:
            contadorPorCurva.setdefault(curva, 0)
            contadorPorCurva[curva] += 1
            nombre = f"Tapa_{curva}"
            esTapa = True
        elif creadoPorCorte:
            contadorResto += 1
            nombre = f"Resto_{contadorResto}"
            esTapa = False
        else:
            contadorHueso += 1
            nombre = f"Hueso_{contadorHueso}"
            esTapa = False
        conservados.append((segId, nombre, esTapa, curva, vol))

    # Sufijos a/b/c cuando una misma curva tiene varias tapas (sutura de por medio)
    sufijos = {}
    for curva, n in contadorPorCurva.items():
        if n > 1:
            sufijos[curva] = 0
    finales = []
    for segId, nombre, esTapa, curva, vol in conservados:
        if esTapa and curva in sufijos:
            letra = chr(ord('a') + sufijos[curva])
            sufijos[curva] += 1
            nombre = f"{nombre}_{letra}"
        finales.append((segId, nombre, esTapa))

    if not finales:
        print("ERROR: no sobrevivio ningun fragmento.")
        slicer.mrmlScene.RemoveNode(editorNode)
        return

    # ---- 8. Rellenar cavidades y nombrar ----
    conservados = []
    for segId, nombre, esTapa in finales:
        arr, nRell = _rellenar_cavidades(mascaras[segId], voxel_mm3, MAX_CAVIDAD_MM3)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            arr.astype(np.uint8), segTrabajo, segId, volumeNode)
        volf = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        segTrabajo.GetSegmentation().GetSegment(segId).SetName(nombre)
        conservados.append((segId, nombre, esTapa))
        origen = "creada por el corte" if len(hijosDe[padreDe[segId]]) > 1 else "hueso preexistente"
        print(f"  {nombre}: {origen}, vol={volf:.4f}cm3 "
              f"({nRell} cavidad(es) rellenada(s))")

    idsExportar = [segId for segId, _, _ in conservados]
    nombresInteriores = set(n for _, n, esT in conservados if esT)
    for segId in ids:
        if segId not in idsExportar:
            segTrabajo.GetSegmentation().RemoveSegment(segId)

    # ---- 9. Exportar a mallas ----
    segTrabajo.GetSegmentation().SetConversionParameter("Smoothing factor", SUAVIZADO_EXPORT)
    segTrabajo.RemoveClosedSurfaceRepresentation()
    segTrabajo.CreateClosedSurfaceRepresentation()

    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    ok = slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
        segTrabajo, idsExportar, folderItemId)
    if not ok:
        print("Fallo la exportacion a modelos. Revisa el Error log de Slicer.")
        slicer.mrmlScene.RemoveNode(editorNode)
        return

    # ---- 10. Control de calidad ----
    print("")
    print("--- CONTROL DE CALIDAD ---")
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(folderItemId, hijos)
    listos, total = 0, 0
    for i in range(hijos.GetNumberOfIds()):
        node = shNode.GetItemDataNode(hijos.GetId(i))
        if node is None or not node.IsA("vtkMRMLModelNode"):
            continue
        nombre = node.GetName()
        if not (nombre.startswith("Tapa_") or nombre.startswith("Resto_")
                or nombre.startswith("Hueso_")):
            continue
        total += 1
        _limpiar_malla(node)
        cerrada, volumen = _control_calidad(node)
        estado = "WATERTIGHT (apto STL)" if cerrada else "ABIERTA (NO apta para STL)"
        tipo = "tapa de hueso" if nombre in nombresInteriores else "resto"
        print(f"  {nombre} [{tipo}]: {estado}  vol={volumen/1000.0:.3f}cm3")
        if nombre in nombresInteriores and node.GetDisplayNode():
            node.GetDisplayNode().SetColor(0.9, 0.25, 0.25)
        if cerrada:
            listos += 1

    # ---- 11. Ocultar lo que estorba ----
    try:
        modeloOriginal = slicer.util.getNode(NOMBRE_MODELO_CRANEO)
        if modeloOriginal.IsA("vtkMRMLModelNode") and modeloOriginal.GetDisplayNode():
            modeloOriginal.GetDisplayNode().SetVisibility(False)
    except Exception:
        pass
    for h in herramientas:
        if h.GetDisplayNode():
            h.GetDisplayNode().SetVisibility(False)
    if segOrigen.GetDisplayNode():
        segOrigen.GetDisplayNode().SetVisibility(False)
    if segTrabajo.GetDisplayNode():
        segTrabajo.GetDisplayNode().SetVisibility(False)

    slicer.mrmlScene.RemoveNode(editorNode)

    print("")
    print(f"Listo: {total} pieza(s), {listos} watertight, "
          f"{len(nombresInteriores)} tapa(s) de hueso (en rojo).")
    print(f"Estan en la carpeta '{NOMBRE_CARPETA_SH}' del panel Data.")
    print("Para exportar STL:  exportar_stl('C:/ruta/de/salida')")


def exportar_stl(carpeta_destino):
    """Exporta a STL todos los fragmentos watertight."""
    import os
    if not os.path.isdir(carpeta_destino):
        os.makedirs(carpeta_destino)
    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(folderItemId, hijos)
    n = 0
    for i in range(hijos.GetNumberOfIds()):
        node = shNode.GetItemDataNode(hijos.GetId(i))
        if node is None or not node.IsA("vtkMRMLModelNode"):
            continue
        cerrada, _ = _control_calidad(node)
        if not cerrada:
            print(f"  {node.GetName()}: NO watertight, no se exporta.")
            continue
        ruta = os.path.join(carpeta_destino, node.GetName() + ".stl")
        slicer.util.saveNode(node, ruta)
        print(f"  {node.GetName()} -> {ruta}")
        n += 1
    print(f"{n} STL exportado(s).")


print("Bloque F v5 (unificado) cargado.")
print(f"  ancho de corte: {ANCHO_CORTE_MM}mm   profundidad: {PROFUNDIDAD_CORTE_MM}mm")
print("Comandos:")
print("  cortar()                    -> genera las cintas y ejecuta el corte completo")
print("  generar_cintas()            -> solo genera las cintas, para inspeccionarlas")
print("  cortar(regenerar=False)     -> corta con las cintas que ya estan en la escena")
print("  exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight")
