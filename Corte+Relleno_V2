# ============================================================
# CRANIOPLAN - BLOQUE F v3 (UNIFICADO)
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
VOLUMEN_MINIMO_CM3     = 0.5    # fragmentos mas chicos se descartan (ruido)
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
        node = _construir_cinta(curva, craneoConNormales, locator,
                                normalesArray, PREFIJO_HERRAMIENTAS)
        if node is not None:
            generadas.append(node)
    return generadas


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

    # ---- 7. Filtrar por volumen y rellenar cavidades ----
    supervivientes = []
    for segId in ids:
        arr = slicer.util.arrayFromSegmentBinaryLabelmap(
            segTrabajo, segId, volumeNode).astype(bool)
        vol_cm3 = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        if vol_cm3 < VOLUMEN_MINIMO_CM3:
            print(f"  descartado (vol={vol_cm3:.2f}cm3 < {VOLUMEN_MINIMO_CM3}cm3)")
            continue

        arr, nRellenadas = _rellenar_cavidades(arr, voxel_mm3, MAX_CAVIDAD_MM3)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            arr.astype(np.uint8), segTrabajo, segId, volumeNode)

        vol_final = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        supervivientes.append((segId, vol_final))
        print(f"  fragmento vol={vol_final:.2f}cm3  ({nRellenadas} cavidad(es) rellenada(s))")

    idsOk = [s for s, _ in supervivientes]
    for segId in ids:
        if segId not in idsOk:
            segTrabajo.GetSegmentation().RemoveSegment(segId)

    if not supervivientes:
        print("ERROR: no sobrevivio ningun fragmento. Baja VOLUMEN_MINIMO_CM3.")
        slicer.mrmlScene.RemoveNode(editorNode)
        return

    # ---- 8. Nombrar por volumen descendente ----
    supervivientes.sort(key=lambda t: t[1], reverse=True)
    for i, (segId, vol) in enumerate(supervivientes, start=1):
        segTrabajo.GetSegmentation().GetSegment(segId).SetName(f"Fragmento_{i}")

    # ---- 9. Exportar a mallas ----
    segTrabajo.GetSegmentation().SetConversionParameter("Smoothing factor", SUAVIZADO_EXPORT)
    segTrabajo.RemoveClosedSurfaceRepresentation()
    segTrabajo.CreateClosedSurfaceRepresentation()

    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    ok = slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
        segTrabajo, [s for s, _ in supervivientes], folderItemId)
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
        if not node.GetName().startswith("Fragmento_"):
            continue
        total += 1
        _limpiar_malla(node)
        cerrada, volumen = _control_calidad(node)
        estado = "WATERTIGHT (apto STL)" if cerrada else "ABIERTA (NO apta para STL)"
        print(f"  {node.GetName()}: {estado}  vol={volumen/1000.0:.2f}cm3")
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
    print(f"Listo: {total} fragmento(s), {listos} watertight.")
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


print("Bloque F v3 (unificado) cargado.")
print(f"  ancho de corte: {ANCHO_CORTE_MM}mm   profundidad: {PROFUNDIDAD_CORTE_MM}mm")
print("Comandos:")
print("  cortar()                    -> genera las cintas y ejecuta el corte completo")
print("  generar_cintas()            -> solo genera las cintas, para inspeccionarlas")
print("  cortar(regenerar=False)     -> corta con las cintas que ya estan en la escena")
print("  exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight")
