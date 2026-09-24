# ============================================================
# CRANIOPLAN - BLOQUE F v10 (UNIFICADO)
#
# Combina lo mejor de las dos ramas:
#   - MOTOR DE CORTE del codigo del colaborador: anillo construido en voxeles
#     barriendo la curva RESAMPLEADA (puntos cada 0.5mm) a lo largo de la
#     normal del hueso SUAVIZADA (promedio movil). Esto da un corte curvo,
#     continuo y sin el 'twist' ni los huecos por aliasing de las versiones
#     con anillos mallados. Es la geometria de corte que se pidio conservar.
#   - CLASIFICACION POR CONECTIVIDAD (no por particion geometrica): tras
#     restar el anillo, se re-etiquetan las islas y se mapea cada hija a su
#     isla madre. Cuando una tapa queda partida por una sutura, cada pedazo
#     sale como un fragmento SEPARADO (esto es lo que la particion geometrica
#     del colaborador fusionaba en uno solo).
#   - IDENTIFICACION INTERIOR/EXTERIOR por el TAPON extruido de la curva:
#     un fragmento es tapa si >50% de sus voxeles caen dentro del prisma del
#     lazo. Funciona con tapas concavas o partidas en varios pedazos.
#
# REGLAS DE CONSERVACION:
#   - Toda pieza INTERIOR al lazo (tapa) se conserva SIN umbral de tamano.
#     Solo se descarta ruido numerico de pocos voxeles (RUIDO_VOXELES).
#   - Piezas EXTERIORES creadas por el corte -> se conservan (Resto_N).
#   - Islas preexistentes que el corte no toco -> Hueso_N (con filtro de
#     volumen minimo para descartar basura de segmentacion).
#
# CURVAS ABIERTAS (v10): ademas de curvas cerradas (lazos que definen una
#   tapa), ahora se admiten curvas ABIERTAS (vtkMRMLMarkupsCurveNode). Una
#   linea abierta no encierra nada: solo separa una placa si la cruza de
#   borde a borde. Si el corte abierto no llega de un borde al otro, NO
#   separa (es geometria, no un bug): se graba la ranura (surco visible) y
#   se avisa por consola. Las piezas que un corte abierto si separa salen
#   como Fragmento_N.
#
# SALIDA: solo modelos (vtkMRMLModelNode), watertight, en la carpeta
#   'CranioPlan_Fragmentos'. Aptos para mover, medir y exportar STL.
#   Las tapas de cada curva cerrada se agrupan en subcarpeta 'Tapas_<curva>'.
#
# FLUJO:
#   Bloque A+B -> confirmar_craneo() -> enviar_a_planner()
#   -> dibujar una o mas curvas cerradas sobre Craneo_Final
#      (Shortest distance on surface + Constrain to Model)
#   -> cargar este archivo -> cortar()
#
# REQUISITO: scipy (ya se usa en el Bloque A).
# ============================================================

import numpy as np
import vtk
from scipy import ndimage


# ==================== CONFIGURACION ====================

# --- Geometria del corte ---
GROSOR_CORTE_MM        = 1.5    # kerf REAL objetivo: ancho de la ranura de corte
PROFUNDIDAD_CORTE_MM   = 4.0    # semi-barrido perpendicular (atraviesa el espesor y basta)
#   OJO: el barrido va de -PROFUNDIDAD a +PROFUNDIDAD. Con hueso de 1-2mm,
#   4.0 (=8mm total) ya sobra. Barridos largos + normal inclinada = corte
#   ancho, porque la linea cruza el hueso en diagonal. No subir sin motivo.
DIST_MUESTREO_MM       = 0.5    # resampleo de la curva: mas chico = mas suave
RADIO_NORMAL_MM        = 3.0    # radio para promediar la normal de la superficie
VENTANA_SUAVIZADO      = 3      # promedio movil de normales: +-3 puntos = 7 puntos

# --- Nodos ---
NOMBRE_MODELO_CRANEO   = "Craneo_Final"       # Model node (para normales del hueso)
NOMBRE_SEGMENTACION    = "Craneo_Automatico"  # nodo de segmentacion del Bloque A
NOMBRE_SEGMENTO_HUESO  = "Craneo_Final"       # segmento dentro de ese nodo
NOMBRE_SEG_TRABAJO     = "CranioPlan_Corte"
NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"

# --- Clasificacion y filtrado ---
RUIDO_VOXELES          = 30     # piezas con menos voxeles que esto = ruido numerico
VOLUMEN_MINIMO_CM3     = 0.10   # solo para islas PREEXISTENTES (no para tapas)
FRACCION_INTERIOR      = 0.5    # fraccion de voxeles dentro del tapon para ser tapa
MAX_CAVIDAD_MM3        = 2000.0 # cavidades internas mayores NO se rellenan
SUAVIZADO_EXPORT       = "0.3"  # suavizado al convertir a malla
UMBRAL_LIMPIEZA_MALLA  = 0.10   # limpieza de islas residuales de marching cubes


# ============================================================
# PARTE 1 - MOTOR DE CORTE (adaptado del codigo del colaborador)
# ============================================================

def _resamplear_puntos(puntosVtk, distanciaMM):
    """
    Recorre la polilinea de la curva y devuelve puntos espaciados
    uniformemente cada distanciaMM. Es lo que da la suavidad: en vez de
    saltar entre los puntos de control, se toma la curva densa de Slicer y
    se la muestrea fino y parejo.
    """
    n = puntosVtk.GetNumberOfPoints()
    if n < 2:
        return []
    original = [np.array(puntosVtk.GetPoint(i)) for i in range(n)]

    resampleados = [original[0]]
    acumulada = 0.0
    anterior = original[0]
    for i in range(1, n):
        actual = original[i]
        segmento = actual - anterior
        largo = np.linalg.norm(segmento)
        while acumulada + largo >= distanciaMM:
            falta = distanciaMM - acumulada
            direccion = segmento / largo if largo > 1e-9 else segmento
            nuevo = anterior + direccion * falta
            resampleados.append(nuevo)
            anterior = nuevo
            segmento = actual - anterior
            largo = np.linalg.norm(segmento)
            acumulada = 0.0
        acumulada += largo
        anterior = actual

    if np.linalg.norm(resampleados[-1] - original[-1]) > 1e-6:
        resampleados.append(original[-1])
    return resampleados


def _quitar_coincidentes(posiciones, esCerrada, toleranciaMM=0.2):
    """Elimina puntos consecutivos casi superpuestos (un tramo de largo cero
    rompe el calculo de la normal y arruina el corte sin avisar)."""
    limpias = []
    for p in posiciones:
        if not limpias or np.linalg.norm(p - limpias[-1]) > toleranciaMM:
            limpias.append(p)
    if esCerrada and len(limpias) > 2:
        if np.linalg.norm(limpias[-1] - limpias[0]) <= toleranciaMM:
            limpias.pop()
    return limpias


def _normal_promedio_en_punto(punto, locator, normalesArray, radioMM):
    """Normal de superficie promediada sobre los vertices dentro de un radio
    (mas estable que la de una sola celda sobre mallas de marching cubes)."""
    ids = vtk.vtkIdList()
    locator.FindPointsWithinRadius(radioMM, punto.tolist(), ids)
    if ids.GetNumberOfIds() == 0:
        idc = locator.FindClosestPoint(punto.tolist())
        if idc < 0:
            return None
        ids.InsertNextId(idc)
    acum = np.zeros(3)
    for j in range(ids.GetNumberOfIds()):
        acum += np.array(normalesArray.GetTuple3(ids.GetId(j)))
    norma = np.linalg.norm(acum)
    if norma < 1e-6:
        idc = locator.FindClosestPoint(punto.tolist())
        if idc < 0:
            return None
        nrm = np.array(normalesArray.GetTuple3(idc))
        nn = np.linalg.norm(nrm)
        return nrm / nn if nn > 1e-9 else None
    return acum / norma


def _anillo_de_corte(boneMask, curvaNode, mallaHueso, volumeNode, esCerrada):
    """
    Construye el anillo de corte DIRECTAMENTE en voxeles: por cada punto de
    la curva resampleada se barre a lo largo de la normal del hueso
    (suavizada) para atravesar todo el espesor, y luego se dilata el grosor
    de hoja. Continuo por construccion: inmune al twist y a los huecos por
    aliasing. Devuelve la mascara booleana de hueso a quitar, o None.
    """
    puntosCurva = curvaNode.GetCurvePointsWorld()
    if puntosCurva is None or puntosCurva.GetNumberOfPoints() < 2:
        return None

    posiciones = _resamplear_puntos(puntosCurva, DIST_MUESTREO_MM)
    posiciones = _quitar_coincidentes(posiciones, esCerrada)
    posiciones = list(posiciones)
    if esCerrada and len(posiciones) >= 3:
        posiciones.append(posiciones[0])
    n = len(posiciones)
    if n < 2:
        return None

    # Normales del hueso
    nf = vtk.vtkPolyDataNormals()
    nf.SetInputData(mallaHueso)
    nf.ComputePointNormalsOn()
    nf.ComputeCellNormalsOff()
    nf.SplittingOff()
    nf.ConsistencyOn()
    nf.AutoOrientNormalsOn()
    nf.Update()
    mallaN = nf.GetOutput()
    normalesArray = mallaN.GetPointData().GetNormals()
    if normalesArray is None:
        return None
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(mallaN)
    locator.BuildLocator()

    # Normal por punto
    normales = [_normal_promedio_en_punto(p, locator, normalesArray, RADIO_NORMAL_MM)
                for p in posiciones]
    validos = [i for i, nr in enumerate(normales) if nr is not None]
    if not validos:
        return None
    for i in range(n):
        if normales[i] is None:
            j = min(validos, key=lambda k: abs(k - i))
            normales[i] = np.array(normales[j])

    # Alinear al sentido medio (evita promediar normales opuestas)
    media = np.zeros(3)
    for i in validos:
        media += normales[i]
    if np.linalg.norm(media) > 1e-9:
        media /= np.linalg.norm(media)
        for i in range(n):
            if float(np.dot(normales[i], media)) < 0:
                normales[i] = -normales[i]

    # Promedio movil a lo largo de la curva (circular si es cerrada)
    normalesSuaves = []
    for i in range(n):
        acum = np.zeros(3)
        for d in range(-VENTANA_SUAVIZADO, VENTANA_SUAVIZADO + 1):
            k = (i + d) % n if esCerrada else min(max(i + d, 0), n - 1)
            acum += normales[k]
        nn = np.linalg.norm(acum)
        normalesSuaves.append(acum / nn if nn > 1e-9 else normales[i])

    # Barrido en voxeles
    rasToIjk = vtk.vtkMatrix4x4()
    volumeNode.GetRASToIJKMatrix(rasToIjk)
    dims = boneMask.shape  # (z, y, x)
    espaciado = volumeNode.GetSpacing()
    pasoMM = max(min(espaciado) * 0.5, 0.1)

    curtain = np.zeros(dims, dtype=bool)

    def marcar(ras):
        ijk = [0.0, 0.0, 0.0, 0.0]
        rasToIjk.MultiplyPoint([float(ras[0]), float(ras[1]), float(ras[2]), 1.0], ijk)
        ii, jj, kk = int(round(ijk[0])), int(round(ijk[1])), int(round(ijk[2]))
        if 0 <= kk < dims[0] and 0 <= jj < dims[1] and 0 <= ii < dims[2]:
            curtain[kk, jj, ii] = True

    for i, p in enumerate(posiciones):
        direccion = normalesSuaves[i]
        t = -PROFUNDIDAD_CORTE_MM
        while t <= PROFUNDIDAD_CORTE_MM:
            marcar(np.asarray(p, dtype=float) + direccion * t)
            t += pasoMM

    if not curtain.any():
        return None

    # Dilatar hasta el kerf objetivo. El barrido ya deja una cortina de ~1 voxel
    # de espesor; la dilatacion la lleva al GROSOR_CORTE_MM pedido. Se usa una
    # estructura de conectividad de CARA (cruz 3D, 6-conexa) en vez de la cubica
    # 26-conexa: la cubica engorda tambien en diagonal y hace el corte mas ancho
    # que el nominal. La cruz crece 1 voxel por iteracion por cara, que es lo que
    # queremos para controlar el ancho.
    #
    # radio en voxeles para alcanzar GROSOR_CORTE_MM de ANCHO TOTAL: como la
    # cortina ya aporta ~1 voxel, dilatamos (kerf - 1 voxel)/2 por lado.
    kerf_vox = GROSOR_CORTE_MM / min(espaciado)
    radioVox = max(0, int(round((kerf_vox - 1.0) / 2.0)))
    if radioVox > 0:
        estructura = ndimage.generate_binary_structure(3, 1)  # 6-conexa (cruz)
        curtainDil = ndimage.binary_dilation(curtain, structure=estructura, iterations=radioVox)
    else:
        curtainDil = curtain

    anillo = np.logical_and(curtainDil, boneMask)
    return anillo if anillo.any() else None


# ============================================================
# PARTE 2 - TAPON DE LA CURVA (identificacion interior/exterior)
# ============================================================

def _tapon_de_curva(curvaNode, longitud):
    """
    Superficie minima que tapa la curva cerrada, extruida en prisma a lo
    largo de su normal media. Sirve para reconocer que fragmento es interior
    del lazo, aun si es concavo o esta partido por una sutura.
    """
    cap = vtk.vtkPolyData()
    try:
        slicer.modules.markups.logic().GetClosedCurveSurfaceArea(curvaNode, cap)
    except Exception as e:
        print(f"  no pude generar el tapon de '{curvaNode.GetName()}': {e}")
        return None
    if cap.GetNumberOfPoints() < 3:
        return None

    pts = np.array([cap.GetPoint(i) for i in range(cap.GetNumberOfPoints())])
    centro = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - centro)
    normal = vt[2] / np.linalg.norm(vt[2])

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
    """{nombreCurva: mascara booleana} rasterizando el tapon de cada curva."""
    curvas = slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode')
    if not curvas:
        return {}

    previo = slicer.mrmlScene.GetFirstNodeByName("CranioPlan_Tapones")
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)

    segTap = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', "CranioPlan_Tapones")
    segTap.CreateDefaultDisplayNodes()
    segTap.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    if segTap.GetDisplayNode():
        segTap.GetDisplayNode().SetVisibility(False)

    temporales, mapaNombre = [], {}
    for curva in curvas:
        poly = _tapon_de_curva(curva, longitud)
        if poly is None:
            continue
        mNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", f"_Tapon_{curva.GetName()}")
        mNode.SetAndObservePolyData(poly)
        mNode.CreateDefaultDisplayNodes()
        mNode.GetDisplayNode().SetVisibility(False)
        temporales.append(mNode)

        antes = set(_ids_de_segmentos(segTap))
        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(mNode, segTap)
        for i in [x for x in _ids_de_segmentos(segTap) if x not in antes]:
            mapaNombre[i] = curva.GetName()

    mascaras = {}
    for segId, nombreCurva in mapaNombre.items():
        try:
            mascaras[nombreCurva] = slicer.util.arrayFromSegmentBinaryLabelmap(
                segTap, segId, volumeNode).astype(bool)
        except Exception as e:
            print(f"  no pude rasterizar el tapon de '{nombreCurva}': {e}")

    for m in temporales:
        slicer.mrmlScene.RemoveNode(m)
    slicer.mrmlScene.RemoveNode(segTap)
    return mascaras


# ============================================================
# PARTE 3 - UTILIDADES DE SEGMENTACION / MALLA
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


def _rellenar_cavidades(arr, voxel_mm3, max_mm3):
    """Rellena solo cavidades internas por debajo de max_mm3 (el diploe),
    sin tocar la cavidad craneal ni el exterior."""
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


def _subcarpeta(shNode, parentItemId, nombre):
    """Devuelve (o crea) una subcarpeta dentro de parentItemId."""
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(parentItemId, hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if (shNode.GetItemName(itemId) == nombre
                and shNode.GetItemDataNode(itemId) is None):
            return itemId
    return shNode.CreateFolderItem(parentItemId, nombre)


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

def cortar():
    # ---- 1. Entradas ----
    try:
        segOrigen = slicer.util.getNode(NOMBRE_SEGMENTACION)
    except Exception:
        print(f"ERROR: no encuentro el nodo de segmentacion '{NOMBRE_SEGMENTACION}'.")
        return
    try:
        craneoModel = slicer.util.getNode(NOMBRE_MODELO_CRANEO)
    except Exception:
        print(f"ERROR: no encuentro el modelo '{NOMBRE_MODELO_CRANEO}'. Corre enviar_a_planner().")
        return
    if not craneoModel.IsA("vtkMRMLModelNode"):
        print(f"ERROR: '{NOMBRE_MODELO_CRANEO}' no es un Model node.")
        return

    volumeNode = _volumen_de_referencia(segOrigen)
    if volumeNode is None:
        print("ERROR: no encuentro el volumen de referencia del CT.")
        return

    segIdHueso = segOrigen.GetSegmentation().GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_HUESO)
    if not segIdHueso:
        print(f"ERROR: no hay segmento '{NOMBRE_SEGMENTO_HUESO}'. Corre confirmar_craneo() primero.")
        return

    curvasCerradas = slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode')
    curvasAbiertas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                      if not n.IsA('vtkMRMLMarkupsClosedCurveNode')]
    curvas = list(curvasCerradas) + list(curvasAbiertas)
    if not curvas:
        print("ERROR: no hay ninguna curva de corte en la escena. Dibuja al menos una.")
        return

    spacing = volumeNode.GetSpacing()
    voxel_mm3 = spacing[0] * spacing[1] * spacing[2]
    piso_kerf = 2.0 * max(spacing)
    print(f"Volumen de referencia: {volumeNode.GetName()}  spacing={spacing}")
    print(f"Curvas cerradas (lazo/tapa): {[c.GetName() for c in curvasCerradas]}")
    print(f"Curvas abiertas (linea):     {[c.GetName() for c in curvasAbiertas]}")
    if GROSOR_CORTE_MM < piso_kerf:
        print(f"ADVERTENCIA: GROSOR_CORTE_MM={GROSOR_CORTE_MM} < piso de {piso_kerf:.2f}mm (2 voxeles).")

    mallaHueso = craneoModel.GetPolyData()

    # ---- 2. Copia de trabajo (no tocamos el nodo del Bloque A) ----
    previo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_SEG_TRABAJO)
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)

    segTrabajo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', NOMBRE_SEG_TRABAJO)
    segTrabajo.CreateDefaultDisplayNodes()
    segTrabajo.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segTrabajo.GetSegmentation().CopySegmentFromSegmentation(
        segOrigen.GetSegmentation(), segIdHueso)
    idHueso = _ids_de_segmentos(segTrabajo)[0]

    # ---- 3. Mascara del hueso e islas madre ANTES de cortar ----
    boneMask = slicer.util.arrayFromSegmentBinaryLabelmap(
        segTrabajo, idHueso, volumeNode).astype(bool)
    etiquetasPadre, nPadres = ndimage.label(boneMask)
    print(f"Islas de hueso antes del corte: {nPadres}")

    # ---- 4. Tapones de cada curva ----
    taponesMasc = _mascaras_de_tapones(volumeNode, 2.0 * PROFUNDIDAD_CORTE_MM)
    print(f"Tapones generados: {list(taponesMasc.keys()) if taponesMasc else 'ninguno'}")

    # ---- 5. Construir el anillo de corte de cada curva y unirlos ----
    # Se guarda el nombre de cada curva ABIERTA para despues avisar si su
    # corte no llego a separar nada (no cruza de borde a borde).
    nombresAbiertas = set(c.GetName() for c in curvasAbiertas)
    anilloTotal = np.zeros_like(boneMask)
    for curva in curvas:
        esCerrada = curva.IsA('vtkMRMLMarkupsClosedCurveNode')
        anillo = _anillo_de_corte(boneMask, curva, mallaHueso, volumeNode, esCerrada=esCerrada)
        if anillo is None:
            print(f"  ADVERTENCIA: no se pudo construir el anillo de '{curva.GetName()}'.")
            continue
        anilloTotal |= anillo
        tipo = "lazo" if esCerrada else "linea"
        print(f"  anillo de '{curva.GetName()}' ({tipo}): {int(np.count_nonzero(anillo))} voxeles")

    if not anilloTotal.any():
        print("ERROR: no se pudo construir ningun anillo de corte.")
        slicer.mrmlScene.RemoveNode(segTrabajo)
        return

    # ---- 6. Restar el anillo: hueso AND NOT anillo ----
    boneCut = np.logical_and(boneMask, np.logical_not(anilloTotal))
    print(f"Anillo total restado: {int(np.count_nonzero(anilloTotal))} voxeles de hueso.")

    # ---- 7. Separar en islas por CONECTIVIDAD ----
    etiquetasHijas, nHijas = ndimage.label(boneCut)
    print(f"Islas tras el corte: {nHijas}")

    # ---- 8. Clasificar cada hija ----
    # Mapa hija -> madre (por etiqueta mayoritaria), y hijas por madre
    hijasDeMadre = {}
    for lab in range(1, nHijas + 1):
        m = (etiquetasHijas == lab)
        nvox = int(np.count_nonzero(m))
        if nvox < RUIDO_VOXELES:
            continue
        etiquetas = etiquetasPadre[m]
        etiquetas = etiquetas[etiquetas > 0]
        madre = int(np.bincount(etiquetas).argmax()) if etiquetas.size else 0
        hijasDeMadre.setdefault(madre, []).append(lab)

    # Fraccion de cada hija dentro de cada tapon -> a que curva es interior
    def curva_interior_de(mask):
        nvox = np.count_nonzero(mask)
        mejor, mejorFrac = None, 0.0
        for nombreCurva, tap in taponesMasc.items():
            if tap.shape != mask.shape:
                continue
            frac = np.count_nonzero(mask & tap) / float(nvox) if nvox else 0.0
            if frac > mejorFrac:
                mejor, mejorFrac = nombreCurva, frac
        return mejor if mejorFrac >= FRACCION_INTERIOR else None

    # Precalcular, por madre, si alguna hija es interior de un tapon (= hubo
    # un lazo cerrado sobre esa madre). Sirve para nombrar el exterior:
    #   - exterior de una madre CON tapa  -> Resto_N (resto del lazo)
    #   - pieza de una madre SIN tapa que se partio -> Fragmento_N (corte abierto)
    madreTieneTapa = {}
    for madre, labs in hijasDeMadre.items():
        madreTieneTapa[madre] = any(
            curva_interior_de(etiquetasHijas == lab) is not None for lab in labs)

    # Construir la lista de piezas a conservar con su nombre
    piezas = []  # (mask, nombre, esTapa, curva, vol)
    contadorTapa = {}   # por curva
    contadorResto = 0
    contadorFragmento = 0
    contadorHueso = 0

    for madre, labs in hijasDeMadre.items():
        creadaPorCorte = len(labs) > 1
        for lab in labs:
            mask = (etiquetasHijas == lab)
            nvox = int(np.count_nonzero(mask))
            vol_cm3 = nvox * voxel_mm3 / 1000.0
            curva = curva_interior_de(mask)

            if curva is not None:
                # TAPA (interior de un lazo cerrado): se conserva SIEMPRE
                contadorTapa.setdefault(curva, 0)
                contadorTapa[curva] += 1
                piezas.append((mask, f"Tapa_{curva}", True, curva, vol_cm3))
            elif creadaPorCorte and madreTieneTapa[madre]:
                # Exterior de un lazo cerrado
                contadorResto += 1
                piezas.append((mask, f"Resto_{contadorResto}", False, None, vol_cm3))
            elif creadaPorCorte:
                # Pieza de un corte ABIERTO que separo la placa
                contadorFragmento += 1
                piezas.append((mask, f"Fragmento_{contadorFragmento}", False, None, vol_cm3))
            else:
                # Isla preexistente que el corte no separo: filtro de volumen minimo
                if vol_cm3 < VOLUMEN_MINIMO_CM3:
                    print(f"  descartado: isla preexistente diminuta, vol={vol_cm3:.4f}cm3")
                    continue
                contadorHueso += 1
                piezas.append((mask, f"Hueso_{contadorHueso}", False, None, vol_cm3))

    if not piezas:
        print("ERROR: no sobrevivio ningun fragmento.")
        slicer.mrmlScene.RemoveNode(segTrabajo)
        return

    # Aviso: cortes abiertos que NO separaron (no cruzan de borde a borde).
    # Heuristica: cada corte abierto que separa produce >=1 Fragmento_. Si hay
    # menos Fragmentos que curvas abiertas, al menos una no separo.
    # (Con un solo corte abierto por corrida, que es el caso normal, es exacto.)
    if nombresAbiertas and contadorFragmento < len(nombresAbiertas):
        sinSeparar = len(nombresAbiertas) - contadorFragmento
        print("")
        print(f"AVISO: {sinSeparar} corte(s) abierto(s) NO separaron el hueso.")
        print("Un corte abierto solo divide si cruza la placa de borde a borde.")
        print("La ranura quedo grabada en el hueso, pero la pieza sigue entera.")
        print("Para dividir, extende la linea hasta el otro borde de la placa.")
        print("")

    # Sufijos a/b/c cuando una misma curva tiene varias tapas (sutura de por medio)
    sufijos = {c: 0 for c, n in contadorTapa.items() if n > 1}
    piezasFinales = []  # (mask, nombre, esTapa, curva, vol)
    for mask, nombre, esTapa, curva, vol in piezas:
        if esTapa and curva in sufijos:
            letra = chr(ord('a') + sufijos[curva])
            sufijos[curva] += 1
            nombre = f"{nombre}_{letra}"
        piezasFinales.append((mask, nombre, esTapa, curva, vol))

    # ---- 9. Materializar cada pieza como segmento, rellenar cavidades ----
    idsPorNombre = {}
    curvaPorNombre = {}   # nombre de pieza -> curva de la que es tapa (o None)
    for mask, nombre, esTapa, curva, vol in piezasFinales:
        arr, nRell = _rellenar_cavidades(mask.copy(), voxel_mm3, MAX_CAVIDAD_MM3)
        sid = segTrabajo.GetSegmentation().AddEmptySegment(nombre)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            arr.astype(np.uint8), segTrabajo, sid, volumeNode)
        idsPorNombre[nombre] = sid
        curvaPorNombre[nombre] = curva if esTapa else None
        volf = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        tipo = "TAPA" if esTapa else "resto/hueso"
        print(f"  {nombre} [{tipo}]: vol={volf:.4f}cm3 ({nRell} cavidad(es) rellenada(s))")

    # Quitar el segmento de hueso original (ya lo reemplazamos por las piezas)
    segTrabajo.GetSegmentation().RemoveSegment(idHueso)

    nombresTapas = set(n for _, n, esT, _, _ in piezasFinales if esT)

    # ---- 10. Exportar a modelos ----
    segTrabajo.GetSegmentation().SetConversionParameter("Smoothing factor", SUAVIZADO_EXPORT)
    segTrabajo.RemoveClosedSurfaceRepresentation()
    segTrabajo.CreateClosedSurfaceRepresentation()

    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    ok = slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
        segTrabajo, list(idsPorNombre.values()), folderItemId)
    if not ok:
        print("Fallo la exportacion a modelos. Revisa el Error log de Slicer.")
        return

    # ---- 11. Control de calidad + resolver el nodo de cada pieza por nombre ----
    # El export nombra cada modelo igual que su segmento. Buscamos el modelo
    # por nombre directamente, en vez de recorrer el arbol (que era fragil).
    print("")
    print("--- CONTROL DE CALIDAD ---")
    nodoPorNombre = {}
    modelos = slicer.util.getNodesByClass("vtkMRMLModelNode")
    listos, total = 0, 0
    for nombre in idsPorNombre.keys():
        # el modelo exportado puede tener el nombre exacto; si hay colision,
        # Slicer le agrega sufijo, asi que buscamos por prefijo exacto.
        candidatos = [m for m in modelos if m.GetName() == nombre
                      and shNode.GetItemParent(shNode.GetItemByDataNode(m)) == folderItemId]
        if not candidatos:
            candidatos = [m for m in modelos if m.GetName() == nombre]
        if not candidatos:
            print(f"  ADVERTENCIA: no encontre el modelo '{nombre}' exportado.")
            continue
        node = candidatos[0]
        nodoPorNombre[nombre] = node
        total += 1
        _limpiar_malla(node)
        cerrada, volumen = _control_calidad(node)
        estado = "WATERTIGHT (apto STL)" if cerrada else "ABIERTA (NO apta para STL)"
        tipo = "tapa de hueso" if nombre in nombresTapas else "resto/hueso"
        print(f"  {nombre} [{tipo}]: {estado}  vol={volumen/1000.0:.4f}cm3")
        if nombre in nombresTapas and node.GetDisplayNode():
            node.GetDisplayNode().SetColor(0.9, 0.25, 0.25)
        if cerrada:
            listos += 1

    # ---- 11b. Agrupar las tapas de cada curva en una subcarpeta ----
    # Usa el nodo directo (nodoPorNombre) y la curva conocida (curvaPorNombre),
    # sin re-parsear strings ni recorrer el arbol.
    itemsPorCurva = {}
    for nombre, node in nodoPorNombre.items():
        curva = curvaPorNombre.get(nombre)
        if curva is None:
            continue  # solo se agrupan las tapas
        itemId = shNode.GetItemByDataNode(node)
        if itemId:
            itemsPorCurva.setdefault(curva, []).append(itemId)

    for curva, items in itemsPorCurva.items():
        sub = _subcarpeta(shNode, folderItemId, f"Tapas_{curva}")
        for itemId in items:
            shNode.SetItemParent(itemId, sub)
    if itemsPorCurva:
        print(f"Tapas agrupadas en subcarpeta(s): "
              f"{['Tapas_' + c for c in itemsPorCurva.keys()]}")

    # ---- 12. Ocultar lo que estorba ----
    if craneoModel.GetDisplayNode():
        craneoModel.GetDisplayNode().SetVisibility(False)
    if segOrigen.GetDisplayNode():
        segOrigen.GetDisplayNode().SetVisibility(False)
    if segTrabajo.GetDisplayNode():
        segTrabajo.GetDisplayNode().SetVisibility(False)

    # Colapsar el nodo de trabajo en el arbol para que no moleste.
    shTrabajo = shNode.GetItemByDataNode(segTrabajo)
    if shTrabajo:
        shNode.SetItemExpanded(shTrabajo, False)

    print("")
    print(f"Listo: {total} pieza(s), {listos} watertight, "
          f"{len(nombresTapas)} tapa(s) de hueso (en rojo).")
    print(f"Estan en la carpeta '{NOMBRE_CARPETA_SH}' del panel Data.")
    print(f"El nodo de trabajo '{NOMBRE_SEG_TRABAJO}' quedo oculto; podes borrarlo con limpiar().")
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


print("Bloque F v10 (unificado) cargado.")
print(f"  grosor de corte: {GROSOR_CORTE_MM}mm   profundidad: {PROFUNDIDAD_CORTE_MM}mm   "
      f"muestreo: {DIST_MUESTREO_MM}mm")
print("Comandos:")
print("  cortar()                    -> corte curvo + separacion por conectividad")
print("  exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight")
print("  limpiar()                   -> borra el nodo de trabajo CranioPlan_Corte")
