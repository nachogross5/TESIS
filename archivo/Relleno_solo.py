# ============================================================
# BLOQUE F v2 - CORTE DE OSTEOTOMIA EN DOMINIO DE VOXELS
#
# Reemplaza el boolean de mallas (vtkBooleanOperationPolyDataFilter), que
# generaba "empapelado": superficies abiertas sueltas en vez de solidos, y
# una misma pieza repartida en 3 fragmentos (tabla externa / interna / pared).
#
# LOGICA:
#   1. Toma el segmento 'Craneo_Final' del Bloque A (voxels, no la malla).
#   2. Rasteriza las cintas 'HerramientaCorte*' generadas por tu script actual.
#   3. Resta voxel a voxel: hueso AND NOT cortes. No hay topologia que romper.
#   4. Split islands -> cada pieza es UN componente conectado, por definicion.
#   5. Rellena las cavidades del diploe con limite de tamano (macizo relleno).
#   6. Exporta a modelos watertight, listos para mover, medir y exportar STL.
#
# NO TOCA el script de la cinta. Ese sigue igual: curva CC -> HerramientaCorte.
# NO TOCA el nodo del Bloque A. Trabaja sobre una copia.
#
# ORDEN DE USO:
#   Bloque A+B v5 -> confirmar_craneo() -> enviar_a_planner()
#   -> curva CC sobre Craneo_Final -> script de la cinta (HerramientaCorte)
#   -> [repetir curva + cinta por cada corte, se acumulan]
#   -> este archivo -> cortar()
# ============================================================

import numpy as np
import vtk
from scipy import ndimage


# ==================== CONFIGURACION ====================

NOMBRE_SEGMENTACION   = "Craneo_Automatico"   # nodo de segmentacion del Bloque A
NOMBRE_SEGMENTO_HUESO = "Craneo_Final"        # segmento dentro de ese nodo
PREFIJO_HERRAMIENTAS  = "HerramientaCorte"    # toma TODOS los modelos que empiecen asi
NOMBRE_SEG_TRABAJO    = "CranioPlan_Corte"    # nodo de trabajo (se recrea cada corrida)
NOMBRE_CARPETA_SH     = "CranioPlan_Fragmentos"

VOLUMEN_MINIMO_CM3    = 0.5    # fragmentos mas chicos que esto se descartan (ruido)
MAX_CAVIDAD_MM3       = 2000.0 # cavidades internas mayores NO se rellenan (protege la cavidad craneal)
SUAVIZADO_EXPORT      = "0.3"  # factor de suavizado al convertir a malla
UMBRAL_LIMPIEZA_MALLA = 0.10   # mismo criterio que tu _limpiar_islas_residuales_malla del v5


# ==================== UTILIDADES ====================

def _volumen_de_referencia(segNode):
    """Recupera el volumen que define la geometria del labelmap."""
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
    """Devuelve (widget, editorNode). Acordate de borrar el editorNode al final."""
    w = slicer.qMRMLSegmentEditorWidget()
    w.setMRMLScene(slicer.mrmlScene)
    eNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
    w.setMRMLSegmentEditorNode(eNode)
    w.setSegmentationNode(segNode)
    w.setSourceVolumeNode(volumeNode)
    return w, eNode


def _rellenar_cavidades(arr, voxel_mm3, max_mm3):
    """
    Rellena SOLO las cavidades internas por debajo de max_mm3.
    Las cavidades que tocan el borde del array (= exterior) y las grandes
    (= cavidad craneal) quedan intactas. Esto da el 'macizo relleno' sin
    convertir el craneo en una bocha maciza.
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

    a_rellenar = []
    for i in range(1, n + 1):
        if i in ids_borde:
            continue
        if conteos[i] * voxel_mm3 <= max_mm3:
            a_rellenar.append(i)

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
    """Igual que tu _limpiar_islas_residuales_malla del Bloque A v5."""
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
    """Devuelve (esWatertight, volumen_mm3). Watertight = apto para STL."""
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


# ==================== PROCESO PRINCIPAL ====================

def cortar():
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

    herramientas = [n for n in slicer.util.getNodesByClass('vtkMRMLModelNode')
                    if n.GetName().startswith(PREFIJO_HERRAMIENTAS)]
    if not herramientas:
        print(f"ERROR: no hay ningun modelo que empiece con '{PREFIJO_HERRAMIENTAS}'.")
        print("Genera al menos una cinta con el script de corte antes de correr esto.")
        return

    spacing = volumeNode.GetSpacing()
    voxel_mm3 = spacing[0] * spacing[1] * spacing[2]
    print(f"Volumen de referencia: {volumeNode.GetName()}  spacing={spacing}")
    print(f"Herramientas de corte encontradas: {[h.GetName() for h in herramientas]}")

    # ---- 2. Nodo de trabajo (no tocamos el del Bloque A) ----
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
    print("Copia de trabajo creada.")

    # ---- 3. Rasterizar las cintas e unirlas en un solo segmento 'Cortes' ----
    idsCortes = []
    for h in herramientas:
        antes = set(_ids_de_segmentos(segTrabajo))
        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(h, segTrabajo)
        nuevos = [i for i in _ids_de_segmentos(segTrabajo) if i not in antes]
        idsCortes.extend(nuevos)
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

    # ---- 4. La resta: hueso AND NOT cortes ----
    widget.setCurrentSegmentID(idHueso)
    widget.setActiveEffectByName("Logical operators")
    efecto = widget.activeEffect()
    efecto.setParameter("Operation", "SUBTRACT")
    efecto.setParameter("ModifierSegmentID", idCortes)
    efecto.self().onApply()
    segTrabajo.GetSegmentation().RemoveSegment(idCortes)
    print("Resta aplicada.")

    # ---- 5. Separar en fragmentos ----
    widget.setCurrentSegmentID(idHueso)
    widget.setActiveEffectByName("Islands")
    efecto = widget.activeEffect()
    efecto.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
    efecto.self().onApply()
    widget.setActiveEffectByName(None)

    ids = _ids_de_segmentos(segTrabajo)
    print(f"Islas generadas por el corte: {len(ids)}")

    # ---- 6. Filtrar por volumen y rellenar cavidades del diploe ----
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

    for segId in ids:
        if segId not in [s for s, _ in supervivientes]:
            segTrabajo.GetSegmentation().RemoveSegment(segId)

    if not supervivientes:
        print("ERROR: no sobrevivio ningun fragmento. Baja VOLUMEN_MINIMO_CM3.")
        slicer.mrmlScene.RemoveNode(editorNode)
        return

    # ---- 7. Nombrar por volumen descendente ----
    supervivientes.sort(key=lambda t: t[1], reverse=True)
    for i, (segId, vol) in enumerate(supervivientes, start=1):
        segTrabajo.GetSegmentation().GetSegment(segId).SetName(f"Fragmento_{i}")

    # ---- 8. Exportar a mallas ----
    segTrabajo.GetSegmentation().SetConversionParameter("Smoothing factor", SUAVIZADO_EXPORT)
    segTrabajo.RemoveClosedSurfaceRepresentation()
    segTrabajo.CreateClosedSurfaceRepresentation()

    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    idsExportar = [s for s, _ in supervivientes]
    ok = slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
        segTrabajo, idsExportar, folderItemId)
    if not ok:
        print("Fallo la exportacion a modelos. Revisa el Error log de Slicer.")
        slicer.mrmlScene.RemoveNode(editorNode)
        return

    # ---- 9. Control de calidad ----
    print("")
    print("--- CONTROL DE CALIDAD ---")
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(folderItemId, hijos)
    listos_para_stl = 0
    total = 0
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
            listos_para_stl += 1

    # ---- 10. Ocultar el craneo original para poder ver los fragmentos ----
    try:
        modeloOriginal = slicer.util.getNode(NOMBRE_SEGMENTO_HUESO)
        if modeloOriginal.IsA("vtkMRMLModelNode") and modeloOriginal.GetDisplayNode():
            modeloOriginal.GetDisplayNode().SetVisibility(False)
    except Exception:
        pass
    if segOrigen.GetDisplayNode():
        segOrigen.GetDisplayNode().SetVisibility(False)
    if segTrabajo.GetDisplayNode():
        segTrabajo.GetDisplayNode().SetVisibility(False)

    slicer.mrmlScene.RemoveNode(editorNode)

    print("")
    print(f"Listo: {total} fragmento(s), {listos_para_stl} watertight.")
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


print("Bloque F v2 cargado.")
print("Comandos:")
print("  cortar()                    -> ejecuta el corte en voxels y genera los fragmentos")
print("  exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight a STL")
