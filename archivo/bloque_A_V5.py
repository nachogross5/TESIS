# ============================================================
# BLOQUE A + B integrados + PUENTE hacia Osteotomy Planner  (v5)
# Parte A: threshold + islands + filtro de bordes + filtro de PROXIMIDAD FISICA
# Parte B: revision manual antes de fusionar (resaltar / eliminar / confirmar)
# Parte PUENTE: exporta Craneo_Final como Model node dentro de una carpeta
#               de Subject Hierarchy, lista para el combo del Osteotomy Planner
#
# CAMBIO v5: se agrega limpieza de islas residuales A NIVEL DE MALLA, dentro de
# enviar_a_planner(). El filtro de islas de v4 (UMBRAL_RELATIVO) opera sobre
# los segmentos (voxels) ANTES de confirmar_craneo(), y funciona bien para eso.
# Pero la conversion de segmento a superficie (ExportSegmentsToModels, marching
# cubes + suavizado) puede generar fragmentos de malla nuevos y chicos que no
# existian como problema a nivel de voxel. Confirmado empiricamente: un craneo
# con 2 islas reales (validadas a nivel segmento) genero 272 regiones conectadas
# en la malla final exportada. Se agrega un segundo filtro, por tamano relativo,
# aplicado sobre la malla ya exportada.
# ============================================================

import numpy as np
from scipy import ndimage

UMBRAL_RELATIVO = 0.05
MARGEN_CENTRALIDAD_MM = 40.0   # distancia fisica de proximidad, no a un centroide
UMBRAL_LIMPIEZA_MALLA = 0.10   # NUEVO v5: regiones de malla con menos de este % del tamano de la mayor se descartan

volumeNodes = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')

if len(volumeNodes) == 0:
    print("No hay ningun volumen cargado. Carga el DICOM primero.")
else:
    volumeNode = volumeNodes[-1]
    print(f"Usando el volumen: {volumeNode.GetName()}")

    segmentationNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode')
    segmentationNode.SetName("Craneo_Automatico")
    segmentationNode.CreateDefaultDisplayNodes()
    segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentId = segmentationNode.GetSegmentation().AddEmptySegment("Hueso")

    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
    segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setSourceVolumeNode(volumeNode)
    segmentEditorWidget.setCurrentSegmentID(segmentId)

    # --- Threshold ---
    segmentEditorWidget.setActiveEffectByName("Threshold")
    thresholdEffect = segmentEditorWidget.activeEffect()
    thresholdEffect.setParameter("MinimumThreshold", "300")
    thresholdEffect.setParameter("MaximumThreshold", "3000")
    thresholdEffect.self().onApply()
    print("Threshold aplicado (300-3000 HU).")

    # --- Islands ---
    segmentEditorWidget.setActiveEffectByName("Islands")
    islandsEffect = segmentEditorWidget.activeEffect()
    islandsEffect.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
    islandsEffect.self().onApply()
    print("Islas separadas en segmentos individuales.")

    segmentacion = segmentationNode.GetSegmentation()
    nSegments = segmentacion.GetNumberOfSegments()
    print(f"Se encontraron {nSegments} islas.")

    spacing = volumeNode.GetSpacing()  # (x,y,z)
    sampling_zyx = (spacing[2], spacing[1], spacing[0])
    voxel_vol_cm3 = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    def toca_borde(maskArr):
        return (
            np.any(maskArr[0, :, :])  or np.any(maskArr[-1, :, :]) or
            np.any(maskArr[:, 0, :])  or np.any(maskArr[:, -1, :]) or
            np.any(maskArr[:, :, 0])  or np.any(maskArr[:, :, -1])
        )

    islas = {}
    for i in range(nSegments):
        segId = segmentacion.GetNthSegmentID(i)
        maskArr = slicer.util.arrayFromSegmentBinaryLabelmap(segmentationNode, segId, volumeNode).astype(bool)
        n_voxels = np.count_nonzero(maskArr)
        if n_voxels == 0:
            continue
        islas[segId] = {
            "mask": maskArr,
            "vol": n_voxels * voxel_vol_cm3,
            "tocaBorde": toca_borde(maskArr),
        }

    nDescartadasPorBorde = sum(1 for d in islas.values() if d["tocaBorde"])
    print(f"{nDescartadasPorBorde} isla(s) descartada(s) por tocar el borde del volumen (camilla/colchoneta/etc).")

    pool = {segId: d for segId, d in islas.items() if not d["tocaBorde"]}

    if len(pool) == 0:
        print("ERROR: todas las islas tocan el borde del volumen. Revisa el FOV o el umbral.")
    else:
        ref_segId = max(pool, key=lambda k: pool[k]["vol"])
        print(f"\nIsla de referencia (mayor volumen, sin tocar borde): vol={pool[ref_segId]['vol']:.2f}cm3")

        aceptadas = {ref_segId}
        mascara_aceptada = pool[ref_segId]["mask"].copy()
        pendientes = {k: v for k, v in pool.items() if k != ref_segId}

        ronda = 0
        while True:
            ronda += 1
            distMap = ndimage.distance_transform_edt(~mascara_aceptada, sampling=sampling_zyx)
            nuevas = []
            for segId, d in pendientes.items():
                dist_min = distMap[d["mask"]].min()
                if dist_min <= MARGEN_CENTRALIDAD_MM:
                    nuevas.append(segId)
            if not nuevas:
                break
            for segId in nuevas:
                mascara_aceptada |= pendientes[segId]["mask"]
                aceptadas.add(segId)
                del pendientes[segId]
            print(f"  Ronda {ronda}: se sumaron {len(nuevas)} isla(s) por proximidad.")

        candidatas = [(segId, pool[segId]["vol"]) for segId in aceptadas]

        volumen_mayor = max(v for _, v in candidatas)
        candidatasFinales = [(segId, vol) for segId, vol in candidatas if vol >= volumen_mayor * UMBRAL_RELATIVO]

        idsCandidatas = set(segId for segId, _ in candidatasFinales)
        for i in range(nSegments - 1, -1, -1):
            segId = segmentacion.GetNthSegmentID(i)
            if segId not in idsCandidatas:
                segmentacion.RemoveSegment(segId)

        print(f"\nQuedaron {len(candidatasFinales)} islas candidatas (pasaron borde + proximidad + tamano).")
        print("Se detiene aca para revision manual.\n")

        # ============================================================
        # PARTE B: revision manual antes de fusionar
        # ============================================================

        islasRevision = []
        coloresOriginales = {}

        print("--- ISLAS CANDIDATAS (revisar antes de confirmar) ---")
        for idx, (segId, vol) in enumerate(sorted(candidatasFinales, key=lambda c: c[1], reverse=True), start=1):
            seg = segmentacion.GetSegment(segId)
            seg.SetName(f"Isla_{idx}")
            coloresOriginales[segId] = seg.GetColor()
            islasRevision.append({"numero": idx, "segId": segId, "vol": vol})
            print(f"  [{idx}] vol={vol:.2f}cm3")
        print("------------------------------------------------------\n")

        segmentationNode.CreateClosedSurfaceRepresentation()
        layoutManager = slicer.app.layoutManager()
        threeDView = layoutManager.threeDWidget(0).threeDView()
        threeDView.resetFocalPoint()

        def _segIdPorNumero(numero):
            for item in islasRevision:
                if item["numero"] == numero:
                    return item["segId"]
            print(f"No existe la isla numero {numero}.")
            return None

        def resaltar(numero):
            segId = _segIdPorNumero(numero)
            if segId is None:
                return
            for item in islasRevision:
                seg = segmentacion.GetSegment(item["segId"])
                if item["segId"] == segId:
                    seg.SetColor(1.0, 0.0, 0.0)
                else:
                    seg.SetColor(*coloresOriginales[item["segId"]])
            print(f"Isla {numero} resaltada en rojo.")

        def eliminar(numero):
            segId = _segIdPorNumero(numero)
            if segId is None:
                return
            segmentacion.RemoveSegment(segId)
            islasRevision[:] = [it for it in islasRevision if it["segId"] != segId]
            print(f"Isla {numero} eliminada. Quedan: {[it['numero'] for it in islasRevision]}")

        def mostrar_todas():
            for item in islasRevision:
                seg = segmentacion.GetSegment(item["segId"])
                seg.SetColor(*coloresOriginales[item["segId"]])
            print("Colores restaurados.")

        def confirmar_craneo():
            if len(islasRevision) == 0:
                print("No queda ninguna isla. Revisa si eliminaste de mas.")
                return
            idsRestantes = [item["segId"] for item in islasRevision]
            if len(idsRestantes) > 1:
                segmentEditorWidget2 = slicer.qMRMLSegmentEditorWidget()
                segmentEditorWidget2.setMRMLScene(slicer.mrmlScene)
                segEditorNode2 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
                segmentEditorWidget2.setMRMLSegmentEditorNode(segEditorNode2)
                segmentEditorWidget2.setSegmentationNode(segmentationNode)
                segmentEditorWidget2.setSourceVolumeNode(volumeNode)

                primero = idsRestantes[0]
                segmentEditorWidget2.setCurrentSegmentID(primero)
                segmentEditorWidget2.setActiveEffectByName("Logical operators")
                logicalEffect = segmentEditorWidget2.activeEffect()
                for otro in idsRestantes[1:]:
                    logicalEffect.setParameter("Operation", "UNION")
                    logicalEffect.setParameter("ModifierSegmentID", otro)
                    logicalEffect.self().onApply()
                    segmentacion.RemoveSegment(otro)
                segmentacion.GetSegment(primero).SetName("Craneo_Final")
                slicer.mrmlScene.RemoveNode(segEditorNode2)
            else:
                segmentacion.GetSegment(idsRestantes[0]).SetName("Craneo_Final")

            segmentationNode.CreateClosedSurfaceRepresentation()
            threeDView.resetFocalPoint()
            print("Craneo final confirmado y fusionado.")
            print("Ahora corre enviar_a_planner() para dejarlo listo para el Osteotomy Planner.")

        # ============================================================
        # PARTE PUENTE: Craneo_Final (segmento) -> Model node en carpeta SH
        # ============================================================

        def _buscar_carpeta(nombre, shNode, parentId):
            children = vtk.vtkIdList()
            shNode.GetItemChildren(parentId, children)
            for i in range(children.GetNumberOfIds()):
                itemId = children.GetId(i)
                if shNode.GetItemName(itemId) == nombre:
                    return itemId
            return None

        def _buscar_modelo_en_carpeta(nombreModelo, shNode, folderItemId):
            children = vtk.vtkIdList()
            shNode.GetItemChildren(folderItemId, children)
            for i in range(children.GetNumberOfIds()):
                itemId = children.GetId(i)
                node = shNode.GetItemDataNode(itemId)
                if node and node.IsA("vtkMRMLModelNode") and node.GetName() == nombreModelo:
                    return node
            return None

        def _limpiar_islas_residuales_malla(modelNode, umbralRelativo=UMBRAL_LIMPIEZA_MALLA):
            """
            NUEVO v5. Corre vtkPolyDataConnectivityFilter sobre la malla ya exportada
            y descarta regiones chicas (residuos de la conversion segmento->superficie).
            No toca el segmento original, solo el Model node exportado.
            """
            poly = modelNode.GetPolyData()
            if poly is None or poly.GetNumberOfPoints() == 0:
                print("Modelo sin geometria, no se puede limpiar.")
                return

            conectividad = vtk.vtkPolyDataConnectivityFilter()
            conectividad.SetInputData(poly)
            conectividad.SetExtractionModeToAllRegions()
            conectividad.ColorRegionsOn()
            conectividad.Update()
            numRegionesOriginal = conectividad.GetNumberOfExtractedRegions()

            if numRegionesOriginal <= 1:
                print(f"Malla de '{modelNode.GetName()}': 1 sola region conectada, no hace falta limpiar.")
                return

            regionArray = conectividad.GetOutput().GetPointData().GetArray("RegionId")
            conteoPorRegion = {}
            for i in range(regionArray.GetNumberOfTuples()):
                rid = int(regionArray.GetTuple1(i))
                conteoPorRegion[rid] = conteoPorRegion.get(rid, 0) + 1

            volumenMayor = max(conteoPorRegion.values())
            regionesAConservar = [rid for rid, cnt in conteoPorRegion.items() if cnt >= volumenMayor * umbralRelativo]

            conectividad.SetExtractionModeToSpecifiedRegions()
            conectividad.InitializeSpecifiedRegionList()
            for rid in regionesAConservar:
                conectividad.AddSpecifiedRegion(rid)
            conectividad.Update()

            limpiador = vtk.vtkCleanPolyData()
            limpiador.SetInputConnection(conectividad.GetOutputPort())
            limpiador.Update()

            modelNode.SetAndObservePolyData(limpiador.GetOutput())
            modelNode.GetDisplayNode().Modified() if modelNode.GetDisplayNode() else None

            print(f"Limpieza de malla en '{modelNode.GetName()}': {numRegionesOriginal} regiones -> "
                  f"{len(regionesAConservar)} conservadas (umbral {umbralRelativo:.0%} del tamano mayor).")
            if numRegionesOriginal - len(regionesAConservar) > 0:
                print(f"  ({numRegionesOriginal - len(regionesAConservar)} region(es) residual(es) descartada(s))")

        def enviar_a_planner(nombre_carpeta="CranioPlan"):
            """
            Exporta el segmento 'Craneo_Final' a un Model node dentro de una carpeta
            de Subject Hierarchy, y limpia islas residuales de la malla resultante.
            Llamar DESPUES de confirmar_craneo().
            """
            segId = segmentacion.GetSegmentIdBySegmentName("Craneo_Final")
            if not segId:
                print("No encontre el segmento 'Craneo_Final'. Corre confirmar_craneo() primero.")
                return

            shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
            sceneItemId = shNode.GetSceneItemID()

            folderItemId = _buscar_carpeta(nombre_carpeta, shNode, sceneItemId)
            if folderItemId is None:
                folderItemId = shNode.CreateFolderItem(sceneItemId, nombre_carpeta)
                print(f"Carpeta '{nombre_carpeta}' creada en Subject Hierarchy.")
            else:
                print(f"Reutilizando carpeta existente '{nombre_carpeta}'.")

            exito = slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
                segmentationNode, [segId], folderItemId)

            if exito:
                print(f"Craneo_Final exportado como modelo dentro de la carpeta '{nombre_carpeta}'.")

                modeloExportado = _buscar_modelo_en_carpeta("Craneo_Final", shNode, folderItemId)
                if modeloExportado is not None:
                    _limpiar_islas_residuales_malla(modeloExportado)
                else:
                    print("ADVERTENCIA: se exporto el modelo pero no lo pude encontrar por nombre "
                          "para limpiarlo. Revisa el nombre del nodo en el panel Data.")

                print("Ahora anda al modulo 'Osteotomy Planner', y en el combo de arriba")
                print(f"selecciona la carpeta '{nombre_carpeta}'. El craneo deberia aparecer")
                print("en el arbol de abajo, listo para seleccionar y empezar a cortar.")
            else:
                print("Fallo la exportacion. Revisa la consola de error de Slicer (Python console -> pestaña Error log).")

        print("Comandos disponibles:")
        print("  resaltar(numero)     -> pinta esa isla de rojo, el resto vuelve a su color")
        print("  eliminar(numero)     -> borra esa isla de la revision (ej: si es algo que se colo)")
        print("  mostrar_todas()      -> restaura colores normales")
        print("  confirmar_craneo()   -> fusiona lo que haya quedado en un Craneo_Final")
        print("  enviar_a_planner()   -> exporta Craneo_Final al formato que usa el Osteotomy Planner,")
        print("                          y limpia automaticamente islas residuales de la malla")

    segmentEditorWidget = None
