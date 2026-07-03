import numpy as np
import logging
import os
from typing import Annotated

import vtk
import qt
import ctk
import pydicom

import slicer
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import (
    parameterNodeWrapper,
    WithinRange,
)
from DICOMLib import DICOMUtils

from slicer import vtkMRMLScalarVolumeNode


#
# CranioPlan
#


class CranioPlan(ScriptedLoadableModule):

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("CranioPlan")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "Craneofacial")]
        self.parent.dependencies = []
        self.parent.contributors = ["Valentino Andri", "Ignacio Gross"]
        self.parent.helpText = _("""
Módulo de planificación prequirúrgica para craneosinostosis.
Trabajo final de grado, Ingeniería Biomédica, FCEFyN-UNC.
""")
        self.parent.acknowledgementText = _("""
Desarrollado en colaboración con el Servicio de Neurocirugía del
Hospital Garrahan, Buenos Aires.
""")


#
# CranioPlanParameterNode
#


@parameterNodeWrapper
class CranioPlanParameterNode:
    """
    Parámetros que se guardan junto con la escena de Slicer.
    estudioCargado - El volumen DICOM que se cargó en el Paso 1.
    """
    estudioCargado: vtkMRMLScalarVolumeNode = None


#
# CranioPlanWidget
#


class CranioPlanWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):

    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self._parameterNode = None
        self._parameterNodeGuiTag = None

        # Estado interno del Bloque B
        # Se llena cuando el algoritmo termina y quedan islas para revisar
        self._segmentationNode = None   # nodo de segmentación activo
        self._islasRevision = []        # lista de dicts: {numero, segId, vol, dist}
        self._coloresOriginales = {}    # segId -> (r,g,b) antes de resaltar
        self._botonesIslas = []         # lista de widgets de fila por isla (para limpiarlos)

    def setup(self) -> None:
        ScriptedLoadableModuleWidget.setup(self)

        self.logic = CranioPlanLogic()

        # -------------------------------------------------------
        # PASO 1 — Cargar estudio
        # -------------------------------------------------------
        pasoUnoCollapsible = ctk.ctkCollapsibleButton()
        pasoUnoCollapsible.text = "Paso 1 — Cargar estudio"
        self.layout.addWidget(pasoUnoCollapsible)
        pasoUnoLayout = qt.QVBoxLayout(pasoUnoCollapsible)

        instrucciones = qt.QLabel(
            "Seleccioná la carpeta con la tomografía del paciente.\n"
            "El sistema va a buscar y cargar el estudio automáticamente."
        )
        instrucciones.setWordWrap(True)
        pasoUnoLayout.addWidget(instrucciones)

        self.botonCargarEstudio = qt.QPushButton("Seleccionar carpeta del estudio...")
        self.botonCargarEstudio.toolTip = "Elegí la carpeta DICOM que te mandaron del Garrahan"
        pasoUnoLayout.addWidget(self.botonCargarEstudio)

        self.etiquetaEstadoCarga = qt.QLabel("Todavía no se cargó ningún estudio.")
        self.etiquetaEstadoCarga.setStyleSheet("color: gray;")
        pasoUnoLayout.addWidget(self.etiquetaEstadoCarga)

        self.botonCargarEstudio.connect("clicked(bool)", self.onBotonCargarEstudioClicked)

        # -------------------------------------------------------
        # PASO 2 — Generar cráneo 3D
        # Sub-paso 2A: botón para correr el algoritmo
        # Sub-paso 2B: panel de revisión de islas (se llena dinámicamente)
        # -------------------------------------------------------
        pasoDosCollapsible = ctk.ctkCollapsibleButton()
        pasoDosCollapsible.text = "Paso 2 — Generar cráneo 3D"
        self.layout.addWidget(pasoDosCollapsible)
        self.pasoDosLayout = qt.QVBoxLayout(pasoDosCollapsible)

        # 2A — Botón principal
        instruccionesPasoDos = qt.QLabel(
            "Con el estudio ya cargado, generá automáticamente el modelo 3D del cráneo.\n"
            "El sistema va a identificar las partes candidatas para que puedas "
            "revisarlas antes de confirmar."
        )
        instruccionesPasoDos.setWordWrap(True)
        self.pasoDosLayout.addWidget(instruccionesPasoDos)

        self.botonGenerarCraneo = qt.QPushButton("Generar cráneo 3D")
        self.botonGenerarCraneo.toolTip = "Aplica segmentación automática sobre el estudio cargado"
        self.pasoDosLayout.addWidget(self.botonGenerarCraneo)

        self.etiquetaEstadoCraneo = qt.QLabel("Todavía no se generó el cráneo.")
        self.etiquetaEstadoCraneo.setStyleSheet("color: gray;")
        self.pasoDosLayout.addWidget(self.etiquetaEstadoCraneo)

        self.botonGenerarCraneo.connect("clicked(bool)", self.onBotonGenerarCraneoClicked)

        # 2B — Panel de revisión (vacío al inicio, se llena al terminar 2A)
        self.separadorRevision = qt.QFrame()
        self.separadorRevision.setFrameShape(qt.QFrame.HLine)
        self.separadorRevision.setStyleSheet("color: #CCCCCC;")
        self.pasoDosLayout.addWidget(self.separadorRevision)
        self.separadorRevision.setVisible(False)

        self.etiquetaRevision = qt.QLabel("Revisión de islas candidatas:")
        self.etiquetaRevision.setStyleSheet("font-weight: bold;")
        self.pasoDosLayout.addWidget(self.etiquetaRevision)
        self.etiquetaRevision.setVisible(False)

        self.etiquetaAyudaRevision = qt.QLabel(
            "Usá 'Resaltar' para ver cada isla en rojo en el visor 3D.\n"
            "Si una isla es la camilla u otra estructura que no es cráneo, eliminala.\n"
            "Cuando estés conforme, presioná 'Confirmar cráneo'."
        )
        self.etiquetaAyudaRevision.setWordWrap(True)
        self.etiquetaAyudaRevision.setStyleSheet("color: gray; font-size: 9px;")
        self.pasoDosLayout.addWidget(self.etiquetaAyudaRevision)
        self.etiquetaAyudaRevision.setVisible(False)

        # Contenedor donde se van a generar los botones de islas dinámicamente
        self.contenedorIslas = qt.QWidget()
        self.layoutIslas = qt.QVBoxLayout(self.contenedorIslas)
        self.layoutIslas.setContentsMargins(0, 0, 0, 0)
        self.pasoDosLayout.addWidget(self.contenedorIslas)
        self.contenedorIslas.setVisible(False)

        # Botón confirmar (siempre al final del panel de revisión)
        self.botonConfirmarCraneo = qt.QPushButton("✓  Confirmar cráneo")
        self.botonConfirmarCraneo.toolTip = "Fusiona las islas restantes en un único Craneo_Final"
        self.botonConfirmarCraneo.setStyleSheet(
            "background-color: #1E7B45; color: white; font-weight: bold; padding: 6px;"
        )
        self.pasoDosLayout.addWidget(self.botonConfirmarCraneo)
        self.botonConfirmarCraneo.setVisible(False)
        self.botonConfirmarCraneo.connect("clicked(bool)", self.onBotonConfirmarCraneoClicked)

        self.etiquetaEstadoConfirmacion = qt.QLabel("")
        self.pasoDosLayout.addWidget(self.etiquetaEstadoConfirmacion)

        # -------------------------------------------------------
        # Espacio para próximos pasos (Paso 3, 4)
        # -------------------------------------------------------
        self.layout.addStretch(1)

        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        self.initializeParameterNode()

    def cleanup(self) -> None:
        self.removeObservers()

    def enter(self) -> None:
        self.initializeParameterNode()

    def exit(self) -> None:
        if self._parameterNode:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
            self._parameterNodeGuiTag = None

    def onSceneStartClose(self, caller, event) -> None:
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event) -> None:
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self) -> None:
        self.setParameterNode(self.logic.getParameterNode())

    def setParameterNode(self, inputParameterNode) -> None:
        if self._parameterNode:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
        self._parameterNode = inputParameterNode
        if self._parameterNode:
            self._parameterNodeGuiTag = None

    # -------------------------------------------------------
    # Paso 1 — handlers
    # -------------------------------------------------------

    def onBotonCargarEstudioClicked(self):
        carpetaSeleccionada = qt.QFileDialog.getExistingDirectory(
            self.parent, "Seleccioná la carpeta del estudio DICOM"
        )
        if not carpetaSeleccionada:
            return

        self.etiquetaEstadoCarga.setText("Cargando estudio, por favor esperá...")
        self.etiquetaEstadoCarga.setStyleSheet("color: orange;")
        slicer.app.processEvents()

        volumenNode, metodoUsado = self.logic.cargarCarpetaDicom(carpetaSeleccionada)

        if volumenNode is None:
            self.etiquetaEstadoCarga.setText(
                "No se encontró ninguna serie volumétrica válida en esa carpeta."
            )
            self.etiquetaEstadoCarga.setStyleSheet("color: red;")
            return

        self._parameterNode.estudioCargado = volumenNode

        if metodoUsado == "hueso_explicito":
            self.etiquetaEstadoCarga.setText(f"Estudio cargado: {volumenNode.GetName()}")
            self.etiquetaEstadoCarga.setStyleSheet("color: green;")
        else:
            self.etiquetaEstadoCarga.setText(
                f"Estudio cargado: {volumenNode.GetName()}\n"
                "(no había serie de \"hueso\" explícita; se usó la mejor "
                "serie volumétrica disponible)"
            )
            self.etiquetaEstadoCarga.setStyleSheet("color: #B8860B;")

    # -------------------------------------------------------
    # Paso 2A — Generar candidatas
    # -------------------------------------------------------

    def onBotonGenerarCraneoClicked(self):
        volumenActual = self._parameterNode.estudioCargado

        if volumenActual is None:
            self.etiquetaEstadoCraneo.setText("Primero tenés que cargar un estudio (Paso 1).")
            self.etiquetaEstadoCraneo.setStyleSheet("color: red;")
            return

        self.etiquetaEstadoCraneo.setText("Generando cráneo 3D, por favor esperá...")
        self.etiquetaEstadoCraneo.setStyleSheet("color: orange;")
        slicer.app.processEvents()

        resultado = self.logic.generarCandidatas(volumenActual)

        if resultado is None:
            self.etiquetaEstadoCraneo.setText(
                "No se pudo generar el cráneo. Revisá el estudio cargado."
            )
            self.etiquetaEstadoCraneo.setStyleSheet("color: red;")
            return

        self._segmentationNode, self._islasRevision, self._coloresOriginales = resultado

        n = len(self._islasRevision)
        self.etiquetaEstadoCraneo.setText(
            f"Se encontraron {n} isla(s) candidata(s). "
            f"{'Revisalas antes de confirmar.' if n > 1 else 'Una sola isla — podés confirmar directamente.'}"
        )
        self.etiquetaEstadoCraneo.setStyleSheet("color: #1F4E79;")

        self._construirPanelRevision()

    def _construirPanelRevision(self):
        """Arma dinámicamente los botones de revisión para cada isla candidata."""

        # Limpiar botones anteriores si hubiera
        while self.layoutIslas.count():
            item = self.layoutIslas.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._botonesIslas = []

        for isla in self._islasRevision:
            numero = isla["numero"]
            vol = isla["vol"]
            dist = isla["dist"]

            filaWidget = qt.QWidget()
            filaLayout = qt.QHBoxLayout(filaWidget)
            filaLayout.setContentsMargins(0, 2, 0, 2)

            etiqueta = qt.QLabel(f"Isla {numero}  |  {vol:.1f} cm³  |  dist. centro: {dist:.0f} mm")
            etiqueta.setStyleSheet("font-size: 10px;")
            filaLayout.addWidget(etiqueta, 2)

            botonResaltar = qt.QPushButton("Resaltar")
            botonResaltar.setStyleSheet("padding: 3px 8px;")
            botonResaltar.setFixedWidth(70)
            filaLayout.addWidget(botonResaltar)

            botonEliminar = qt.QPushButton("Eliminar")
            botonEliminar.setStyleSheet(
                "padding: 3px 8px; background-color: #C0392B; color: white;"
            )
            botonEliminar.setFixedWidth(70)
            filaLayout.addWidget(botonEliminar)

            self.layoutIslas.addWidget(filaWidget)
            self._botonesIslas.append({
                "numero": numero,
                "filaWidget": filaWidget,
                "botonResaltar": botonResaltar,
                "botonEliminar": botonEliminar,
            })

            # Conectamos con lambdas que capturan el número de isla
            botonResaltar.connect(
                "clicked(bool)",
                lambda _, n=numero: self.onResaltarIsla(n)
            )
            botonEliminar.connect(
                "clicked(bool)",
                lambda _, n=numero: self.onEliminarIsla(n)
            )

        # Mostrar el panel de revisión
        self.separadorRevision.setVisible(True)
        self.etiquetaRevision.setVisible(True)
        self.etiquetaAyudaRevision.setVisible(True)
        self.contenedorIslas.setVisible(True)
        self.botonConfirmarCraneo.setVisible(True)
        self.etiquetaEstadoConfirmacion.setText("")

    # -------------------------------------------------------
    # Paso 2B — Revisión de islas
    # -------------------------------------------------------

    def onResaltarIsla(self, numero):
        """Pinta la isla seleccionada de rojo, las demás vuelven a su color original."""
        segmentacion = self._segmentationNode.GetSegmentation()
        for isla in self._islasRevision:
            seg = segmentacion.GetSegment(isla["segId"])
            if seg is None:
                continue
            if isla["numero"] == numero:
                seg.SetColor(1.0, 0.2, 0.2)
            else:
                color = self._coloresOriginales.get(isla["segId"], (0.5, 0.5, 0.5))
                seg.SetColor(*color)

    def onEliminarIsla(self, numero):
        """Elimina la isla del segmento y la quita del panel de revisión."""
        isla = next((i for i in self._islasRevision if i["numero"] == numero), None)
        if isla is None:
            return

        segmentacion = self._segmentationNode.GetSegmentation()
        segmentacion.RemoveSegment(isla["segId"])

        self._islasRevision = [i for i in self._islasRevision if i["numero"] != numero]
        self._coloresOriginales.pop(isla["segId"], None)

        # Quitar la fila de botones de la UI
        fila = next((b for b in self._botonesIslas if b["numero"] == numero), None)
        if fila:
            fila["filaWidget"].setVisible(False)
            self._botonesIslas = [b for b in self._botonesIslas if b["numero"] != numero]

        n = len(self._islasRevision)
        self.etiquetaEstadoCraneo.setText(
            f"Isla {numero} eliminada. Quedan {n} isla(s)."
        )
        self.etiquetaEstadoCraneo.setStyleSheet("color: #B8860B;")

    def onBotonConfirmarCraneoClicked(self):
        """Fusiona las islas restantes en un único Craneo_Final."""
        if not self._islasRevision:
            self.etiquetaEstadoConfirmacion.setText(
                "No quedan islas. Volvé a generar el cráneo."
            )
            self.etiquetaEstadoConfirmacion.setStyleSheet("color: red;")
            return

        self.etiquetaEstadoConfirmacion.setText("Confirmando, por favor esperá...")
        self.etiquetaEstadoConfirmacion.setStyleSheet("color: orange;")
        slicer.app.processEvents()

        volumenActual = self._parameterNode.estudioCargado
        exito = self.logic.confirmarCraneo(
            self._segmentationNode,
            self._islasRevision,
            volumenActual
        )

        if not exito:
            self.etiquetaEstadoConfirmacion.setText("Error al confirmar el cráneo.")
            self.etiquetaEstadoConfirmacion.setStyleSheet("color: red;")
            return

        self.etiquetaEstadoConfirmacion.setText("¡Cráneo final confirmado!")
        self.etiquetaEstadoConfirmacion.setStyleSheet(
            "color: #1E7B45; font-weight: bold;"
        )

        # Ocultar el panel de revisión
        self.separadorRevision.setVisible(False)
        self.etiquetaRevision.setVisible(False)
        self.etiquetaAyudaRevision.setVisible(False)
        self.contenedorIslas.setVisible(False)
        self.botonConfirmarCraneo.setVisible(False)

        # Centrar la vista 3D sobre el resultado
        layoutManager = slicer.app.layoutManager()
        threeDWidget = layoutManager.threeDWidget(0)
        threeDWidget.threeDView().resetFocalPoint()


#
# CranioPlanLogic
#


class CranioPlanLogic(ScriptedLoadableModuleLogic):
    """Funciones de procesamiento del módulo CranioPlan, sin interfaz."""

    def __init__(self) -> None:
        ScriptedLoadableModuleLogic.__init__(self)

    def getParameterNode(self):
        return CranioPlanParameterNode(super().getParameterNode())

    # ============================================================
    # BLOQUE C — Identificación automática de la serie de hueso
    # ============================================================

    def identificarSerieDeHueso(self, db, seriesUIDs):
        """
        Devuelve (seriesUID, metodoUsado) o (None, None).
        metodoUsado: "hueso_explicito" | "fallback_volumetrico"

        Regla en cascada:
        1. Descartar planos fijos (axial/coronal/sagittal) y sin SliceThickness.
        2. Preferir las que digan "hueso"/"bone" en la descripción.
        3. Si no hay ninguna con esa palabra, usar la de menor SliceThickness
           entre todas las volumétricas válidas (fallback).

        Ajustada el 28/06/2026 tras confirmar con 2 casos reales del Garrahan
        que no siempre hay una serie etiquetada explícitamente como "hueso".
        """
        palabrasClaveHueso = ["hueso", "bone"]
        palabrasClavePlanoFijo = ["axial", "coronal", "sagittal"]
        volumetricasValidas = []

        for seriesUID in seriesUIDs:
            archivos = db.filesForSeries(seriesUID)
            if not archivos:
                continue
            ds = pydicom.dcmread(archivos[0], stop_before_pixels=True)
            descripcion = str(getattr(ds, "SeriesDescription", "")).lower()

            if any(p in descripcion for p in palabrasClavePlanoFijo):
                continue
            espesorCorte = getattr(ds, "SliceThickness", None)
            if espesorCorte is None:
                continue

            tieneHueso = any(p in descripcion for p in palabrasClaveHueso)
            volumetricasValidas.append((seriesUID, float(espesorCorte), len(archivos), tieneHueso))

        if not volumetricasValidas:
            return None, None

        candidatasIdeales = [v for v in volumetricasValidas if v[3]]
        if candidatasIdeales:
            candidatasIdeales.sort(key=lambda c: (c[1], -c[2]))
            return candidatasIdeales[0][0], "hueso_explicito"

        volumetricasValidas.sort(key=lambda c: (c[1], -c[2]))
        return volumetricasValidas[0][0], "fallback_volumetrico"

    def cargarCarpetaDicom(self, rutaCarpeta):
        """
        Importa una carpeta DICOM y carga la serie elegida por
        identificarSerieDeHueso.
        Devuelve (volumenNode, metodoUsado) o (None, None).
        """
        with DICOMUtils.TemporaryDICOMDatabase() as db:
            DICOMUtils.importDicom(rutaCarpeta, db)

            patientUIDs = db.patients()
            if not patientUIDs:
                return None, None

            studyUIDs = db.studiesForPatient(patientUIDs[0])
            if not studyUIDs:
                return None, None

            seriesUIDs = db.seriesForStudy(studyUIDs[0])
            if not seriesUIDs:
                return None, None

            serieElegidaUID, metodoUsado = self.identificarSerieDeHueso(db, seriesUIDs)
            if serieElegidaUID is None:
                return None, None

            loadedNodeIDs = DICOMUtils.loadSeriesByUID([serieElegidaUID])
            for nodeID in loadedNodeIDs:
                node = slicer.mrmlScene.GetNodeByID(nodeID)
                if node and node.IsA("vtkMRMLScalarVolumeNode"):
                    return node, metodoUsado

        return None, None

    # ============================================================
    # BLOQUE A + B — Segmentación automática y revisión manual
    # ============================================================

    def generarCandidatas(self, volumeNode):
        """
        BLOQUE A: Corre threshold + islands + filtro de centralidad/tamaño.
        A diferencia de generarCraneoDesdeVolumen (que fusionaba automáticamente),
        esta función se DETIENE antes de fusionar y devuelve las islas candidatas
        para que el Widget permita la revisión manual del Bloque B.

        Devuelve una tupla (segmentationNode, islasRevision, coloresOriginales),
        o None si algo falla.

        islasRevision: lista de dicts {numero, segId, vol, dist}
        coloresOriginales: dict {segId: (r, g, b)}
        """
        if volumeNode is None:
            return None

        UMBRAL_RELATIVO = 0.05
        MARGEN_CENTRALIDAD_MM = 40.0

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

        # Threshold 300-3000 HU
        segmentEditorWidget.setActiveEffectByName("Threshold")
        thresholdEffect = segmentEditorWidget.activeEffect()
        thresholdEffect.setParameter("MinimumThreshold", "300")
        thresholdEffect.setParameter("MaximumThreshold", "3000")
        thresholdEffect.self().onApply()

        # Islands: separar regiones conectadas
        segmentEditorWidget.setActiveEffectByName("Islands")
        islandsEffect = segmentEditorWidget.activeEffect()
        islandsEffect.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
        islandsEffect.self().onApply()

        segmentacion = segmentationNode.GetSegmentation()
        nSegments = segmentacion.GetNumberOfSegments()

        if nSegments == 0:
            segmentEditorWidget = None
            return None

        # Info geométrica para IJK -> RAS
        ijkToRas = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijkToRas)

        def ijk_a_ras(i, j, k):
            p = ijkToRas.MultiplyPoint([i, j, k, 1])
            return p[0], p[1], p[2]

        spacing = volumeNode.GetSpacing()
        voxel_vol_cm3 = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

        datos = []
        for i in range(nSegments):
            segId = segmentacion.GetNthSegmentID(i)
            arr = slicer.util.arrayFromSegmentBinaryLabelmap(segmentationNode, segId, volumeNode)
            n_voxels = np.count_nonzero(arr)
            if n_voxels == 0:
                continue
            vol_cm3 = n_voxels * voxel_vol_cm3
            zs, ys, xs = np.where(arr > 0)
            cx, cy, cz = ijk_a_ras(xs.mean(), ys.mean(), zs.mean())
            datos.append((segId, vol_cm3, cx, cy, cz))

        if not datos:
            segmentEditorWidget = None
            return None

        datos.sort(key=lambda d: d[1], reverse=True)
        ref_segId, ref_vol, ref_cx, ref_cy, ref_cz = datos[0]

        candidatas = []
        for segId, vol, cx, cy, cz in datos:
            dist_xy = ((cx - ref_cx) ** 2 + (cy - ref_cy) ** 2) ** 0.5
            if dist_xy <= MARGEN_CENTRALIDAD_MM:
                candidatas.append((segId, vol, dist_xy))

        if not candidatas:
            segmentEditorWidget = None
            return None

        volumen_mayor = max(c[1] for c in candidatas)
        candidatasFinales = [c for c in candidatas if c[1] >= volumen_mayor * UMBRAL_RELATIVO]

        # Eliminar del segmento todo lo que no sea candidata final
        idsCandidatas = set(c[0] for c in candidatasFinales)
        for i in range(nSegments - 1, -1, -1):
            segId = segmentacion.GetNthSegmentID(i)
            if segId not in idsCandidatas:
                segmentacion.RemoveSegment(segId)

        # Preparar estructura de revisión para el Bloque B (Widget)
        islasRevision = []
        coloresOriginales = {}

        for idx, (segId, vol, dist) in enumerate(
            sorted(candidatasFinales, key=lambda c: c[1], reverse=True), start=1
        ):
            seg = segmentacion.GetSegment(segId)
            seg.SetName(f"Isla_{idx}")
            coloresOriginales[segId] = seg.GetColor()
            islasRevision.append({"numero": idx, "segId": segId, "vol": vol, "dist": dist})

        segmentationNode.CreateClosedSurfaceRepresentation()
        segmentEditorWidget = None

        return segmentationNode, islasRevision, coloresOriginales

    def confirmarCraneo(self, segmentationNode, islasRevision, volumeNode):
        """
        BLOQUE B: fusiona las islas que quedaron en islasRevision en un único
        segmento llamado "Craneo_Final". Equivale al confirmar_craneo() del
        script de consola de Nacho.

        Devuelve True si tuvo éxito, False si algo falló.
        """
        if not islasRevision:
            return False

        segmentacion = segmentationNode.GetSegmentation()
        idsRestantes = [isla["segId"] for isla in islasRevision]

        if len(idsRestantes) > 1:
            segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
            segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
            segEditorNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
            segmentEditorWidget.setMRMLSegmentEditorNode(segEditorNode)
            segmentEditorWidget.setSegmentationNode(segmentationNode)
            segmentEditorWidget.setSourceVolumeNode(volumeNode)

            primero = idsRestantes[0]
            segmentEditorWidget.setCurrentSegmentID(primero)
            segmentEditorWidget.setActiveEffectByName("Logical operators")
            logicalEffect = segmentEditorWidget.activeEffect()

            for otro in idsRestantes[1:]:
                logicalEffect.setParameter("Operation", "UNION")
                logicalEffect.setParameter("ModifierSegmentID", otro)
                logicalEffect.self().onApply()
                segmentacion.RemoveSegment(otro)

            segmentacion.GetSegment(primero).SetName("Craneo_Final")
            slicer.mrmlScene.RemoveNode(segEditorNode)
            segmentEditorWidget = None
        else:
            segmentacion.GetSegment(idsRestantes[0]).SetName("Craneo_Final")

        segmentationNode.CreateClosedSurfaceRepresentation()
        return True


#
# CranioPlanTest
#


class CranioPlanTest(ScriptedLoadableModuleTest):

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.delayDisplay("Sin tests automatizados por ahora.")



