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
        self._segmentationNode = None
        self._islasRevision = []
        self._coloresOriginales = {}
        self._botonesIslas = []

        # Estado interno del Bloque F (planificación de osteotomías)
        self._modeloCraneoActual = None
        self._fragmentosActuales = []  # hueso vigente: se actualiza tras cada corte
        self._curvaCorteActual = None
        self._curvaEsCerrada = False  # se fija al trazar, según el checkbox
        self._observadorCurvaTag = None  # para seguir los puntos en tiempo real
        self.placeWidgetCorte = None  # qSlicerMarkupsPlaceWidget, se crea en setup()

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
        # -------------------------------------------------------
        pasoDosCollapsible = ctk.ctkCollapsibleButton()
        pasoDosCollapsible.text = "Paso 2 — Generar cráneo 3D"
        self.layout.addWidget(pasoDosCollapsible)
        self.pasoDosLayout = qt.QVBoxLayout(pasoDosCollapsible)

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

        self.contenedorIslas = qt.QWidget()
        self.layoutIslas = qt.QVBoxLayout(self.contenedorIslas)
        self.layoutIslas.setContentsMargins(0, 0, 0, 0)
        self.pasoDosLayout.addWidget(self.contenedorIslas)
        self.contenedorIslas.setVisible(False)

        self.botonConfirmarCraneo = qt.QPushButton("Confirmar cráneo")
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
        # PASO 3 — Planificar osteotomía
        # -------------------------------------------------------
        pasoTresCollapsible = ctk.ctkCollapsibleButton()
        pasoTresCollapsible.text = "Paso 3 — Planificar osteotomía"
        self.layout.addWidget(pasoTresCollapsible)
        pasoTresLayout = qt.QVBoxLayout(pasoTresCollapsible)

        instruccionesPasoTres = qt.QLabel(
            "Trazá la línea de corte sobre la superficie del cráneo, con clicks "
            "izquierdos siguiendo el camino deseado.\n"
            "No hace falta cerrar la línea ni volver al punto inicial: cuando "
            "termines, presioná directamente 'Finalizar trazado'.\n"
            "El sistema calcula automáticamente la profundidad necesaria para "
            "atravesar el hueso en cada punto."
        )
        instruccionesPasoTres.setWordWrap(True)
        pasoTresLayout.addWidget(instruccionesPasoTres)

        self.botonTrazarCorte = qt.QPushButton("Trazar línea de corte")
        pasoTresLayout.addWidget(self.botonTrazarCorte)

        self.checkCurvaCerrada = qt.QCheckBox("Curva cerrada (para aislar una región de hueso)")
        self.checkCurvaCerrada.setToolTip(
            "Sin tildar: línea abierta. Solo separa si sus dos extremos llegan\n"
            "a un borde real del cráneo (una órbita, el foramen magnum, etc.).\n"
            "Tildado: lazo cerrado. Aísla siempre la región que encierra,\n"
            "aunque esté en el medio del hueso, sin tocar ningún borde.\n\n"
            "Con esta opción tildada, NO cierres el lazo a mano: colocá los\n"
            "puntos del contorno y frená. El sistema une el último con el\n"
            "primero automáticamente."
        )
        pasoTresLayout.addWidget(self.checkCurvaCerrada)

        self.botonFinalizarTrazado = qt.QPushButton("Finalizar trazado")
        self.botonFinalizarTrazado.enabled = False  # se activa al entrar en modo trazado
        pasoTresLayout.addWidget(self.botonFinalizarTrazado)

        filaGrosor = qt.QHBoxLayout()
        etiquetaGrosor = qt.QLabel("Grosor de la osteotomía (mm):")
        filaGrosor.addWidget(etiquetaGrosor)
        self.spinGrosor = qt.QDoubleSpinBox()
        self.spinGrosor.setRange(0.1, 5.0)
        self.spinGrosor.setSingleStep(0.1)
        self.spinGrosor.setValue(1.0)
        self.spinGrosor.setToolTip(
            "Ancho de la hoja/sierra real usada en la cirugía.\n"
            "Valor pendiente de confirmar con el equipo del Garrahan."
        )
        filaGrosor.addWidget(self.spinGrosor)
        pasoTresLayout.addLayout(filaGrosor)

        self.botonGenerarCorte = qt.QPushButton("Generar corte")
        self.botonGenerarCorte.enabled = False
        pasoTresLayout.addWidget(self.botonGenerarCorte)

        self.etiquetaEstadoCorte = qt.QLabel("Todavía no se planificó ningún corte.")
        self.etiquetaEstadoCorte.setStyleSheet("color: gray;")
        pasoTresLayout.addWidget(self.etiquetaEstadoCorte)

        self.botonTrazarCorte.connect("clicked(bool)", self.onBotonTrazarCorteClicked)
        self.botonFinalizarTrazado.connect("clicked(bool)", self.onBotonFinalizarTrazadoClicked)
        self.botonGenerarCorte.connect("clicked(bool)", self.onBotonGenerarCorteClicked)

        # Widget oficial de Slicer para manejar la colocación de puntos.
        # Se usa "headless" (sin sus botones propios, ocultos) porque ya
        # tenemos nuestros propios botones en español; solo aprovechamos
        # su mecanismo interno, que es más confiable que armar a mano la
        # conexión entre clicks del mouse y el nodo activo.
        self.placeWidgetCorte = slicer.qSlicerMarkupsPlaceWidget()
        self.placeWidgetCorte.setMRMLScene(slicer.mrmlScene)
        self.placeWidgetCorte.buttonsVisible = False
        pasoTresLayout.addWidget(self.placeWidgetCorte)
        self.placeWidgetCorte.hide()

        # -------------------------------------------------------
        # Espacio para próximos pasos (Paso 4)
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

            etiqueta = qt.QLabel(f"Isla {numero}  |  {vol:.1f} cm3  |  dist. centro: {dist:.0f} mm")
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

            botonResaltar.connect(
                "clicked(bool)",
                lambda _, n=numero: self.onResaltarIsla(n)
            )
            botonEliminar.connect(
                "clicked(bool)",
                lambda _, n=numero: self.onEliminarIsla(n)
            )

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
        isla = next((i for i in self._islasRevision if i["numero"] == numero), None)
        if isla is None:
            return

        segmentacion = self._segmentationNode.GetSegmentation()
        segmentacion.RemoveSegment(isla["segId"])

        self._islasRevision = [i for i in self._islasRevision if i["numero"] != numero]
        self._coloresOriginales.pop(isla["segId"], None)

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

        self.etiquetaEstadoConfirmacion.setText("Cráneo final confirmado.")
        self.etiquetaEstadoConfirmacion.setStyleSheet(
            "color: #1E7B45; font-weight: bold;"
        )

        self.separadorRevision.setVisible(False)
        self.etiquetaRevision.setVisible(False)
        self.etiquetaAyudaRevision.setVisible(False)
        self.contenedorIslas.setVisible(False)
        self.botonConfirmarCraneo.setVisible(False)

        # Exportamos automáticamente a modelo 3D, para que el Paso 3
        # (planificar osteotomía) ya tenga con qué trabajar sin que
        # el usuario tenga que pasar por el módulo Segmentations.
        self._modeloCraneoActual = self.logic.exportarSegmentoAModelo(
            self._segmentationNode, nombreSegmento="Craneo_Final"
        )
        # El hueso vigente arranca siendo el cráneo entero. Cada corte
        # va a reemplazar esta lista por los fragmentos resultantes, así
        # el corte siguiente opera sobre el hueso ya cortado.
        self._fragmentosActuales = (
            [self._modeloCraneoActual] if self._modeloCraneoActual else []
        )

        # A partir de acá se trabaja sobre el MODELO, no sobre la
        # segmentación. Si dejamos las dos visibles, Slicer renderiza el
        # cráneo dos veces en cada movimiento del mouse — es una de las
        # causas principales del lag al rotar la vista 3D. Ocultar la
        # segmentación no pierde nada: sus datos siguen en la escena.
        displaySegmentacion = self._segmentationNode.GetDisplayNode()
        if displaySegmentacion is not None:
            displaySegmentacion.SetVisibility(False)

        layoutManager = slicer.app.layoutManager()
        threeDWidget = layoutManager.threeDWidget(0)
        threeDWidget.threeDView().resetFocalPoint()

    # -------------------------------------------------------
    # Paso 3 — Planificar osteotomía
    # -------------------------------------------------------

    def onBotonTrazarCorteClicked(self):
        if not self._fragmentosActuales:
            self.etiquetaEstadoCorte.setText("Primero confirmá el cráneo (Paso 2).")
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return

        # Si había una curva anterior a medio trazar, la limpiamos
        self._quitarObservadorCurva()
        if self._curvaCorteActual is not None:
            try:
                slicer.mrmlScene.RemoveNode(self._curvaCorteActual)
            except Exception:
                pass
            self._curvaCorteActual = None

        # Creamos la curva de osteotomía.
        # NOTA: vtkMRMLMarkupsCurveNode (abierta) nace abierta por
        # definición, sin setter de "closed" — solo GetCurveClosed()
        # para consultar. Para una curva CERRADA usamos directamente la
        # clase vtkMRMLMarkupsClosedCurveNode, que sí cierra el lazo
        # automáticamente. La Logic ya maneja ambos casos por igual en
        # _construirParedDeCorte() según GetCurveClosed().
        claseNodo = (
            "vtkMRMLMarkupsClosedCurveNode"
            if self.checkCurvaCerrada.checked
            else "vtkMRMLMarkupsCurveNode"
        )
        curvaNode = slicer.mrmlScene.AddNewNodeByClass(claseNodo, "Osteotomia_1")
        curvaNode.CreateDefaultDisplayNodes()
        displayNode = curvaNode.GetDisplayNode()
        if displayNode is not None:
            displayNode.SetSelectedColor(1.0, 0.2, 0.2)   # rojo al estar seleccionada
            displayNode.SetColor(1.0, 0.4, 0.4)           # rojo claro
            displayNode.SetGlyphScale(2.5)                # tamaño de los puntos
            displayNode.SetLineThickness(0.5)             # grosor de la línea que une los puntos
            displayNode.SetPropertiesLabelVisibility(False)  # oculta el texto de propiedades
            displayNode.SetPointLabelsVisibility(False)       # oculta las etiquetas 5-2, 5-3, etc.

        self._curvaCorteActual = curvaNode
        # Guardamos si es cerrada AHORA (al trazarla), no al generar el
        # corte: si el usuario toca el checkbox en el medio, lo que vale
        # es con qué tipo de curva se trazó realmente.
        self._curvaEsCerrada = bool(self.checkCurvaCerrada.checked)

        # Observador: cada vez que se agrega/mueve un punto, refrescamos
        # la etiqueta de estado con la cantidad de puntos colocados.
        self._observadorCurvaTag = curvaNode.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self._onPuntoAgregadoALaCurva
        )

        # Delegamos la colocación de puntos al widget oficial de Slicer
        # (qSlicerMarkupsPlaceWidget), en vez de manejar nosotros mismos
        # los nodos de interacción/selección. Este widget garantiza que
        # el nodo pasado a setCurrentNode quede correctamente activo
        # mientras dura la colocación, evitando que los puntos se
        # pierdan o queden en otro nodo.
        self.placeWidgetCorte.setCurrentNode(curvaNode)
        self.placeWidgetCorte.setPlaceModePersistency(True)
        self.placeWidgetCorte.setPlaceModeEnabled(True)

        self.botonFinalizarTrazado.enabled = True
        self.botonGenerarCorte.enabled = False
        self.etiquetaEstadoCorte.setText(
            "Modo de trazado activo. Hacé click sobre el cráneo en el visor 3D "
            "para colocar los puntos del corte (0 hasta ahora)."
        )
        self.etiquetaEstadoCorte.setStyleSheet("color: orange;")

    def _onPuntoAgregadoALaCurva(self, caller, event):
        """Refresca la etiqueta de estado en tiempo real al colocar puntos."""
        if self._curvaCorteActual is None:
            return
        n = self._curvaCorteActual.GetNumberOfControlPoints()
        self.etiquetaEstadoCorte.setText(
            f"Modo de trazado activo. Puntos colocados: {n}. "
            "Cuando termines, presioná 'Finalizar trazado'."
        )
        self.etiquetaEstadoCorte.setStyleSheet("color: orange;")

    def _quitarObservadorCurva(self):
        """Saca el observador de la curva, si estaba puesto."""
        if self._observadorCurvaTag is not None and self._curvaCorteActual is not None:
            try:
                self._curvaCorteActual.RemoveObserver(self._observadorCurvaTag)
            except Exception:
                pass
        self._observadorCurvaTag = None

    def onBotonFinalizarTrazadoClicked(self):
        # Apagamos el modo de colocación a través del widget oficial —
        # es lo que garantiza que el nodo quede correctamente finalizado
        # (en vez de tocar interactionNode a mano, que fue lo que fallaba).
        self.placeWidgetCorte.setPlaceModeEnabled(False)
        self._quitarObservadorCurva()

        numeroPuntos = 0 if self._curvaCorteActual is None else self._curvaCorteActual.GetNumberOfControlPoints()

        if numeroPuntos < 2:
            # Diagnóstico: si igual da menos de 2, listamos todas las
            # curvas de la escena para saber dónde terminaron los puntos.
            print("CranioPlan: DIAGNÓSTICO — curvas presentes en la escena:")
            todasLasCurvas = slicer.util.getNodesByClass("vtkMRMLMarkupsCurveNode")
            for c in todasLasCurvas:
                print(f"  - {c.GetName()} (ID {c.GetID()}): {c.GetNumberOfControlPoints()} puntos")

            self.etiquetaEstadoCorte.setText(
                f"Necesitás al menos 2 puntos para trazar el corte (tenés {numeroPuntos}). "
                "No hace falta cerrar la curva ni volver al punto inicial: colocá los "
                "puntos y presioná 'Finalizar trazado' directamente. "
                "Volvé a presionar 'Trazar línea de corte' e intentá de nuevo."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return

        self.botonFinalizarTrazado.enabled = False
        self.botonGenerarCorte.enabled = True
        self.etiquetaEstadoCorte.setText(
            f"Línea de corte lista ({numeroPuntos} puntos). "
            "Ajustá el grosor y presioná 'Generar corte'."
        )
        self.etiquetaEstadoCorte.setStyleSheet("color: #1F4E79;")

    def onBotonGenerarCorteClicked(self):
        if not self._fragmentosActuales:
            self.etiquetaEstadoCorte.setText("Primero confirmá el cráneo (Paso 2).")
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return

        self.etiquetaEstadoCorte.setText("Calculando el corte, por favor esperá...")
        self.etiquetaEstadoCorte.setStyleSheet("color: orange;")
        slicer.app.processEvents()

        grosorMM = self.spinGrosor.value
        volumenActual = self._parameterNode.estudioCargado

        # Cortamos sobre el hueso VIGENTE (los fragmentos que dejó el
        # corte anterior), no sobre el cráneo original. Así los cortes
        # se acumulan en vez de pisarse entre sí.
        fragmentosPrevios = list(self._fragmentosActuales)

        fragmentosNuevos = self.logic.generarOsteotomia(
            fragmentosPrevios,
            self._curvaCorteActual,
            volumenActual,
            grosorMM=grosorMM,
        )

        if not fragmentosNuevos:
            self.etiquetaEstadoCorte.setText(
                "No se pudo calcular el corte. Revisá que la curva esté "
                "bien trazada sobre el hueso, con al menos 2 puntos."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return

        # El corte salió bien: los fragmentos previos ya no representan
        # el hueso actual, así que los sacamos de la escena.
        for viejo in fragmentosPrevios:
            if viejo is not None:
                slicer.mrmlScene.RemoveNode(viejo)

        self._fragmentosActuales = fragmentosNuevos
        self._modeloCraneoActual = fragmentosNuevos[0]  # el mayor: cráneo restante

        extraidos = len(fragmentosNuevos) - 1

        if extraidos >= 1:
            self.etiquetaEstadoCorte.setText(
                f"Corte realizado. Hueso actual: {len(fragmentosNuevos)} pieza(s) — "
                f"el cráneo restante (color hueso) y {extraidos} fragmento(s) "
                "extraído(s), resaltados en color.\n"
                "Podés trazar otro corte sobre el resultado."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: green;")
        else:
            self.etiquetaEstadoCorte.setText(
                "El corte se calculó, pero no separó ninguna pieza nueva. "
                "Si usaste una línea abierta, tené en cuenta que solo separa "
                "si sus extremos llegan a un borde del hueso: para aislar una "
                "región en el medio del cráneo, usá 'Curva cerrada'."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: #B8860B;")

        # Preparar para el próximo corte
        self.botonGenerarCorte.enabled = False
        self._curvaCorteActual = None

        layoutManager = slicer.app.layoutManager()
        layoutManager.threeDWidget(0).threeDView().resetFocalPoint()


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
        BLOQUE A: threshold + islands + filtro de centralidad/tamaño.
        Se DETIENE antes de fusionar y devuelve las islas candidatas
        para que el Widget permita la revisión manual del Bloque B.
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

        segmentEditorWidget.setActiveEffectByName("Threshold")
        thresholdEffect = segmentEditorWidget.activeEffect()
        thresholdEffect.setParameter("MinimumThreshold", "300")
        thresholdEffect.setParameter("MaximumThreshold", "3000")
        thresholdEffect.self().onApply()

        segmentEditorWidget.setActiveEffectByName("Islands")
        islandsEffect = segmentEditorWidget.activeEffect()
        islandsEffect.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
        islandsEffect.self().onApply()

        segmentacion = segmentationNode.GetSegmentation()
        nSegments = segmentacion.GetNumberOfSegments()

        if nSegments == 0:
            segmentEditorWidget = None
            return None

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

        idsCandidatas = set(c[0] for c in candidatasFinales)
        for i in range(nSegments - 1, -1, -1):
            segId = segmentacion.GetNthSegmentID(i)
            if segId not in idsCandidatas:
                segmentacion.RemoveSegment(segId)

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
        BLOQUE B: fusiona las islas restantes en un único segmento
        llamado "Craneo_Final".
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

    # ============================================================
    # BLOQUE F — Planificación de osteotomías (corte propio,
    # sin depender del Osteotomy Planner de KitwareMedical)
    #
    # Decisión de diseño (04/07/2026): se descarta el Curve Cut nativo
    # de Dynamic Modeler (y por extensión el Osteotomy Planner) porque,
    # probado con casos reales del Garrahan:
    #   (a) confunde islas naturalmente desconectadas (suturas
    #       craneales abiertas) con resultado de un corte.
    #   (b) no garantiza atravesar el espesor real del hueso en cada
    #       punto de la curva.
    #
    # Se construye una pared de corte propia, con profundidad y
    # orientación calculadas por rayo desde la normal de superficie en
    # cada punto, y se resta del cráneo con una resta booleana real.
    # Al ser geometría real (no clasificación superficial), una isla
    # que la pared no toca queda intacta.
    #
    # La profundidad se calcula automáticamente por paciente y por
    # punto. El grosor (ancho de la hoja) es un dato clínico y queda
    # expuesto como parámetro editable, pendiente de confirmar con
    # el instrumental real del Garrahan.
    # ============================================================

    def _decimarMalla(self, polyData, reduccion):
        """
        Reduce la cantidad de triángulos de una malla, preservando la
        topología (no abre agujeros ni separa piezas).

        Marching cubes sobre una TC de 0.5 mm genera cientos de miles de
        triángulos: mucho más detalle del que se necesita para planificar
        un corte, y suficiente para que la vista 3D se vuelva lenta al
        rotar. El error geométrico que introduce la decimación queda muy
        por debajo del tamaño de voxel del estudio, así que no afecta la
        precisión de la planificación.

        reduccion: fracción de triángulos a eliminar (0.0 a 1.0).
        0.0 desactiva la decimación.
        """
        if not reduccion or reduccion <= 0.0 or polyData is None:
            return polyData
        if polyData.GetNumberOfCells() == 0:
            return polyData

        decimador = vtk.vtkDecimatePro()
        decimador.SetInputData(polyData)
        decimador.SetTargetReduction(reduccion)
        decimador.PreserveTopologyOn()
        decimador.Update()

        resultado = vtk.vtkPolyData()
        resultado.DeepCopy(decimador.GetOutput())

        if resultado.GetNumberOfPoints() == 0:
            return polyData  # la decimación falló; devolvemos la original
        return resultado

    def exportarSegmentoAModelo(self, segmentationNode, nombreSegmento="Craneo_Final",
                                  reduccionMalla=0.7):
        """
        Exporta un segmento puntual a un vtkMRMLModelNode, sin pasar
        por el módulo Segmentations a mano.

        La conversión de segmento a superficie cerrada suele dejar
        decenas o cientos de fragmentos diminutos de ruido (voxeles
        aislados que pasaron el filtro de tamaño relativo del Bloque A
        por muy poco, sin estar realmente soldados al cráneo). Acá nos
        quedamos únicamente con la pieza conectada más grande, que es
        el cráneo real.

        reduccionMalla (0.0 a 1.0): fracción de triángulos a eliminar.
        Marching cubes sobre una TC de 0.5 mm genera cientos de miles de
        triángulos — mucho más detalle del que se necesita para
        planificar un corte, y suficiente para que la vista 3D se vuelva
        lenta e inusable al rotar el modelo. Con 0.5 se elimina la mitad
        de los triángulos, preservando la topología (PreserveTopologyOn)
        y manteniendo el error geométrico muy por debajo del voxel del
        estudio, así que no afecta la precisión de la planificación.
        Poner 0.0 desactiva la decimación.
        """
        segmentacion = segmentationNode.GetSegmentation()
        segId = segmentacion.GetSegmentIdBySegmentName(nombreSegmento)
        if not segId:
            return None

        segmentationNode.CreateClosedSurfaceRepresentation()
        polyDataBruto = vtk.vtkPolyData()
        segmentationNode.GetClosedSurfaceRepresentation(segId, polyDataBruto)

        if polyDataBruto is None or polyDataBruto.GetNumberOfPoints() == 0:
            return None

        conectividad = vtk.vtkPolyDataConnectivityFilter()
        conectividad.SetInputData(polyDataBruto)
        conectividad.SetExtractionModeToLargestRegion()
        conectividad.Update()

        limpiar = vtk.vtkCleanPolyData()
        limpiar.SetInputConnection(conectividad.GetOutputPort())
        limpiar.Update()

        polyData = limpiar.GetOutput()
        if polyData.GetNumberOfPoints() == 0:
            return None

        trianguloAntes = polyData.GetNumberOfCells()
        polyData = self._decimarMalla(polyData, reduccionMalla)

        print(
            f"CranioPlan: malla del cráneo — {trianguloAntes} triángulos antes de "
            f"decimar, {polyData.GetNumberOfCells()} después."
        )

        modeloNode = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLModelNode', nombreSegmento + '_Modelo'
        )
        modeloNode.SetAndObservePolyData(polyData)
        modeloNode.CreateDefaultDisplayNodes()
        modeloNode.GetDisplayNode().SetColor(0.9, 0.8, 0.6)
        modeloNode.GetDisplayNode().SetScalarVisibility(False)
        return modeloNode

    def _resamplearPuntos(self, puntosOriginales, distanciaMuestreoMM):
        """
        Recorre la polilínea de una curva y devuelve puntos (numpy
        arrays) espaciados uniformemente cada distanciaMuestreoMM.
        Implementación propia, sin depender de utilidades internas
        de Slicer que puedan variar entre versiones.
        """
        n = puntosOriginales.GetNumberOfPoints()
        if n < 2:
            return []

        original = [np.array(puntosOriginales.GetPoint(i)) for i in range(n)]

        resampleados = [original[0]]
        distanciaAcumulada = 0.0
        puntoAnterior = original[0]

        for i in range(1, n):
            puntoActual = original[i]
            segmento = puntoActual - puntoAnterior
            largoSegmento = np.linalg.norm(segmento)

            while distanciaAcumulada + largoSegmento >= distanciaMuestreoMM:
                falta = distanciaMuestreoMM - distanciaAcumulada
                direccion = segmento / largoSegmento if largoSegmento > 1e-9 else segmento
                nuevoPunto = puntoAnterior + direccion * falta
                resampleados.append(nuevoPunto)
                puntoAnterior = nuevoPunto
                segmento = puntoActual - puntoAnterior
                largoSegmento = np.linalg.norm(segmento)
                distanciaAcumulada = 0.0

            distanciaAcumulada += largoSegmento
            puntoAnterior = puntoActual

        if np.linalg.norm(resampleados[-1] - original[-1]) > 1e-6:
            resampleados.append(original[-1])

        return resampleados

    def _medirNormalYEspesorLocal(self, punto, cellLocator, mallaConNormales,
                                    distanciaMaximaBusquedaMM=12.0):
        """
        Devuelve (normalHaciaAfuera, espesorLocalMM) en un punto sobre
        la superficie del cráneo, midiendo el espesor con un rayo
        lanzado hacia adentro desde la normal de superficie.
        """
        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        dist2 = vtk.mutable(0.0)
        puntoMasCercano = [0.0, 0.0, 0.0]

        cellLocator.FindClosestPoint(punto.tolist(), puntoMasCercano, cellId, subId, dist2)

        normales = mallaConNormales.GetCellData().GetNormals()
        if normales is None:
            print("CranioPlan DIAGNÓSTICO: la malla no tiene normales por celda calculadas.")
            return None, None

        numeroCeldas = mallaConNormales.GetNumberOfCells()
        if int(cellId) < 0 or int(cellId) >= numeroCeldas:
            print(
                f"CranioPlan DIAGNÓSTICO: cellId fuera de rango ({int(cellId)} de "
                f"{numeroCeldas}) en punto {punto.tolist()}"
            )
            return None, None

        normal = np.array(normales.GetTuple(int(cellId)))
        normaLongitud = np.linalg.norm(normal)
        if normaLongitud < 1e-9:
            print(
                f"CranioPlan DIAGNÓSTICO: normal degenerada (longitud {normaLongitud}) "
                f"en punto {punto.tolist()}, celda {int(cellId)}"
            )
            return None, None
        normal = normal / normaLongitud

        puntoNp = np.array(puntoMasCercano)
        inicioRayo = puntoNp + normal * 1.0
        finRayo = puntoNp - normal * distanciaMaximaBusquedaMM

        t = vtk.mutable(0.0)
        xInterseccion = [0.0, 0.0, 0.0]
        pcoords = [0.0, 0.0, 0.0]
        subIdRayo = vtk.mutable(0)

        huboInterseccion = cellLocator.IntersectWithLine(
            inicioRayo.tolist(), finRayo.tolist(), 0.01,
            t, xInterseccion, pcoords, subIdRayo
        )

        if huboInterseccion:
            espesor = np.linalg.norm(np.array(xInterseccion) - puntoNp)
        else:
            espesor = distanciaMaximaBusquedaMM

        return normal, espesor

    def _quitarPuntosCoincidentes(self, posiciones, esCerrada, toleranciaMM=0.2):
        """
        Elimina puntos consecutivos que están (casi) en la misma
        posición. Un tramo de longitud cero entre dos puntos rompe el
        cálculo de la tangente (división por cero -> NaN) y arruina
        toda la pared de corte, sin disparar ningún error explícito.

        Aparecen en dos situaciones:
          - El usuario colocó un punto encima de otro (por ejemplo,
            cerrando el lazo a mano cuando la curva ya es cerrada).
          - GetCurvePointsWorld() de una curva CERRADA devuelve el lazo
            completo, repitiendo el punto inicial al final para cerrarlo.
            Ese punto repetido genera un cuadrilátero degenerado al
            construir la pared.
        """
        limpias = []
        for p in posiciones:
            if not limpias or np.linalg.norm(p - limpias[-1]) > toleranciaMM:
                limpias.append(p)

        # En una curva cerrada, el tramo n-1 -> 0 se construye igual, así
        # que si el último punto coincide con el primero sobra.
        if esCerrada and len(limpias) > 2:
            if np.linalg.norm(limpias[-1] - limpias[0]) <= toleranciaMM:
                limpias.pop()

        return limpias

    def _construirParedDeCorte(self, curvaNode, mallaHueso, grosorMM,
                                 margenSeguridadMM, esCerrada,
                                 distanciaMuestreoMM=2.0,
                                 margenExteriorMM=2.0):
        """
        Construye la pared de corte (sólido delgado) que sigue la
        curva de osteotomía, atravesando el espesor real del hueso
        en cada punto, con el grosor de hoja configurado.

        mallaHueso es un vtkPolyData (no un nodo): la unión de todos los
        fragmentos óseos vigentes. Así, tras varios cortes, la pared se
        calcula contra el hueso tal como está en ese momento, no contra
        el cráneo original.

        esCerrada se recibe como PARÁMETRO EXPLÍCITO. Ver la nota en
        generarOsteotomia sobre por qué no se consulta al nodo.
        """
        puntosCurvaOriginal = curvaNode.GetCurvePointsWorld()
        if puntosCurvaOriginal is None or puntosCurvaOriginal.GetNumberOfPoints() < 2:
            return None

        posiciones = self._resamplearPuntos(puntosCurvaOriginal, distanciaMuestreoMM)
        posiciones = self._quitarPuntosCoincidentes(posiciones, esCerrada)
        n = len(posiciones)

        minimoPuntos = 3 if esCerrada else 2
        if n < minimoPuntos:
            print(
                f"CranioPlan DIAGNÓSTICO: solo {n} punto(s) útiles tras resamplear y "
                f"quitar coincidentes (mínimo {minimoPuntos} para una curva "
                f"{'cerrada' if esCerrada else 'abierta'})."
            )
            return None
        print(
            f"CranioPlan DIAGNÓSTICO: curva {'cerrada' if esCerrada else 'abierta'} "
            f"resampleada a {n} puntos útiles."
        )

        normalesFilter = vtk.vtkPolyDataNormals()
        normalesFilter.SetInputData(mallaHueso)
        normalesFilter.ComputeCellNormalsOn()
        normalesFilter.ComputePointNormalsOff()
        normalesFilter.AutoOrientNormalsOn()
        normalesFilter.ConsistencyOn()
        normalesFilter.Update()
        mallaConNormales = normalesFilter.GetOutput()

        cellLocator = vtk.vtkCellLocator()
        cellLocator.SetDataSet(mallaConNormales)
        cellLocator.BuildLocator()

        normales = []
        espesores = []
        for i in range(n):
            normal, espesor = self._medirNormalYEspesorLocal(
                posiciones[i], cellLocator, mallaConNormales
            )
            if normal is None:
                print(
                    f"CranioPlan DIAGNÓSTICO: falló el cálculo de normal/espesor en el "
                    f"punto {i} de {n} (posición {posiciones[i].tolist()}). "
                    "Se aborta la pared de corte."
                )
                return None
            normales.append(normal)
            espesores.append(espesor)

        laterales = []
        for i in range(n):
            if esCerrada:
                # En un lazo, el punto anterior al primero es el último,
                # y el siguiente al último es el primero.
                tangente = posiciones[(i + 1) % n] - posiciones[(i - 1) % n]
            elif i == 0:
                tangente = posiciones[1] - posiciones[0]
            elif i == n - 1:
                tangente = posiciones[n - 1] - posiciones[n - 2]
            else:
                tangente = posiciones[i + 1] - posiciones[i - 1]

            normaTangente = np.linalg.norm(tangente)
            if normaTangente < 1e-9:
                # Dos puntos coincidentes que se escaparon del filtro.
                # Sin esta guarda, la división de abajo produce NaN y
                # toda la pared queda con coordenadas inválidas.
                print(
                    f"CranioPlan DIAGNÓSTICO: tangente degenerada en el punto {i} de {n} "
                    "(puntos coincidentes). Se aborta la pared de corte."
                )
                return None
            tangente = tangente / normaTangente

            lateral = np.cross(tangente, normales[i])
            normaLateral = np.linalg.norm(lateral)
            if normaLateral < 1e-6:
                referencia = np.array([1.0, 0.0, 0.0])
                if abs(np.dot(referencia, normales[i])) > 0.9:
                    referencia = np.array([0.0, 1.0, 0.0])
                lateral = np.cross(referencia, normales[i])
                normaLateral = np.linalg.norm(lateral)
            lateral = lateral / normaLateral
            laterales.append(lateral)

        puntosSolido = vtk.vtkPoints()
        PROFUNDIDAD_MAXIMA_MM = 8.0  # generoso para hueso craneal pediátrico;
        # evita que un resguardo de medición (cuando el rayo no encuentra
        # la tabla interna) dispare una pared desproporcionada respecto
        # al tamaño real de la curva, lo que la haría autointersectarse
        # en lazos chicos y romper la resta booleana.
        for i in range(n):
            p = posiciones[i]
            normal = normales[i]
            lateral = laterales[i]
            profundidad = min(espesores[i] + margenSeguridadMM, PROFUNDIDAD_MAXIMA_MM)
            mitadGrosor = grosorMM / 2.0

            topLeft = p + normal * profundidad + lateral * mitadGrosor
            topRight = p + normal * profundidad - lateral * mitadGrosor
            bottomRight = p - normal * profundidad - lateral * mitadGrosor
            bottomLeft = p - normal * profundidad + lateral * mitadGrosor

            puntosSolido.InsertNextPoint(topLeft.tolist())
            puntosSolido.InsertNextPoint(topRight.tolist())
            puntosSolido.InsertNextPoint(bottomRight.tolist())
            puntosSolido.InsertNextPoint(bottomLeft.tolist())

        triangulos = vtk.vtkCellArray()

        def indice(i, esquina):
            return 4 * i + esquina

        def agregarCuad(a, b, c, d):
            t1 = vtk.vtkTriangle()
            t1.GetPointIds().SetId(0, a)
            t1.GetPointIds().SetId(1, b)
            t1.GetPointIds().SetId(2, c)
            triangulos.InsertNextCell(t1)
            t2 = vtk.vtkTriangle()
            t2.GetPointIds().SetId(0, a)
            t2.GetPointIds().SetId(1, c)
            t2.GetPointIds().SetId(2, d)
            triangulos.InsertNextCell(t2)

        rangoSegmentos = range(n) if esCerrada else range(n - 1)

        for i in rangoSegmentos:
            j = (i + 1) % n
            agregarCuad(indice(i, 0), indice(i, 1), indice(j, 1), indice(j, 0))
            agregarCuad(indice(i, 1), indice(i, 2), indice(j, 2), indice(j, 1))
            agregarCuad(indice(i, 2), indice(i, 3), indice(j, 3), indice(j, 2))
            agregarCuad(indice(i, 3), indice(i, 0), indice(j, 0), indice(j, 3))

        if not esCerrada:
            agregarCuad(indice(0, 0), indice(0, 1), indice(0, 2), indice(0, 3))
            ultimo = n - 1
            agregarCuad(indice(ultimo, 3), indice(ultimo, 2), indice(ultimo, 1), indice(ultimo, 0))

        paredBruta = vtk.vtkPolyData()
        paredBruta.SetPoints(puntosSolido)
        paredBruta.SetPolys(triangulos)

        limpiar = vtk.vtkCleanPolyData()
        limpiar.SetInputData(paredBruta)
        limpiar.Update()

        triangular = vtk.vtkTriangleFilter()
        triangular.SetInputConnection(limpiar.GetOutputPort())
        triangular.Update()

        corregirNormales = vtk.vtkPolyDataNormals()
        corregirNormales.SetInputConnection(triangular.GetOutputPort())
        corregirNormales.ConsistencyOn()
        corregirNormales.AutoOrientNormalsOn()
        corregirNormales.SplittingOff()
        corregirNormales.Update()

        return corregirNormales.GetOutput()

    def generarOsteotomia(self, modelosEntrada, curvaNode, volumeNode,
                            grosorMM=1.0, margenSeguridadMM=3.0,
                            volumenMinimoFragmentoMM3=50.0,
                            reduccionMallaFragmentos=0.7):
        """
        Ejecuta un corte de osteotomía sobre el hueso vigente.

        modelosEntrada: LISTA de vtkMRMLModelNode con los fragmentos
        óseos actuales. En el primer corte es una lista de un solo
        elemento (el cráneo completo). En los cortes siguientes son los
        fragmentos que dejó el corte anterior.

        CORTES ENCADENADOS — corrección (11/07/2026): antes esta función
        recibía siempre el modelo del cráneo ORIGINAL, así que el segundo
        corte re-cortaba el cráneo intacto y descartaba el resultado del
        primero. Ahora recibe los fragmentos vigentes, los fusiona en un
        único volumen de hueso, aplica el corte sobre esa unión y vuelve
        a separar en piezas. Así cada corte se acumula sobre el anterior.

        MOTOR DE CORTE: resta VOLUMÉTRICA (por voxeles) con el Segment
        Editor, no booleano de mallas. El filtro booleano de mallas de
        VTK (vtkBooleanOperationPolyDataFilter) falló sistemáticamente
        sobre mallas craneales reales (errores de vtkDelaunay2D y de
        vtkIntersectionPolyDataFilter). La resta por voxeles no usa
        triangulación de intersecciones, así que no puede fallar por esa
        vía. Contrapartida: la precisión queda limitada al voxel del
        estudio (0.5 mm), más fino que el grosor de hoja que se
        planifica (~1 mm), por lo que no es una limitación práctica.

        Devuelve la lista de vtkMRMLModelNode resultantes (los fragmentos
        óseos tras este corte), o lista vacía si el corte no se pudo
        calcular.
        """
        if not modelosEntrada or curvaNode is None or volumeNode is None:
            print("CranioPlan DIAGNÓSTICO: falta el hueso, la curva o el volumen.")
            return []

        # --- Determinar si la curva es cerrada ---
        # FUENTE DE VERDAD: la clase real del nodo. Ni GetCurveClosed()
        # ni el estado del checkbox resultaron confiables en Slicer 5.10
        # (ambos reportaron "abierta" para un lazo creado y dibujado como
        # vtkMRMLMarkupsClosedCurveNode). IsA() consulta la jerarquía de
        # clases de VTK, así que refleja qué tipo de nodo se creó.
        esCerrada = bool(curvaNode.IsA("vtkMRMLMarkupsClosedCurveNode"))
        print(
            f"CranioPlan DIAGNÓSTICO: curva {curvaNode.GetClassName()} -> "
            f"{'CERRADA' if esCerrada else 'ABIERTA'}. "
            f"Hueso de entrada: {len(modelosEntrada)} fragmento(s)."
        )

        # --- Fusionar los fragmentos vigentes en una sola malla ---
        # Sirve para dos cosas: medir normales/espesor para la pared, y
        # cargar todo el hueso como un único segmento a cortar.
        fusionar = vtk.vtkAppendPolyData()
        for modelo in modelosEntrada:
            if modelo is not None and modelo.GetPolyData() is not None:
                fusionar.AddInputData(modelo.GetPolyData())
        fusionar.Update()
        mallaHueso = fusionar.GetOutput()

        if mallaHueso is None or mallaHueso.GetNumberOfPoints() == 0:
            print("CranioPlan DIAGNÓSTICO: la malla de hueso de entrada está vacía.")
            return []

        paredDeCorte = self._construirParedDeCorte(
            curvaNode, mallaHueso, grosorMM, margenSeguridadMM, esCerrada
        )
        if paredDeCorte is None:
            print("CranioPlan DIAGNÓSTICO: no se pudo construir la pared de corte (ver arriba).")
            return []
        print(
            f"CranioPlan DIAGNÓSTICO: pared de corte con "
            f"{paredDeCorte.GetNumberOfPoints()} puntos y "
            f"{paredDeCorte.GetNumberOfCells()} celdas."
        )

        # --- Segmentación temporal donde hacemos el corte por voxeles ---
        segCorte = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLSegmentationNode', 'CranioPlan_Corte_Temporal'
        )
        segCorte.CreateDefaultDisplayNodes()
        segCorte.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        idHueso = segCorte.AddSegmentFromClosedSurfaceRepresentation(
            mallaHueso, "Hueso", [0.9, 0.8, 0.6]
        )
        idPared = segCorte.AddSegmentFromClosedSurfaceRepresentation(
            paredDeCorte, "ParedDeCorte", [1.0, 0.2, 0.2]
        )

        if not idHueso or not idPared:
            print("CranioPlan DIAGNÓSTICO: no se pudo importar el hueso o la pared a la segmentación.")
            slicer.mrmlScene.RemoveNode(segCorte)
            return []

        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segEditorNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
        segmentEditorWidget.setMRMLSegmentEditorNode(segEditorNode)
        segmentEditorWidget.setSegmentationNode(segCorte)
        segmentEditorWidget.setSourceVolumeNode(volumeNode)

        # --- Resta volumétrica: hueso MENOS pared de corte ---
        segmentEditorWidget.setCurrentSegmentID(idHueso)
        segmentEditorWidget.setActiveEffectByName("Logical operators")
        efectoLogico = segmentEditorWidget.activeEffect()
        efectoLogico.setParameter("Operation", "SUBTRACT")
        efectoLogico.setParameter("ModifierSegmentID", idPared)
        efectoLogico.self().onApply()

        segCorte.GetSegmentation().RemoveSegment(idPared)

        # --- Separar el resultado en piezas conectadas ---
        segmentEditorWidget.setCurrentSegmentID(idHueso)
        segmentEditorWidget.setActiveEffectByName("Islands")
        efectoIslas = segmentEditorWidget.activeEffect()
        efectoIslas.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
        efectoIslas.self().onApply()

        segmentacion = segCorte.GetSegmentation()
        numeroPiezas = segmentacion.GetNumberOfSegments()

        # --- Medir el volumen REAL de cada pieza ---
        espaciado = volumeNode.GetSpacing()
        volumenVoxelMM3 = espaciado[0] * espaciado[1] * espaciado[2]

        piezas = []
        for i in range(numeroPiezas):
            segId = segmentacion.GetNthSegmentID(i)
            arr = slicer.util.arrayFromSegmentBinaryLabelmap(segCorte, segId, volumeNode)
            voxeles = int(np.count_nonzero(arr))
            if voxeles == 0:
                continue
            volumenMM3 = voxeles * volumenVoxelMM3
            piezas.append((volumenMM3, voxeles, segId))

        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segEditorNode)

        if not piezas:
            print("CranioPlan DIAGNÓSTICO: no quedó ninguna pieza con contenido tras el corte.")
            slicer.mrmlScene.RemoveNode(segCorte)
            return []

        piezas.sort(key=lambda p: p[0], reverse=True)

        print(f"CranioPlan DIAGNÓSTICO: tras la resta quedaron {len(piezas)} pieza(s):")
        for volumenMM3, voxeles, segId in piezas:
            estado = "SE CONSERVA" if volumenMM3 >= volumenMinimoFragmentoMM3 else "se descarta (ruido)"
            print(f"    - {volumenMM3 / 1000.0:.2f} cm3 ({voxeles} voxeles) -> {estado}")

        # UMBRAL ABSOLUTO, no relativo. Un umbral relativo al fragmento
        # más grande (p. ej. 1%) descarta por error el disco de una
        # osteotomía, que es legítimamente chico frente al resto del
        # cráneo. Lo que distingue ruido de hueso real es un volumen
        # mínimo absoluto, no su proporción respecto del cráneo entero.
        piezasValidas = [p for p in piezas if p[0] >= volumenMinimoFragmentoMM3]

        if not piezasValidas:
            print("CranioPlan DIAGNÓSTICO: ninguna pieza superó el volumen mínimo.")
            slicer.mrmlScene.RemoveNode(segCorte)
            return []

        # --- Exportar cada fragmento a su propio modelo 3D ---
        segCorte.CreateClosedSurfaceRepresentation()

        # El fragmento MÁS GRANDE conserva el color natural del hueso;
        # los fragmentos extraídos (más chicos) se resaltan con colores
        # distintos, para que se vea de un vistazo qué se cortó.
        COLOR_HUESO = (0.9, 0.8, 0.6)
        COLORES_EXTRAIDOS = [
            (0.55, 0.75, 0.85),  # celeste
            (0.75, 0.85, 0.55),  # verde claro
            (0.85, 0.55, 0.75),  # rosado
            (0.95, 0.75, 0.45),  # naranja suave
        ]

        fragmentos = []
        for idx, (volumenMM3, voxeles, segId) in enumerate(piezasValidas):
            polyDataFragmento = vtk.vtkPolyData()
            segCorte.GetClosedSurfaceRepresentation(segId, polyDataFragmento)
            if polyDataFragmento.GetNumberOfPoints() == 0:
                continue

            copia = vtk.vtkPolyData()
            copia.DeepCopy(polyDataFragmento)

            # Los fragmentos salen de la segmentación a resolución completa
            # (marching cubes). Sin decimarlos, cada corte vuelve a llenar
            # la escena de mallas pesadas y el lag reaparece.
            copia = self._decimarMalla(copia, reduccionMallaFragmentos)

            if idx == 0:
                nombre = "Craneo_restante"
                color = COLOR_HUESO
            else:
                nombre = f"Fragmento_extraido_{idx}"
                color = COLORES_EXTRAIDOS[(idx - 1) % len(COLORES_EXTRAIDOS)]

            nombreUnico = slicer.mrmlScene.GenerateUniqueName(nombre)
            nodo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', nombreUnico)
            nodo.SetAndObservePolyData(copia)
            nodo.CreateDefaultDisplayNodes()
            nodo.GetDisplayNode().SetColor(*color)
            nodo.GetDisplayNode().SetScalarVisibility(False)
            fragmentos.append(nodo)

            print(
                f"CranioPlan: {nombreUnico} — {volumenMM3 / 1000.0:.2f} cm3, "
                f"{copia.GetNumberOfCells()} triángulos."
            )

        slicer.mrmlScene.RemoveNode(segCorte)

        print(f"CranioPlan DIAGNÓSTICO: quedaron {len(fragmentos)} fragmento(s) óseo(s).")
        return fragmentos


#
# CranioPlanTest
#


class CranioPlanTest(ScriptedLoadableModuleTest):

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.delayDisplay("Sin tests automatizados por ahora.")