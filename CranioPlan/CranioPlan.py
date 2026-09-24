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


# Cartel de versión: se imprime en consola al correr los bloques.
# Sirve para confirmar de un vistazo QUÉ versión está realmente cargada
# en Slicer (después de un Reload), y no depender de suponerlo.
CRANIOPLAN_VERSION = "2026-07-26-j (flap por partición geométrica del prisma del lazo: siempre separa)"


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

        # Estado interno del Bloque F (planificación de osteotomías).
        #
        # MODELO MENTAL (correcto): en todo momento hay UN "cráneo
        # restante" (todo el hueso que todavía no se extrajo) y una lista
        # de fragmentos ya extraídos, que quedan aparte. Cada corte opera
        # SOLO sobre el cráneo restante; los fragmentos ya extraídos no se
        # vuelven a cortar. Tras cada corte, el cráneo restante se
        # reemplaza por el nuevo restante y el/los flap(s) recién
        # separados se agregan a la lista de extraídos.
        self._craneoRestante = None        # un vtkMRMLModelNode
        self._fragmentosExtraidos = []     # lista de vtkMRMLModelNode

        self._curvaCorteActual = None
        self._curvaEsCerrada = False  # se fija al trazar, según el checkbox
        self._observadorCurvaTag = None  # para seguir los puntos en tiempo real
        self.placeWidgetCorte = None  # qSlicerMarkupsPlaceWidget, se crea en setup()

    def setup(self) -> None:
        ScriptedLoadableModuleWidget.setup(self)

        print(f"CranioPlan {CRANIOPLAN_VERSION}: módulo cargado.")

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
            "Las marcadas como 'ALEJADA' no tocan el borde del volumen ni están cerca\n"
            "de la masa principal: pueden ser hueso real separado por una sutura\n"
            "abierta, o ruido. Revisalas con atención antes de decidir.\n"
            "Nada se elimina automáticamente salvo la camilla (toca el borde) y el\n"
            "ruido de pocos vóxeles: lo que quede acá, si lo confirmás, se conserva.\n"
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
        self.etiquetaEstadoConfirmacion.setWordWrap(True)
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
        self.etiquetaEstadoCorte.setWordWrap(True)
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
        nAlejadas = sum(1 for isla in self._islasRevision if isla.get("alejada"))
        mensaje = (
            f"Se encontraron {n} isla(s) candidata(s). "
            f"{'Revisalas antes de confirmar.' if n > 1 else 'Una sola isla — podés confirmar directamente.'}"
        )
        if nAlejadas:
            mensaje += (
                f"\n{nAlejadas} de ellas está(n) marcada(s) como ALEJADA: revisalas "
                "con atención, puede ser hueso real separado por una sutura abierta."
            )
        self.etiquetaEstadoCraneo.setText(mensaje)
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
            alejada = isla.get("alejada", False)

            filaWidget = qt.QWidget()
            filaLayout = qt.QHBoxLayout(filaWidget)
            filaLayout.setContentsMargins(0, 2, 0, 2)

            if alejada:
                etiqueta = qt.QLabel(
                    f"Isla {numero}  |  {vol:.1f} cm3  |  {dist:.0f} mm de la masa "
                    "principal — ALEJADA, revisar"
                )
                etiqueta.setStyleSheet("font-size: 10px; color: #B8860B; font-weight: bold;")
            else:
                etiqueta = qt.QLabel(f"Isla {numero}  |  {vol:.1f} cm3  |  dist: {dist:.0f} mm")
                etiqueta.setStyleSheet("font-size: 10px;")
            # addWidget con factor de stretch: en Slicer 5.10 el argumento
            # va POSICIONAL. Pasarlo como kwarg (stretch=2) falla en runtime.
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
        # El cráneo restante arranca siendo el cráneo entero confirmado.
        self._craneoRestante = self.logic.exportarSegmentoAModelo(
            self._segmentationNode, nombreSegmento="Craneo_Final"
        )
        self._fragmentosExtraidos = []

        # Un cráneo pediátrico tiene suturas abiertas, así que el hueso
        # puede venir ya en varias piezas desconectadas antes de cortar
        # nada. Avisarlo evita confundir una pieza preexistente con un
        # fragmento creado por el corte — que era justamente el bug del
        # planificador. Ahora el corte lo tiene en cuenta explícitamente.
        if self._craneoRestante is not None:
            piezasIniciales = self.logic.contarPiezasConectadas(
                self._craneoRestante.GetPolyData()
            )
            if piezasIniciales > 1:
                self.etiquetaEstadoConfirmacion.setText(
                    f"Cráneo final confirmado.\nAviso: el hueso ya viene en "
                    f"{piezasIniciales} piezas desconectadas antes de cortar "
                    "(normal en cráneos pediátricos con suturas abiertas). "
                    "El planificador de cortes lo tiene en cuenta y NO las "
                    "confunde con fragmentos de una osteotomía."
                )
                self.etiquetaEstadoConfirmacion.setStyleSheet(
                    "color: #B8860B; font-weight: bold;"
                )

        # A partir de acá se trabaja sobre el MODELO, no sobre la
        # segmentación. Si dejamos las dos visibles, Slicer renderiza el
        # cráneo dos veces en cada movimiento del mouse — una de las causas
        # del lag al rotar la vista 3D. Ocultar la segmentación no pierde
        # nada: sus datos siguen en la escena.
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
        if self._craneoRestante is None:
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
        # automáticamente. La Logic maneja ambos casos por igual según
        # la clase real del nodo.
        claseNodo = (
            "vtkMRMLMarkupsClosedCurveNode"
            if self.checkCurvaCerrada.checked
            else "vtkMRMLMarkupsCurveNode"
        )
        # Nombre único: GenerateUniqueName agrega sufijo (_1, _2, ...) si el
        # nombre ya existe, así no quedan varios nodos "Osteotomia" con el
        # mismo nombre visible en el panel Data tras varios cortes.
        nombreCurva = slicer.mrmlScene.GenerateUniqueName("Osteotomia")
        curvaNode = slicer.mrmlScene.AddNewNodeByClass(claseNodo, nombreCurva)
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
        self._curvaEsCerrada = bool(self.checkCurvaCerrada.checked)

        self._observadorCurvaTag = curvaNode.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self._onPuntoAgregadoALaCurva
        )

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
        self.placeWidgetCorte.setPlaceModeEnabled(False)
        self._quitarObservadorCurva()

        numeroPuntos = 0 if self._curvaCorteActual is None else self._curvaCorteActual.GetNumberOfControlPoints()

        if numeroPuntos < 2:
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
        if self._craneoRestante is None:
            self.etiquetaEstadoCorte.setText("Primero confirmá el cráneo (Paso 2).")
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return
        if self._curvaCorteActual is None:
            self.etiquetaEstadoCorte.setText(
                "No hay una curva de corte activa. Trazá una línea primero."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            return

        self.etiquetaEstadoCorte.setText("Calculando el corte, por favor esperá...")
        self.etiquetaEstadoCorte.setStyleSheet("color: orange;")
        slicer.app.processEvents()

        grosorMM = self.spinGrosor.value
        volumenActual = self._parameterNode.estudioCargado

        # El corte opera SOLO sobre el cráneo restante actual. Hacemos una
        # copia de respaldo de su geometría por si el corte falla, y otra
        # copia como entrada para la Logic.
        respaldoPD = vtk.vtkPolyData()
        respaldoPD.DeepCopy(self._craneoRestante.GetPolyData())

        mallaRemanente = vtk.vtkPolyData()
        mallaRemanente.DeepCopy(self._craneoRestante.GetPolyData())

        # Quitamos el nodo viejo del restante ANTES de llamar a la Logic:
        # así el nuevo "Craneo_restante" queda con nombre limpio (sin
        # sufijo "_1" feo por colisión de nombres).
        slicer.mrmlScene.RemoveNode(self._craneoRestante)
        self._craneoRestante = None

        indiceInicial = len(self._fragmentosExtraidos) + 1

        resultado = self.logic.generarOsteotomia(
            mallaRemanente,
            self._curvaCorteActual,
            volumenActual,
            indiceInicialFragmento=indiceInicial,
            grosorMM=grosorMM,
        )

        if not resultado or resultado.get("restante") is None:
            # El corte falló: reconstruimos el restante desde el respaldo
            # para no dejar la escena sin cráneo.
            self._craneoRestante = self.logic.crearModeloDesdePolyData(
                respaldoPD, "Craneo_restante", (0.9, 0.8, 0.6)
            )
            self.etiquetaEstadoCorte.setText(
                "No se pudo calcular el corte. El cráneo restante quedó intacto. "
                "Revisá la consola de Python para ver el diagnóstico detallado."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: red;")
            # dejamos habilitado 'Generar corte' por si quiere reintentar
            return

        self._craneoRestante = resultado["restante"]
        nuevosFragmentos = resultado["fragmentos"]
        self._fragmentosExtraidos.extend(nuevosFragmentos)

        piezasCreadas = resultado["piezasCreadas"]
        piezasAntes = resultado["piezasAntes"]

        if piezasCreadas >= 1:
            self.etiquetaEstadoCorte.setText(
                f"Corte realizado. Se extrajo {piezasCreadas} fragmento(s) nuevo(s), "
                "resaltado(s) en color; el resto quedó como cráneo restante (un solo "
                "modelo, aunque tenga varias placas separadas por suturas).\n"
                f"Total de fragmentos extraídos hasta ahora: {len(self._fragmentosExtraidos)}. "
                "Podés trazar otro corte sobre el cráneo restante."
            )
            self.etiquetaEstadoCorte.setStyleSheet("color: green;")
        else:
            self.etiquetaEstadoCorte.setText(
                "El corte se calculó, pero no separó ningún fragmento nuevo del cráneo "
                "restante.\n"
                "Si usaste una línea abierta, recordá que solo separa si sus extremos "
                "llegan a un borde del hueso. Para aislar una región en el medio, usá "
                "'Curva cerrada'."
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

    def _tocaBordeDelVolumen(self, mascara):
        """True si la mascara toca alguna cara del volumen."""
        return bool(
            np.any(mascara[0, :, :]) or np.any(mascara[-1, :, :]) or
            np.any(mascara[:, 0, :]) or np.any(mascara[:, -1, :]) or
            np.any(mascara[:, :, 0]) or np.any(mascara[:, :, -1])
        )

    def generarCandidatas(self, volumeNode, umbralMinimoHU=300, umbralMaximoHU=3000,
                            margenProximidadMM=40.0, umbralRelativo=0.05):
        """
        BLOQUE A: threshold + islands + filtro de borde + filtro de
        proximidad física + filtro de tamaño relativo.

        Esta es la lógica ORIGINAL de Nacho (v5), que dejaba una revisión
        limpia con solo las placas reales del cráneo. Se restaura tal cual
        porque el intento anterior (mandar TODO a revisión, sin filtros)
        llenaba el panel de 24 islas de motitas de 0.1 cm³ y, peor, ensuciaba
        la malla de Craneo_Final: al unir 24 pedazos, la exportación
        generaba cientos de regiones y el Bloque F recibía una malla
        fragmentada imposible de cortar. Filtrar acá mantiene todo limpio.

        Los tres filtros automáticos:
          (a) BORDE: se descarta toda isla que toque el borde del volumen
              (camilla, colchoneta, soportes; el cráneo nunca toca el borde).
          (b) PROXIMIDAD FÍSICA: partiendo de la isla mayor sin tocar borde,
              se agregan iterativamente las que estén a <= margenProximidadMM
              de la masa ya aceptada (transformada de distancia euclídea).
              Esto incorpora los huesos separados por suturas abiertas, que
              están CERCA aunque su centroide quede lejos. Las que quedan más
              lejos que el margen se descartan (ruido, camilla parcial, etc.).
          (c) TAMAÑO RELATIVO: entre las aceptadas, se descartan las menores
              a umbralRelativo del volumen de la mayor.

        NOTA sobre el bug del hueso posterior (26/07/2026): la porción
        occipital que "desaparecía al confirmar" NO se perdía acá — sobrevive
        a estos tres filtros (es una placa grande y cercana). Se perdía en la
        exportación a modelo, por una limpieza de malla con umbral RELATIVO
        que borraba placas chicas legítimas. Ese punto se corrigió aparte,
        en exportarSegmentoAModelo (ahora usa umbral absoluto). Por eso acá se
        puede restaurar la lógica de Nacho sin reintroducir aquel bug.

        Se DETIENE antes de fusionar y devuelve las piezas candidatas para la
        revisión manual del Bloque B.
        """
        print(f"CranioPlan {CRANIOPLAN_VERSION}: Bloque A (generarCandidatas).")

        if volumeNode is None:
            return None

        try:
            from scipy import ndimage
        except ImportError:
            print("CranioPlan: scipy no disponible; el filtro de proximidad lo necesita. "
                  "Se usará solo el filtro de borde y el de tamaño relativo.")
            ndimage = None

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
        thresholdEffect.setParameter("MinimumThreshold", str(umbralMinimoHU))
        thresholdEffect.setParameter("MaximumThreshold", str(umbralMaximoHU))
        thresholdEffect.self().onApply()

        segmentEditorWidget.setActiveEffectByName("Islands")
        islandsEffect = segmentEditorWidget.activeEffect()
        islandsEffect.setParameter("Operation", "SPLIT_ISLANDS_TO_SEGMENTS")
        islandsEffect.self().onApply()

        segmentacion = segmentationNode.GetSegmentation()
        nSegments = segmentacion.GetNumberOfSegments()
        print(f"CranioPlan: threshold {umbralMinimoHU}-{umbralMaximoHU} HU -> {nSegments} isla(s).")

        if nSegments == 0:
            segmentEditorWidget = None
            slicer.mrmlScene.RemoveNode(segmentEditorNode)
            return None

        espaciado = volumeNode.GetSpacing()
        muestreoZYX = (espaciado[2], espaciado[1], espaciado[0])
        volumenVoxelCM3 = (espaciado[0] * espaciado[1] * espaciado[2]) / 1000.0

        islas = {}
        for i in range(nSegments):
            segId = segmentacion.GetNthSegmentID(i)
            mascara = slicer.util.arrayFromSegmentBinaryLabelmap(
                segmentationNode, segId, volumeNode
            ).astype(bool)
            nVoxeles = int(np.count_nonzero(mascara))
            if nVoxeles == 0:
                continue
            islas[segId] = {
                "mask": mascara,
                "vol": nVoxeles * volumenVoxelCM3,
                "tocaBorde": self._tocaBordeDelVolumen(mascara),
            }

        nPorBorde = sum(1 for d in islas.values() if d["tocaBorde"])
        print(
            f"CranioPlan: {nPorBorde} isla(s) descartada(s) por tocar el borde del "
            "volumen (camilla, colchoneta, soportes)."
        )

        candidatasPool = {k: v for k, v in islas.items() if not v["tocaBorde"]}

        if not candidatasPool:
            print(
                "CranioPlan: todas las islas tocan el borde del volumen. "
                "Revisa el campo de vision del estudio o los umbrales HU."
            )
            segmentEditorWidget = None
            slicer.mrmlScene.RemoveNode(segmentEditorNode)
            slicer.mrmlScene.RemoveNode(segmentationNode)
            return None

        refSegId = max(candidatasPool, key=lambda k: candidatasPool[k]["vol"])
        print(
            f"CranioPlan: isla de referencia (mayor volumen sin tocar borde): "
            f"{candidatasPool[refSegId]['vol']:.2f} cm3."
        )

        # --- Filtro de proximidad física (descarta las lejanas) ---
        aceptadas = {refSegId}
        distanciasPorSegId = {refSegId: 0.0}
        mascaraAceptada = candidatasPool[refSegId]["mask"].copy()
        pendientes = {k: v for k, v in candidatasPool.items() if k != refSegId}

        if ndimage is not None and pendientes:
            ronda = 0
            while pendientes:
                ronda += 1
                mapaDistancia = ndimage.distance_transform_edt(
                    ~mascaraAceptada, sampling=muestreoZYX
                )
                nuevas = {}
                for segId, d in pendientes.items():
                    distMin = float(mapaDistancia[d["mask"]].min())
                    if distMin <= margenProximidadMM:
                        nuevas[segId] = distMin
                if not nuevas:
                    break
                for segId, distMin in nuevas.items():
                    mascaraAceptada |= pendientes[segId]["mask"]
                    aceptadas.add(segId)
                    distanciasPorSegId[segId] = distMin
                    del pendientes[segId]
                print(f"CranioPlan:   ronda {ronda}: +{len(nuevas)} isla(s) por proximidad.")
            if pendientes:
                print(
                    f"CranioPlan: {len(pendientes)} isla(s) descartada(s) por estar a más "
                    f"de {margenProximidadMM:.0f} mm de la masa principal (lejos del cráneo)."
                )
        elif pendientes:
            # Sin scipy no se puede medir proximidad; se aceptan todas las
            # no-borde y el filtro de tamaño relativo hace la limpieza gruesa.
            for segId in list(pendientes.keys()):
                aceptadas.add(segId)
                distanciasPorSegId[segId] = 0.0

        # --- Filtro de tamaño relativo (descarta las chicas) ---
        volumenMayor = max(candidatasPool[segId]["vol"] for segId in aceptadas)
        candidatasFinales = [
            segId for segId in aceptadas
            if candidatasPool[segId]["vol"] >= volumenMayor * umbralRelativo
        ]
        nDescartadasTamano = len(aceptadas) - len(candidatasFinales)
        if nDescartadasTamano:
            print(
                f"CranioPlan: {nDescartadasTamano} isla(s) descartada(s) por tamaño "
                f"(< {umbralRelativo:.0%} de la mayor, {volumenMayor:.2f} cm3)."
            )

        idsCandidatas = set(candidatasFinales)
        for i in range(nSegments - 1, -1, -1):
            segId = segmentacion.GetNthSegmentID(i)
            if segId not in idsCandidatas:
                segmentacion.RemoveSegment(segId)

        print(
            f"CranioPlan: {len(candidatasFinales)} pieza(s) candidata(s) "
            "(pasaron borde + proximidad + tamaño). Se detiene para revisión manual."
        )

        islasRevision = []
        coloresOriginales = {}
        ordenadas = sorted(
            candidatasFinales, key=lambda s: candidatasPool[s]["vol"], reverse=True
        )
        for idx, segId in enumerate(ordenadas, start=1):
            vol = candidatasPool[segId]["vol"]
            dist = distanciasPorSegId.get(segId, 0.0)
            seg = segmentacion.GetSegment(segId)
            seg.SetName(f"Pieza_{idx}")
            coloresOriginales[segId] = seg.GetColor()
            islasRevision.append({
                "numero": idx, "segId": segId, "vol": vol, "dist": dist,
                "alejada": False,
            })
            print(f"CranioPlan:   [{idx}] {vol:.2f} cm3 ({dist:.0f} mm de la masa principal).")

        segmentationNode.CreateClosedSurfaceRepresentation()
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

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
    #   (a) confunde islas naturalmente desconectadas (suturas craneales
    #       abiertas) con resultado de un corte.
    #   (b) no garantiza atravesar el espesor real del hueso en cada
    #       punto de la curva.
    #
    # SEGUIMIENTO DE IDENTIDAD (corrección de fondo 26/07/2026): la clave
    # para distinguir "cráneo restante" de "fragmento extraído" NO es el
    # tamaño (ranking por volumen), sino la IDENTIDAD de cada pieza. Ver
    # el docstring de generarOsteotomia.
    # ============================================================

    def contarPiezasConectadas(self, polyData):
        """Cantidad de componentes conectados de una malla."""
        if polyData is None or polyData.GetNumberOfPoints() == 0:
            return 0
        conectividad = vtk.vtkPolyDataConnectivityFilter()
        conectividad.SetInputData(polyData)
        conectividad.SetExtractionModeToAllRegions()
        conectividad.Update()
        return int(conectividad.GetNumberOfExtractedRegions())

    def _limpiarRuidoMalla(self, polyData, minimoPuntos=200):
        """
        Quita SOLO las regiones conectadas que son ruido de marching cubes
        (menos de minimoPuntos puntos), preservando cualquier pieza ósea
        sustancial aunque esté desconectada del resto (huesos separados por
        suturas abiertas).

        Diferencia clave con la versión anterior (que causaba el bug de
        pérdida de hueso): el umbral es ABSOLUTO en cantidad de puntos,
        calibrado al ruido de la conversión segmento->malla, NO relativo al
        tamaño de la pieza mayor. Un umbral relativo (p. ej. 10% de la
        mayor) borra un plato occipital chico legítimo por el solo hecho de
        ser más chico que la bóveda; un umbral absoluto bajo solo elimina
        specks de pocos puntos y deja intacta cualquier placa real.
        """
        if polyData is None or polyData.GetNumberOfPoints() == 0:
            return polyData

        conectividad = vtk.vtkPolyDataConnectivityFilter()
        conectividad.SetInputData(polyData)
        conectividad.SetExtractionModeToAllRegions()
        conectividad.ColorRegionsOn()
        conectividad.Update()
        numeroRegiones = int(conectividad.GetNumberOfExtractedRegions())

        if numeroRegiones <= 1:
            limpiar = vtk.vtkCleanPolyData()
            limpiar.SetInputData(polyData)
            limpiar.Update()
            salida = vtk.vtkPolyData()
            salida.DeepCopy(limpiar.GetOutput())
            return salida

        arrayRegiones = conectividad.GetOutput().GetPointData().GetArray("RegionId")
        if arrayRegiones is None:
            return polyData

        conteo = {}
        for i in range(arrayRegiones.GetNumberOfTuples()):
            rid = int(arrayRegiones.GetTuple1(i))
            conteo[rid] = conteo.get(rid, 0) + 1

        aConservar = [rid for rid, c in conteo.items() if c >= minimoPuntos]
        if not aConservar:
            # Todo cae bajo el umbral: conservamos la mayor para no vaciar.
            aConservar = [max(conteo, key=conteo.get)]

        conectividad.SetExtractionModeToSpecifiedRegions()
        conectividad.InitializeSpecifiedRegionList()
        for rid in aConservar:
            conectividad.AddSpecifiedRegion(rid)
        conectividad.Update()

        limpiar = vtk.vtkCleanPolyData()
        limpiar.SetInputConnection(conectividad.GetOutputPort())
        limpiar.Update()

        salida = vtk.vtkPolyData()
        salida.DeepCopy(limpiar.GetOutput())

        print(
            f"CranioPlan: limpieza de malla — {numeroRegiones} regiones -> "
            f"{len(aConservar)} conservada(s) (umbral absoluto {minimoPuntos} puntos)."
        )
        return salida

    def _decimarMalla(self, polyData, reduccion):
        """
        Reduce la cantidad de triángulos de una malla, preservando la
        topología (no abre agujeros ni separa piezas).

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

    def crearModeloDesdePolyData(self, polyData, nombre, color):
        """
        Crea un vtkMRMLModelNode a partir de un vtkPolyData, con nombre
        único (GenerateUniqueName) y color/visualización estándar.
        """
        if polyData is None or polyData.GetNumberOfPoints() == 0:
            return None
        nombreUnico = slicer.mrmlScene.GenerateUniqueName(nombre)
        nodo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', nombreUnico)
        nodo.SetAndObservePolyData(polyData)
        nodo.CreateDefaultDisplayNodes()
        nodo.GetDisplayNode().SetColor(*color)
        nodo.GetDisplayNode().SetScalarVisibility(False)
        return nodo

    def exportarSegmentoAModelo(self, segmentationNode, nombreSegmento="Craneo_Final",
                                  reduccionMalla=0.7):
        """
        Exporta un segmento a un vtkMRMLModelNode.

        Corrección 26/07/2026: la limpieza de malla ya NO usa un umbral
        RELATIVO (que borraba placas óseas chicas legítimas, causando la
        desaparición del hueso occipital al confirmar). Ahora usa
        _limpiarRuidoMalla con umbral absoluto, que solo quita specks de
        marching cubes y conserva todas las placas reales aunque estén
        desconectadas por suturas abiertas.

        reduccionMalla (0.0 a 1.0): fracción de triángulos a eliminar para
        aligerar el render. 0.0 desactiva la decimación.
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

        polyData = self._limpiarRuidoMalla(polyDataBruto, minimoPuntos=200)
        if polyData is None or polyData.GetNumberOfPoints() == 0:
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
        Recorre la polilínea de una curva y devuelve puntos (numpy arrays)
        espaciados uniformemente cada distanciaMuestreoMM. Implementación
        propia, sin depender de utilidades internas de Slicer.
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

    def _normalPromedioEnPunto(self, punto, pointLocator, normalesArray, radioMM=3.0):
        """
        Normal de superficie promediada sobre los vértices dentro de un
        radio, en vez de tomar la de una sola celda (más estable sobre
        mallas de marching cubes).
        """
        listaIds = vtk.vtkIdList()
        pointLocator.FindPointsWithinRadius(radioMM, punto.tolist(), listaIds)

        if listaIds.GetNumberOfIds() == 0:
            idCercano = pointLocator.FindClosestPoint(punto.tolist())
            if idCercano < 0:
                return None
            listaIds.InsertNextId(idCercano)

        acumulado = np.zeros(3)
        for j in range(listaIds.GetNumberOfIds()):
            acumulado += np.array(normalesArray.GetTuple3(listaIds.GetId(j)))

        norma = np.linalg.norm(acumulado)
        if norma < 1e-6:
            idCercano = pointLocator.FindClosestPoint(punto.tolist())
            if idCercano < 0:
                return None
            normal = np.array(normalesArray.GetTuple3(idCercano))
            norma = np.linalg.norm(normal)
            return normal / norma if norma > 1e-9 else None

        return acumulado / norma

    def _medirEspesorLocal(self, punto, normal, cellLocator, distanciaMaximaMM=12.0):
        """
        Espesor de hueso bajo un punto: lanza un rayo hacia adentro
        siguiendo la normal y mide dónde choca con la tabla interna.
        """
        inicioRayo = punto + normal * 1.0   # 1 mm afuera, evita auto-interseccion
        finRayo = punto - normal * distanciaMaximaMM

        t = vtk.mutable(0.0)
        xInterseccion = [0.0, 0.0, 0.0]
        pcoords = [0.0, 0.0, 0.0]
        subId = vtk.mutable(0)

        hubo = cellLocator.IntersectWithLine(
            inicioRayo.tolist(), finRayo.tolist(), 0.01,
            t, xInterseccion, pcoords, subId
        )
        if hubo:
            return float(np.linalg.norm(np.array(xInterseccion) - punto))
        return distanciaMaximaMM

    def _quitarPuntosCoincidentes(self, posiciones, esCerrada, toleranciaMM=0.2):
        """
        Elimina puntos consecutivos que están (casi) en la misma posición.
        Un tramo de longitud cero rompe el cálculo de la tangente (NaN) y
        arruina la pared de corte sin disparar ningún error.
        """
        limpias = []
        for p in posiciones:
            if not limpias or np.linalg.norm(p - limpias[-1]) > toleranciaMM:
                limpias.append(p)

        if esCerrada and len(limpias) > 2:
            if np.linalg.norm(limpias[-1] - limpias[0]) <= toleranciaMM:
                limpias.pop()

        return limpias

    def _construirParedDeCorte(self, curvaNode, mallaHueso, grosorMM,
                                 margenSeguridadMM, esCerrada,
                                 distanciaMuestreoMM=1.0,
                                 radioNormalMM=3.0):
        """
        Construye la pared de corte (sólido delgado) que sigue la curva de
        osteotomía, atravesando el espesor real del hueso en cada punto, con
        el grosor de hoja configurado.

        mallaHueso es un vtkPolyData: el hueso remanente actual (puede tener
        varias piezas conexas). La pared se calcula contra ese hueso.
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
        normalesFilter.ComputePointNormalsOn()
        normalesFilter.ComputeCellNormalsOff()
        normalesFilter.SplittingOff()
        normalesFilter.ConsistencyOn()
        normalesFilter.AutoOrientNormalsOn()
        normalesFilter.Update()
        mallaConNormales = normalesFilter.GetOutput()

        normalesArray = mallaConNormales.GetPointData().GetNormals()
        if normalesArray is None:
            print("CranioPlan DIAGNOSTICO: no se pudieron calcular normales del hueso.")
            return None

        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(mallaConNormales)
        pointLocator.BuildLocator()

        cellLocator = vtk.vtkCellLocator()
        cellLocator.SetDataSet(mallaConNormales)
        cellLocator.BuildLocator()

        normales = []
        espesores = []
        for i in range(n):
            normal = self._normalPromedioEnPunto(
                posiciones[i], pointLocator, normalesArray, radioNormalMM
            )
            if normal is None:
                print(
                    f"CranioPlan DIAGNOSTICO: no se pudo calcular la normal en el "
                    f"punto {i} de {n}. Se aborta la pared de corte."
                )
                return None
            normales.append(normal)
            espesores.append(
                self._medirEspesorLocal(posiciones[i], normal, cellLocator)
            )

        lateralesCrudos = []
        for i in range(n):
            if esCerrada:
                tangente = posiciones[(i + 1) % n] - posiciones[(i - 1) % n]
            elif i == 0:
                tangente = posiciones[1] - posiciones[0]
            elif i == n - 1:
                tangente = posiciones[n - 1] - posiciones[n - 2]
            else:
                tangente = posiciones[i + 1] - posiciones[i - 1]

            normaTangente = np.linalg.norm(tangente)
            if normaTangente < 1e-9:
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
            lateralesCrudos.append(lateral / normaLateral)

        # Corrección de continuidad de signo (evita el twist de la pared).
        laterales = [lateralesCrudos[0]]
        for i in range(1, n):
            actual = lateralesCrudos[i]
            if float(np.dot(actual, laterales[i - 1])) < 0:
                actual = -actual
            laterales.append(actual)

        if esCerrada and float(np.dot(laterales[-1], laterales[0])) < 0:
            print(
                "CranioPlan DIAGNOSTICO: el lazo tiene un twist impar. La pared puede "
                "quedar irregular en el punto de cierre; revisa que la curva no se cruce."
            )

        puntosSolido = vtk.vtkPoints()
        PROFUNDIDAD_MAXIMA_MM = 8.0
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

    def _anilloDeCorteEnVoxeles(self, boneMask, curvaNode, mallaHueso, volumeNode,
                                 grosorMM, margenSeguridadMM, esCerrada,
                                 distanciaMuestreoMM=0.5, radioNormalMM=3.0):
        """
        Construye el anillo (curva cerrada) o canal (curva abierta) de corte
        DIRECTAMENTE en el espacio de vóxeles, en vez de mallar un tubo fino
        y esperar que la vóxelización lo llene. Por cada punto de la curva se
        barre a lo largo de la normal del hueso (para atravesar todo el
        espesor) y luego se dilata el resultado el grosor de hoja. Así el
        anillo es continuo por construcción: inmune al 'twist' del lazo y a
        los huecos por aliasing, que eran la causa de que el disco no se
        separara.

        Devuelve una máscara booleana (misma forma que boneMask) con los
        vóxeles de HUESO a quitar, o None si no se pudo construir.
        """
        from scipy import ndimage

        puntosCurva = curvaNode.GetCurvePointsWorld()
        if puntosCurva is None or puntosCurva.GetNumberOfPoints() < 2:
            return None

        posiciones = self._resamplearPuntos(puntosCurva, distanciaMuestreoMM)
        posiciones = self._quitarPuntosCoincidentes(posiciones, esCerrada)
        posiciones = list(posiciones)
        if esCerrada and len(posiciones) >= 3:
            posiciones.append(posiciones[0])  # cerrar el lazo explícitamente
        n = len(posiciones)
        if n < 2:
            return None
        print(
            f"CranioPlan DIAGNÓSTICO: curva {'cerrada' if esCerrada else 'abierta'} "
            f"resampleada a {n} puntos para el barrido en vóxeles."
        )

        # Normales del hueso (para barrer a lo largo del espesor).
        normalesFilter = vtk.vtkPolyDataNormals()
        normalesFilter.SetInputData(mallaHueso)
        normalesFilter.ComputePointNormalsOn()
        normalesFilter.ComputeCellNormalsOff()
        normalesFilter.SplittingOff()
        normalesFilter.ConsistencyOn()
        normalesFilter.AutoOrientNormalsOn()
        normalesFilter.Update()
        mallaN = normalesFilter.GetOutput()
        normalesArray = mallaN.GetPointData().GetNormals()
        if normalesArray is None:
            print("CranioPlan DIAGNÓSTICO: no se pudieron calcular normales del hueso.")
            return None
        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(mallaN)
        pointLocator.BuildLocator()

        rasToIjk = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(rasToIjk)

        dims = boneMask.shape  # (z, y, x)
        espaciado = volumeNode.GetSpacing()
        pasoMM = max(min(espaciado) * 0.5, 0.1)
        PROFUNDIDAD_MM = 8.0  # perpendicular al hueso: no ensancha el corte,
                              #  solo garantiza atravesar el espesor completo

        # --- Dirección de extrusión por punto: normal local SUAVIZADA ---
        # Extruir a lo largo de la normal local (perpendicular al hueso)
        # mantiene el corte fino y pegado a la línea trazada. Pero las normales
        # crudas cerca de huecos/zonas finas tienen picos que producían gubias
        # anchas (sobre todo en curvas abiertas). Se suavizan con un promedio
        # móvil a lo largo de la curva: se conserva la curvatura general y se
        # matan los picos. Reemplaza al 'eje global' anterior, que cortaba de
        # más en las zonas inclinadas.
        normalesPunto = []
        for p in posiciones:
            normalesPunto.append(
                self._normalPromedioEnPunto(p, pointLocator, normalesArray, radioNormalMM)
            )

        indicesValidos = [i for i, nrm in enumerate(normalesPunto) if nrm is not None]
        if not indicesValidos:
            print("CranioPlan DIAGNÓSTICO: no se pudo calcular ninguna normal sobre la curva.")
            return None
        # Rellenar faltantes con la normal válida más cercana.
        for i in range(n):
            if normalesPunto[i] is None:
                j = min(indicesValidos, key=lambda k: abs(k - i))
                normalesPunto[i] = np.array(normalesPunto[j])

        # Alinear todas al sentido medio (evita promediar normales opuestas).
        media = np.zeros(3)
        for i in indicesValidos:
            media = media + normalesPunto[i]
        if np.linalg.norm(media) > 1e-9:
            media = media / np.linalg.norm(media)
            for i in range(n):
                if float(np.dot(normalesPunto[i], media)) < 0:
                    normalesPunto[i] = -normalesPunto[i]

        # Promedio móvil a lo largo de la curva (circular si es cerrada).
        ventana = 3  # a cada lado => 7 puntos ~ 3 mm
        normalesSuaves = []
        for i in range(n):
            acum = np.zeros(3)
            for d in range(-ventana, ventana + 1):
                if esCerrada:
                    k = (i + d) % n
                else:
                    k = min(max(i + d, 0), n - 1)
                acum = acum + normalesPunto[k]
            norma = np.linalg.norm(acum)
            normalesSuaves.append(acum / norma if norma > 1e-9 else normalesPunto[i])

        curtain = np.zeros(dims, dtype=bool)

        def marcar(rasXYZ):
            ijk = [0.0, 0.0, 0.0, 0.0]
            rasToIjk.MultiplyPoint([float(rasXYZ[0]), float(rasXYZ[1]), float(rasXYZ[2]), 1.0], ijk)
            ii = int(round(ijk[0]))
            jj = int(round(ijk[1]))
            kk = int(round(ijk[2]))
            if 0 <= kk < dims[0] and 0 <= jj < dims[1] and 0 <= ii < dims[2]:
                curtain[kk, jj, ii] = True

        for i, p in enumerate(posiciones):
            direccion = normalesSuaves[i]
            t = -PROFUNDIDAD_MM
            while t <= PROFUNDIDAD_MM:
                marcar(np.asarray(p, dtype=float) + direccion * t)
                t += pasoMM

        if not curtain.any():
            return None

        # Dilatar el grosor de hoja (kerf). Estructura 3x3x3 (26-conexa): es
        # la que garantiza un anillo CONTINUO sin huecos diagonales sobre la
        # trayectoria curva. La 6-conexa dejaba fugas de 1 vóxel y el disco no
        # se separaba. El ancho lateral queda ~1.5 mm; el corte ya es preciso
        # porque la extrusión es perpendicular al hueso (sin inclinación).
        radioVox = max(1, int(round((grosorMM / 2.0) / min(espaciado))))
        estructura = np.ones((3, 3, 3), dtype=bool)
        curtainDil = ndimage.binary_dilation(curtain, structure=estructura, iterations=radioVox)

        anillo = np.logical_and(curtainDil, boneMask)
        if not anillo.any():
            return None
        return anillo

    def _centroideRASDeMascara(self, mascara, volumeNode):
        """Centroide en coordenadas RAS de una máscara booleana (en la
        geometría del volumeNode). Devuelve np.array([x, y, z]) o None."""
        indices = np.argwhere(mascara)  # columnas z, y, x
        if indices.shape[0] == 0:
            return None
        cen = indices.mean(axis=0)  # [zc, yc, xc]
        ijk = [float(cen[2]), float(cen[1]), float(cen[0]), 1.0]  # x, y, z, 1
        m = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(m)
        ras = [0.0, 0.0, 0.0, 0.0]
        m.MultiplyPoint(ijk, ras)
        return np.array(ras[:3])

    def _puntoDentroDelLazo(self, puntoRAS, puntosLazoRAS):
        """
        True si puntoRAS cae DENTRO del lazo cerrado definido por
        puntosLazoRAS (lista de np.array Nx3), proyectando ambos al plano de
        mejor ajuste del lazo y haciendo un test punto-en-polígono 2D.

        Se usa para identificar, SIN depender del tamaño, qué piezas quedaron
        encerradas por una osteotomía de curva cerrada (el flap) frente a las
        que quedaron afuera (el resto del cráneo). Es robusto aunque el flap
        sea muy chico, que es donde el criterio por volumen fallaba.
        """
        pts = np.asarray(puntosLazoRAS, dtype=float)
        if pts.shape[0] < 3:
            return False
        c0 = pts.mean(axis=0)
        # Los dos vectores singulares de mayor varianza definen el plano del
        # lazo; el tercero es la normal (que no usamos).
        _, _, vh = np.linalg.svd(pts - c0)
        u = vh[0]
        v = vh[1]
        poligono = np.array([[(p - c0).dot(u), (p - c0).dot(v)] for p in pts])
        q = np.array([(puntoRAS - c0).dot(u), (puntoRAS - c0).dot(v)])

        dentro = False
        n = len(poligono)
        j = n - 1
        for i in range(n):
            xi, yi = poligono[i]
            xj, yj = poligono[j]
            if ((yi > q[1]) != (yj > q[1])) and \
               (q[0] < (xj - xi) * (q[1] - yi) / (yj - yi + 1e-12) + xi):
                dentro = not dentro
            j = i
        return dentro

    def _puntosDentroPoligono2D(self, P, poligono):
        """Vectorizado: para P (Nx2) devuelve un bool array indicando si cada
        punto cae dentro del polígono (Mx2), por ray casting."""
        x = P[:, 0]
        y = P[:, 1]
        dentro = np.zeros(len(P), dtype=bool)
        n = len(poligono)
        j = n - 1
        for i in range(n):
            xi, yi = poligono[i]
            xj, yj = poligono[j]
            cond = ((yi > y) != (yj > y)) & \
                   (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-12) + xi)
            dentro = np.logical_xor(dentro, cond)
            j = i
        return dentro

    def _particionPorPrismaDelLazo(self, boneMask, loopRAS, volumeNode, slabMM=25.0):
        """
        Parte el hueso (boneMask) en (flapMask, restanteMask) según el PRISMA
        del lazo cerrado: un vóxel pertenece al flap si su proyección cae
        DENTRO del polígono del lazo (proyectado a su plano de mejor ajuste) Y
        está a menos de ±slabMM del plano. Todo lo demás va al restante.

        Es una partición puramente GEOMÉTRICA: no usa conectividad, así que
        SIEMPRE separa correctamente el hueso encerrado por el lazo, aunque el
        anillo de corte no lo haya desconectado (que es lo que fallaba cerca de
        los agujeros/suturas abiertas, donde quedaban puentes de hueso).

        La franja ±slabMM evita que el prisma (infinito) capture hueso de la
        pared opuesta del cráneo que quede sobre el mismo eje.
        """
        pts = np.asarray(loopRAS, dtype=float)
        if pts.shape[0] < 3:
            return None, None
        c0 = pts.mean(axis=0)
        _, _, vh = np.linalg.svd(pts - c0)
        u = vh[0]
        v = vh[1]
        eje = vh[2]
        poligono = np.column_stack([(pts - c0) @ u, (pts - c0) @ v])

        indices = np.argwhere(boneMask)  # N x 3 (z, y, x)
        if indices.shape[0] == 0:
            return None, None

        M = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(M)
        Mnp = np.array([[M.GetElement(r, c) for c in range(4)] for r in range(4)])
        # IJK homogéneo por vóxel: (i=x, j=y, k=z, 1)
        ijk = np.column_stack([
            indices[:, 2], indices[:, 1], indices[:, 0], np.ones(len(indices))
        ]).astype(float)
        ras = (ijk @ Mnp.T)[:, :3]

        rel = ras - c0
        s = rel @ eje
        a = rel @ u
        b = rel @ v
        enSlab = np.abs(s) <= slabMM
        dentroPoli = self._puntosDentroPoligono2D(np.column_stack([a, b]), poligono)
        esFlap = enSlab & dentroPoli

        flapMask = np.zeros_like(boneMask)
        restMask = np.zeros_like(boneMask)
        zc = indices[:, 0]
        yc = indices[:, 1]
        xc = indices[:, 2]
        flapMask[zc[esFlap], yc[esFlap], xc[esFlap]] = True
        noFlap = ~esFlap
        restMask[zc[noFlap], yc[noFlap], xc[noFlap]] = True
        return flapMask, restMask

    def generarOsteotomia(self, mallaRemanente, curvaNode, volumeNode,
                            indiceInicialFragmento=1,
                            grosorMM=1.0, margenSeguridadMM=3.0,
                            volumenMinimoFragmentoMM3=300.0,
                            reduccionMallaFragmentos=0.7):
        """
        Ejecuta un corte de osteotomía sobre el cráneo remanente.

        mallaRemanente: vtkPolyData con TODO el hueso que todavía no se
        extrajo (el cráneo restante actual). Puede tener varias piezas
        conexas (placas separadas por suturas abiertas). El corte se aplica
        solo sobre este hueso.

        SEGUIMIENTO DE IDENTIDAD — corrección de fondo (26/07/2026):
        --------------------------------------------------------------
        El bug de raíz de todas las versiones anteriores era decidir "qué es
        cráneo restante" y "qué es fragmento extraído" por TAMAÑO: la pieza
        más grande = restante, todas las demás = extraídas. Eso está mal
        porque un cráneo pediátrico ya viene partido en varias placas
        grandes ANTES de cortar (suturas abiertas). El ranking por tamaño
        tomaba esas placas naturales y las llamaba "fragmentos extraídos",
        aunque el corte nunca las tocó.

        La solución correcta es rastrear la IDENTIDAD de cada pieza:
          1. Se voxeliza el hueso de entrada y se etiquetan sus componentes
             conexos ANTES de cortar (scipy.ndimage.label -> "piezas madre").
          2. Se resta la pared de corte (solo QUITA vóxeles) y se re-separan
             las islas -> "piezas hija".
          3. Como restar solo quita vóxeles, cada hija es subconjunto de
             exactamente una madre. Se mapea cada hija a su madre por la
             etiqueta mayoritaria bajo su máscara.
          4. Por cada madre:
               - si produjo <=1 hija real -> el corte NO la separó: va
                 ENTERA al cráneo restante (aunque no sea la más grande).
               - si produjo >=2 hijas reales -> el corte SÍ la separó: la
                 mayor queda en el restante, las otras son flaps extraídos,
                 y las esquirlas (< volumenMinimoFragmentoMM3) vuelven al
                 restante (conserva el volumen óseo).
          5. Craneo_restante = unión (vtkAppendPolyData) de todas las piezas
             del bucket restante, en UN solo modelo (puede tener varias
             placas). Un modelo Fragmento_extraido por flap.

        Así, con un único corte de lazo cerrado sobre una placa, el
        resultado es exactamente: 1 fragmento extraído + el resto del cráneo
        (todas las demás placas + el remanente de la placa cortada), sin
        importar cuántas placas naturales haya ni cuál sea la más grande.

        indiceInicialFragmento: número con el que arranca la numeración de
        los fragmentos de ESTE corte (para que en cortes sucesivos los
        nombres sigan la cuenta: Fragmento_extraido_1, _2, _3, ...).

        volumenMinimoFragmentoMM3: umbral absoluto (300 mm³ = 0.3 cm³ por
        defecto) para distinguir un flap real de una esquirla de la
        vóxelización. Provisional: falta validar con más casos del Garrahan
        que no descarte un fragmento pediátrico legítimamente chico.

        MOTOR DE CORTE: resta volumétrica (por vóxeles) con el Segment
        Editor, no booleano de mallas (vtkBooleanOperationPolyDataFilter
        falla sistemáticamente en mallas craneales reales). Precisión
        limitada al vóxel (0.5 mm), más fino que el grosor de hoja (~1 mm).

        Devuelve un dict:
          {"restante": modelNode, "fragmentos": [modelNode, ...],
           "piezasCreadas": int, "piezasAntes": int}
        o None si el corte no se pudo calcular.
        """
        print(f"CranioPlan {CRANIOPLAN_VERSION}: Bloque F (generarOsteotomia).")

        if mallaRemanente is None or curvaNode is None or volumeNode is None:
            print("CranioPlan DIAGNÓSTICO: falta el hueso, la curva o el volumen.")
            return None
        if mallaRemanente.GetNumberOfPoints() == 0:
            print("CranioPlan DIAGNÓSTICO: la malla de hueso remanente está vacía.")
            return None

        try:
            from scipy import ndimage
        except ImportError:
            print(
                "CranioPlan DIAGNÓSTICO: scipy no está disponible. El seguimiento de "
                "identidad de piezas lo necesita; sin él no se puede distinguir un flap "
                "real de una placa natural. Se aborta el corte."
            )
            return None

        esCerrada = bool(curvaNode.IsA("vtkMRMLMarkupsClosedCurveNode"))
        print(
            f"CranioPlan DIAGNÓSTICO: curva {curvaNode.GetClassName()} -> "
            f"{'CERRADA' if esCerrada else 'ABIERTA'}."
        )

        piezasAntesMalla = self.contarPiezasConectadas(mallaRemanente)
        print(
            f"CranioPlan: hueso remanente de entrada — {piezasAntesMalla} pieza(s) "
            "conexa(s) según la malla."
        )

        # --- Segmentación temporal donde se hace el corte por vóxeles ---
        segCorte = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLSegmentationNode', 'CranioPlan_Corte_Temporal'
        )
        segCorte.CreateDefaultDisplayNodes()
        segCorte.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        idHueso = segCorte.AddSegmentFromClosedSurfaceRepresentation(
            mallaRemanente, "Hueso", [0.9, 0.8, 0.6]
        )
        if not idHueso:
            print("CranioPlan DIAGNÓSTICO: no se pudo importar el hueso a la segmentación.")
            slicer.mrmlScene.RemoveNode(segCorte)
            return None

        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segEditorNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
        segmentEditorWidget.setMRMLSegmentEditorNode(segEditorNode)
        segmentEditorWidget.setSegmentationNode(segCorte)
        segmentEditorWidget.setSourceVolumeNode(volumeNode)

        # === PASO 1 de identidad: etiquetar las piezas madre ANTES de cortar ===
        boneArrAntes = slicer.util.arrayFromSegmentBinaryLabelmap(segCorte, idHueso, volumeNode)
        boneMaskAntes = boneArrAntes.astype(bool)  # astype copia; no es vista viva
        labeledAntes, nMadres = ndimage.label(boneMaskAntes)
        print(f"CranioPlan: {nMadres} pieza(s) madre etiquetada(s) por vóxeles ANTES del corte.")

        # === PASO 2: construir el anillo de corte EN VÓXELES y restarlo ===
        # No se malla un tubo delgado (eso producía el 'twist' del lazo y
        # huecos por aliasing que dejaban el corte incompleto). El anillo se
        # arma barriendo la curva a lo largo de la normal del hueso y
        # dilatándola el grosor de hoja: continuo por construcción.
        anilloCorte = self._anilloDeCorteEnVoxeles(
            boneMaskAntes, curvaNode, mallaRemanente, volumeNode,
            grosorMM, margenSeguridadMM, esCerrada
        )
        if anilloCorte is None:
            print("CranioPlan DIAGNÓSTICO: no se pudo construir el anillo de corte en vóxeles.")
            segmentEditorWidget = None
            slicer.mrmlScene.RemoveNode(segEditorNode)
            slicer.mrmlScene.RemoveNode(segCorte)
            return None

        nVoxAnillo = int(np.count_nonzero(anilloCorte))
        print(
            f"CranioPlan DIAGNÓSTICO: anillo de corte = {nVoxAnillo} vóxeles de hueso a quitar "
            f"(grosor {grosorMM:.1f} mm)."
        )

        boneCut = np.logical_and(boneMaskAntes, np.logical_not(anilloCorte)).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            boneCut, segCorte, idHueso, volumeNode
        )

        # === PASO 3: descomponer el resultado del corte ===
        # Ya no se usa el Segment Editor para separar por islas: la
        # clasificación se hace sobre las máscaras de vóxeles (numpy), más
        # robusto. Para curva CERRADA se parte por geometría (prisma del lazo);
        # para ABIERTA, por identidad de piezas conexas.
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segEditorNode)

        espaciado = volumeNode.GetSpacing()
        volumenVoxelMM3 = espaciado[0] * espaciado[1] * espaciado[2]
        boneCutMask = boneCut.astype(bool)

        flapMasks = []        # una máscara booleana por flap extraído
        restanteMask = None

        if esCerrada:
            # ---- CLASIFICACIÓN GEOMÉTRICA POR EL PRISMA DEL LAZO ----
            # El flap es el hueso ENCERRADO por el lazo (dentro del polígono
            # proyectado y dentro de ±slab del plano). Partición pura por
            # geometría: NO depende de que el anillo haya desconectado el disco
            # por conectividad (que fallaba cerca de agujeros dejando puentes).
            # Por eso SIEMPRE separa el flap del resto.
            puntosLazo = curvaNode.GetCurvePointsWorld()
            loopRAS = None
            if puntosLazo is not None and puntosLazo.GetNumberOfPoints() >= 3:
                loopRAS = np.array([puntosLazo.GetPoint(i)
                                    for i in range(puntosLazo.GetNumberOfPoints())])

            if loopRAS is None:
                print("CranioPlan DIAGNÓSTICO: no pude leer los puntos del lazo; todo al restante.")
                restanteMask = boneCutMask
            else:
                flapMask, restMask = self._particionPorPrismaDelLazo(
                    boneCutMask, loopRAS, volumeNode
                )
                nFlap = int(np.count_nonzero(flapMask)) if flapMask is not None else 0
                nRest = int(np.count_nonzero(restMask)) if restMask is not None else 0
                print(
                    "CranioPlan DIAGNÓSTICO: partición geométrica por el prisma del lazo -> "
                    f"flap {nFlap * volumenVoxelMM3 / 1000.0:.2f} cm3 ({nFlap} vóx), "
                    f"restante {nRest * volumenVoxelMM3 / 1000.0:.2f} cm3."
                )
                UMBRAL_FLAP_VOXELES = 30
                if flapMask is not None and nFlap >= UMBRAL_FLAP_VOXELES:
                    flapMasks.append(flapMask)
                    restanteMask = restMask
                else:
                    print("CranioPlan DIAGNÓSTICO: casi no hay hueso dentro del lazo; 0 flaps.")
                    restanteMask = boneCutMask
        else:
            # ---- CURVA ABIERTA: identidad por conectividad ----
            # Una línea abierta solo separa hueso si cruza una placa de lado a
            # lado. Se etiquetan las piezas conexas tras el corte y se comparan
            # con las madre: si una madre se partió en >=2 piezas reales, la
            # mayor queda en el restante y las otras son flaps.
            labeledCut, nHijas = ndimage.label(boneCutMask)
            restanteMask = np.zeros_like(boneCutMask)
            UMBRAL_RUIDO_VOXELES = 30
            porMadre = {}
            for lab in range(1, nHijas + 1):
                m = (labeledCut == lab)
                vox = int(np.count_nonzero(m))
                if vox < UMBRAL_RUIDO_VOXELES:
                    restanteMask = np.logical_or(restanteMask, m)  # ruido -> restante
                    continue
                etiquetas = labeledAntes[m]
                etiquetas = etiquetas[etiquetas > 0]
                madre = int(np.bincount(etiquetas).argmax()) if etiquetas.size else 0
                porMadre.setdefault(madre, []).append((vox, m))

            print("CranioPlan DIAGNÓSTICO: resultado del corte por pieza madre:")
            for madre in sorted(porMadre.keys()):
                piezas = sorted(porMadre[madre], key=lambda x: x[0], reverse=True)
                reales = [(v, m) for (v, m) in piezas
                          if v * volumenVoxelMM3 >= volumenMinimoFragmentoMM3]
                esquirlas = [(v, m) for (v, m) in piezas
                             if v * volumenVoxelMM3 < volumenMinimoFragmentoMM3]
                if len(reales) <= 1:
                    for (v, m) in piezas:
                        restanteMask = np.logical_or(restanteMask, m)
                    detalle = ", ".join(f"{v * volumenVoxelMM3 / 1000.0:.2f}" for (v, m) in piezas)
                    print(f"    madre {madre}: NO separada ({len(reales)} real; {detalle} cm3) -> restante.")
                else:
                    restanteMask = np.logical_or(restanteMask, reales[0][1])
                    for (v, m) in reales[1:]:
                        flapMasks.append(m)
                    for (v, m) in esquirlas:
                        restanteMask = np.logical_or(restanteMask, m)
                    vols = ", ".join(f"{v * volumenVoxelMM3 / 1000.0:.2f}" for (v, m) in reales)
                    print(
                        f"    madre {madre}: SEPARADA en {len(reales)} reales ({vols} cm3). "
                        f"Mayor al restante; {len(reales) - 1} extraída(s)."
                    )

        # === PASO 4: exportar cada máscara (restante y flaps) a un modelo ===
        def _mascaraAModelo(mascara, nombre, color):
            if mascara is None or not mascara.any():
                return None
            segTmp = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLSegmentationNode', 'CranioPlan_tmp_export'
            )
            segTmp.CreateDefaultDisplayNodes()
            segTmp.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
            sid = segTmp.GetSegmentation().AddEmptySegment("m")
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                mascara.astype(np.uint8), segTmp, sid, volumeNode
            )
            segTmp.CreateClosedSurfaceRepresentation()
            pd = vtk.vtkPolyData()
            segTmp.GetClosedSurfaceRepresentation(sid, pd)
            ok = pd.GetNumberOfPoints() > 0
            pdcopy = vtk.vtkPolyData()
            if ok:
                pdcopy.DeepCopy(pd)
            slicer.mrmlScene.RemoveNode(segTmp)
            if not ok:
                return None
            pdcopy = self._decimarMalla(pdcopy, reduccionMallaFragmentos)
            return self.crearModeloDesdePolyData(pdcopy, nombre, color)

        slicer.mrmlScene.RemoveNode(segCorte)

        restanteModel = _mascaraAModelo(restanteMask, "Craneo_restante", (0.9, 0.8, 0.6))
        if restanteModel is None:
            print("CranioPlan DIAGNÓSTICO: el cráneo restante quedó vacío tras el corte.")
            return None

        COLORES_EXTRAIDOS = [
            (0.55, 0.75, 0.85),  # celeste
            (0.75, 0.85, 0.55),  # verde claro
            (0.85, 0.55, 0.75),  # rosado
            (0.95, 0.75, 0.45),  # naranja suave
        ]
        fragmentosExtraidos = []
        for k, fmask in enumerate(flapMasks):
            indice = indiceInicialFragmento + k
            nombre = f"Fragmento_extraido_{indice}"
            color = COLORES_EXTRAIDOS[(indice - 1) % len(COLORES_EXTRAIDOS)]
            modelo = _mascaraAModelo(fmask, nombre, color)
            if modelo is not None:
                fragmentosExtraidos.append(modelo)
                print(f"CranioPlan: {modelo.GetName()} creado.")

        print(
            f"CranioPlan: corte terminado. {len(fragmentosExtraidos)} fragmento(s) "
            "extraído(s) por este corte."
        )
        print("=" * 60)

        return {
            "restante": restanteModel,
            "fragmentos": fragmentosExtraidos,
            "piezasCreadas": len(fragmentosExtraidos),
            "piezasAntes": nMadres,
        }


#
# CranioPlanTest
#


class CranioPlanTest(ScriptedLoadableModuleTest):

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.delayDisplay("Sin tests automatizados por ahora.")