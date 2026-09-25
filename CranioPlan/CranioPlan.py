# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - MODULO DE 3D SLICER
# Planificacion prequirurgica de craneosinostosis
#
# Trabajo final de grado, Ingenieria Biomedica, FCEFyN-UNC.
# Valentino Andri e Ignacio Gross.
# En colaboracion con el Servicio de Neurocirugia del Hospital Garrahan.
# ============================================================
#
# COMO ESTA ORGANIZADO ESTE ARCHIVO
# ---------------------------------
# Tres clases, como pide Slicer:
#
#   CranioPlan        - registro del modulo (titulo, categoria, ayuda)
#   CranioPlanWidget  - la interfaz: botones, paneles, mensajes al medico
#   CranioPlanLogic   - orquesta los bloques. NO tiene ni un widget adentro.
#
# Los algoritmos NO estan aca: viven en CranioPlanLib/, un archivo por bloque.
# La Logic solo los llama en orden y traduce sus resultados. Esa separacion es
# la que permite que Nacho cambie el corte sin tocar la interfaz y que la
# interfaz cambie sin tocar los algoritmos.
#
#   CranioPlanLib/Comun.py   - nombres de nodos compartidos
#   CranioPlanLib/BloqueC.py - elegir y cargar la serie DICOM correcta
#   CranioPlanLib/BloqueA.py - preparar el craneo (v7.2)
#   CranioPlanLib/BloqueF.py - cortar (v17)
#   CranioPlanLib/BloqueG.py - reacomodar las piezas (v1.1)
#   CranioPlanLib/PuenteFG.py - nombres y atributos de las piezas (F -> G)
#
# EL FLUJO, DE PUNTA A PUNTA
# --------------------------
#   Paso 1  carpeta del paciente  -> Bloque C -> volumen cargado
#   Paso 2  generar craneo        -> Bloque A -> revision -> "Craneo_Final"
#   Paso 3  lineas de corte       -> Bloque F -> piezas separadas
#   Paso 4  reacomodar            -> Bloque G -> armado nuevo + medidas
#   Paso 5  exportar              -> STL para el molde / guia de corte
#
# Cada paso arranca deshabilitado y se habilita cuando el anterior termino. No
# es decoracion: casi todos los errores de la version anterior venian de correr
# un bloque sin que el anterior hubiera dejado lo que ese bloque busca.
# ============================================================

import os
import sys
import time

import qt
import ctk
import slicer
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import parameterNodeWrapper
from slicer import vtkMRMLScalarVolumeNode

# --- La libreria vive al lado de este archivo -------------------------------
# Slicer agrega al sys.path la carpeta del modulo, pero no siempre y no en
# todas las versiones. Agregarla explicitamente evita el "ModuleNotFoundError:
# CranioPlanLib" que aparece solo en algunas instalaciones.
_CARPETA_MODULO = os.path.dirname(os.path.abspath(__file__))
if _CARPETA_MODULO not in sys.path:
    sys.path.insert(0, _CARPETA_MODULO)

import CranioPlanLib                                   # noqa: E402
from CranioPlanLib import Comun                        # noqa: E402
from CranioPlanLib import BloqueC, BloqueA, BloqueF, BloqueG   # noqa: E402
from CranioPlanLib import PuenteFG                          # noqa: E402


# Se imprime en la consola al cargar o recargar el modulo. Actualizar la fecha
# y la version de cada bloque cada vez que se integra una version nueva.
CRANIOPLAN_VERSION = "2026-09-25 - Bloque A v7.2 + Corte v17 + Bloque G v1.1"

# --- Colores de los carteles de estado ---
GRIS    = "color: #777777;"
AZUL    = "color: #1F4E79;"
VERDE   = "color: #1E7B45; font-weight: bold;"
AMBAR   = "color: #B8860B; font-weight: bold;"
ROJO    = "color: #C0392B; font-weight: bold;"
NARANJA = "color: #E67E22;"

ESTILO_BOTON_PRINCIPAL = (
    "background-color: #1E7B45; color: white; font-weight: bold; padding: 7px;")
ESTILO_BOTON_ACCION = (
    "background-color: #1F4E79; color: white; font-weight: bold; padding: 7px;")
ESTILO_BOTON_QUITAR = (
    "padding: 3px 8px; background-color: #C0392B; color: white;")


# ============================================================
# REGISTRO DEL MODULO
# ============================================================
class CranioPlan(ScriptedLoadableModule):

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("CranioPlan")
        self.parent.categories = [translate("qSlicerAbstractCoreModule",
                                            "Craneofacial")]
        self.parent.dependencies = []
        self.parent.contributors = ["Valentino Andri", "Ignacio Gross"]
        self.parent.helpText = _("""
Planificacion prequirurgica de craneosinostosis a partir de la tomografia del
paciente: carga del estudio, reconstruccion 3D del craneo, trazado de las
osteotomias, reacomodamiento de las piezas y exportacion para el molde.

Trabajo final de grado, Ingenieria Biomedica, FCEFyN - UNC.
""")
        self.parent.acknowledgementText = _("""
Desarrollado junto al Servicio de Neurocirugia del Hospital de Pediatria SAMIC
"Prof. Dr. Juan P. Garrahan", Buenos Aires.
""")


# ============================================================
# PARAMETROS QUE SE GUARDAN CON LA ESCENA
# ============================================================
@parameterNodeWrapper
class CranioPlanParameterNode:
    """estudioCargado: el volumen de la tomografia que se cargo en el Paso 1."""
    estudioCargado: vtkMRMLScalarVolumeNode = None


# ============================================================
# INTERFAZ
# ============================================================
class CranioPlanWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):

    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self._parameterNode = None
        self._parameterNodeGuiTag = None

        # Estado de la interfaz
        self._seriesOrdenadas = []      # ranking de series del Paso 1
        self._filasPiezas = []          # filas del panel de revision (Paso 2)
        self._curvas = []               # curvas de corte trazadas (Paso 3)
        self._curvaEnCurso = None
        self._observadorCurva = None
        self._filasFragmentos = []      # filas del panel de piezas (Paso 4)
        self._piezaSeleccionada = None
        self.placeWidget = None

    # --------------------------------------------------------
    # CONSTRUCCION DE LA INTERFAZ
    # --------------------------------------------------------
    def setup(self) -> None:
        ScriptedLoadableModuleWidget.setup(self)

        # Recargar la libreria en cada Reload del modulo. Sin esto, tocar
        # BloqueF.py y apretar Reload deja corriendo el codigo viejo, porque
        # Python cachea los modulos ya importados.
        try:
            CranioPlanLib.recargar()
        except Exception as e:
            print("CranioPlan: no se pudo recargar la libreria: %s" % e)

        self.logic = CranioPlanLogic(log=self._log)
        print("CranioPlan %s: modulo cargado." % CRANIOPLAN_VERSION)

        self._construirPaso1()
        self._construirPaso2()
        self._construirPaso3()
        self._construirPaso4()
        self._construirPaso5()
        self._construirConsola()

        self.layout.addStretch(1)

        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent,
                         self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent,
                         self.onSceneEndClose)
        self.initializeParameterNode()

        self._habilitar(self.paso2, False)
        self._habilitar(self.paso3, False)
        self._habilitar(self.paso4, False)
        self._habilitar(self.paso5, False)

    # ---------- helpers de construccion ----------
    def _titulo(self, texto):
        etiqueta = qt.QLabel(texto)
        etiqueta.setStyleSheet("font-weight: bold; margin-top: 4px;")
        return etiqueta

    def _ayuda(self, texto):
        etiqueta = qt.QLabel(texto)
        etiqueta.setWordWrap(True)
        etiqueta.setStyleSheet("color: #555555;")
        return etiqueta

    def _estado(self, texto):
        etiqueta = qt.QLabel(texto)
        etiqueta.setWordWrap(True)
        etiqueta.setStyleSheet(GRIS)
        return etiqueta

    def _separador(self):
        linea = qt.QFrame()
        linea.setFrameShape(qt.QFrame.HLine)
        linea.setStyleSheet("color: #CCCCCC;")
        return linea

    def _habilitar(self, collapsible, habilitado):
        collapsible.setEnabled(bool(habilitado))
        if habilitado:
            collapsible.collapsed = False

    def _poner(self, etiqueta, texto, estilo=GRIS):
        etiqueta.setText(texto)
        etiqueta.setStyleSheet(estilo)
        slicer.app.processEvents()

    def _esperando(self, etiqueta, texto):
        self._poner(etiqueta, texto, NARANJA)
        slicer.app.setOverrideCursor(qt.Qt.WaitCursor)
        slicer.app.processEvents()

    def _listo(self):
        slicer.app.restoreOverrideCursor()

    # ========================================================
    # PASO 1 - CARGAR LA TOMOGRAFIA
    # ========================================================
    def _construirPaso1(self):
        self.paso1 = ctk.ctkCollapsibleButton()
        self.paso1.text = "Paso 1 - Cargar la tomografia del paciente"
        self.layout.addWidget(self.paso1)
        caja = qt.QVBoxLayout(self.paso1)

        caja.addWidget(self._ayuda(
            "Elegí la carpeta del estudio o el .zip tal cual lo bajaste del "
            "portal del Garrahan.\n"
            "El sistema mira todas las series y propone la más apta para "
            "reconstruir el hueso. Si hiciera falta, la podés cambiar en la "
            "lista de abajo."))

        fila = qt.QHBoxLayout()
        self.botonCarpeta = qt.QPushButton("Elegir carpeta...")
        self.botonCarpeta.setStyleSheet(ESTILO_BOTON_ACCION)
        self.botonZip = qt.QPushButton("Elegir .zip...")
        self.botonZip.setStyleSheet(ESTILO_BOTON_ACCION)
        fila.addWidget(self.botonCarpeta)
        fila.addWidget(self.botonZip)
        caja.addLayout(fila)

        self.etiquetaSeries = self._titulo("Serie que se va a usar:")
        caja.addWidget(self.etiquetaSeries)
        self.etiquetaSeries.setVisible(False)

        self.comboSeries = qt.QComboBox()
        self.comboSeries.toolTip = (
            "Series aptas para segmentar hueso, de mejor a peor.\n"
            "La primera, marcada con *, es la que el sistema recomienda.")
        caja.addWidget(self.comboSeries)
        self.comboSeries.setVisible(False)

        self.botonCargarSerie = qt.QPushButton("Cargar la serie elegida")
        self.botonCargarSerie.setStyleSheet(ESTILO_BOTON_PRINCIPAL)
        caja.addWidget(self.botonCargarSerie)
        self.botonCargarSerie.setVisible(False)

        self.estadoCarga = self._estado("Todavia no se cargo ningun estudio.")
        caja.addWidget(self.estadoCarga)

        # Las series descartadas y el motivo de cada una. Colapsado: es
        # informacion de control, no algo que el medico tenga que leer siempre.
        # Pero cuando el sistema elige mal, esto es lo primero que hay que
        # mirar, y tenerlo a un click evita ir a buscar la consola de Python.
        self.detalleSeries = ctk.ctkCollapsibleGroupBox()
        self.detalleSeries.title = "Ver todas las series del estudio"
        self.detalleSeries.collapsed = True
        cajaDetalle = qt.QVBoxLayout(self.detalleSeries)
        self.textoSeries = qt.QTextEdit()
        self.textoSeries.setReadOnly(True)
        self.textoSeries.setMinimumHeight(130)
        self.textoSeries.setStyleSheet("font-family: monospace; font-size: 11px;")
        cajaDetalle.addWidget(self.textoSeries)
        caja.addWidget(self.detalleSeries)
        self.detalleSeries.setVisible(False)

        self.botonCarpeta.connect("clicked(bool)", self.onSeleccionarCarpeta)
        self.botonZip.connect("clicked(bool)", self.onSeleccionarZip)
        self.botonCargarSerie.connect("clicked(bool)", self.onCargarSerie)

    def onSeleccionarCarpeta(self):
        ruta = qt.QFileDialog.getExistingDirectory(
            self.parent, "Elegí la carpeta del estudio DICOM")
        if ruta:
            self._analizarEstudio(ruta)

    def onSeleccionarZip(self):
        ruta = qt.QFileDialog.getOpenFileName(
            self.parent, "Elegí el .zip del estudio", "", "Archivos ZIP (*.zip)")
        if ruta:
            self._analizarEstudio(ruta)

    def _analizarEstudio(self, ruta):
        self._esperando(self.estadoCarga,
                        "Leyendo el estudio. Con un .zip grande esto tarda un "
                        "rato...")
        try:
            resultado = self.logic.analizarEstudio(ruta)
        finally:
            self._listo()

        if not resultado or not resultado.get("ordenadas"):
            self.comboSeries.setVisible(False)
            self.etiquetaSeries.setVisible(False)
            self.botonCargarSerie.setVisible(False)
            if resultado:
                self._mostrarListaDeSeries(resultado)
            self._poner(self.estadoCarga,
                        "No encontré ninguna serie volumétrica apta en esa "
                        "carpeta.\nAbrí 'Ver todas las series del estudio': ahí "
                        "está el motivo por el que se descartó cada una.", ROJO)
            return

        self._seriesOrdenadas = resultado["ordenadas"]
        self.comboSeries.clear()
        for i, s in enumerate(self._seriesOrdenadas):
            self.comboSeries.addItem(BloqueC.etiquetaDeSerie(s, recomendada=(i == 0)))
        self.comboSeries.setCurrentIndex(0)

        self.etiquetaSeries.setVisible(True)
        self.comboSeries.setVisible(True)
        self.botonCargarSerie.setVisible(True)
        self._mostrarListaDeSeries(resultado)

        primera = self._seriesOrdenadas[0]
        nombre = primera["series_desc"] or ("n#%s" % primera["series_num"])
        empatadas = resultado.get("empatadas") or []

        if empatadas:
            # El algoritmo no puede desempatar dos reconstrucciones del mismo
            # protocolo. En vez de elegir una callado, lo dice.
            self._poner(self.estadoCarga,
                        "Hay %d series con el mismo nombre y casi la misma "
                        "cantidad de cortes. No puedo saber cuál querés: "
                        "elegí vos en la lista de arriba y presíoná 'Cargar la "
                        "serie elegida'." % len(empatadas), AMBAR)
            return

        aviso = ""
        if resultado["nEstudios"] > 1:
            aviso += ("\nOJO: la carpeta tiene %d estudios de fechas distintas. "
                      "Verificá en la lista que la serie sea la del estudio "
                      "correcto." % resultado["nEstudios"])
        if resultado["nPacientes"] > 1:
            aviso += ("\nOJO: hay %d pacientes distintos en la carpeta."
                      % resultado["nPacientes"])

        self._poner(self.estadoCarga,
                    "%d serie(s) apta(s). Recomendada: %s\n(%s)%s\n"
                    "Revisá y presíoná 'Cargar la serie elegida'."
                    % (len(self._seriesOrdenadas), nombre, primera["razon"], aviso),
                    AMBAR if aviso else AZUL)

    def onCargarSerie(self):
        idx = self.comboSeries.currentIndex
        if idx < 0 or idx >= len(self._seriesOrdenadas):
            return
        serie = self._seriesOrdenadas[idx]

        self._esperando(self.estadoCarga, "Cargando la serie, esperá...")
        try:
            volumen = self.logic.cargarSerie(serie["series_uid"])
        finally:
            self._listo()

        if volumen is None:
            self._poner(self.estadoCarga,
                        "No se pudo cargar esa serie. Mirá el detalle técnico "
                        "de abajo.", ROJO)
            return

        self._parameterNode.estudioCargado = volumen
        self._poner(self.estadoCarga,
                    "Cargado: %s\n(%s)" % (volumen.GetName(), serie["razon"]),
                    VERDE if serie.get("esHueso") else AMBAR)

        self._habilitar(self.paso2, True)
        self._poner(self.estadoCraneo, "Todo listo para generar el craneo.", GRIS)

    def _mostrarListaDeSeries(self, resultado):
        """Tabla de control: que series había y por que se descarto cada una."""
        lineas = ["%-4s %-38s %6s %8s  %s"
                  % ("n#", "SERIE", "CORTES", "PASO mm", "ESTADO")]
        for i, s in enumerate(resultado.get("ordenadas", [])):
            sp = "%.2f" % s["spacing_z_real"] if s.get("spacing_z_real") else "?"
            estado = ("ELEGIDA - " + s["razon"]) if i == 0 else s["razon"]
            lineas.append("%-4s %-38s %6d %8s  %s"
                          % (s.get("series_num", "?"),
                             (s.get("series_desc") or "")[:38],
                             s.get("n_imagenes", 0), sp, estado))
        descartadas = resultado.get("descartadas", [])
        if descartadas:
            lineas.append("")
            lineas.append("DESCARTADAS:")
            for s in descartadas:
                lineas.append("%-4s %-38s %6d %8s  %s"
                              % (s.get("series_num", "?"),
                                 (s.get("series_desc") or "")[:38],
                                 s.get("n_imagenes", 0), "-",
                                 s.get("motivo", "")))
        self.textoSeries.setPlainText("\n".join(lineas))
        self.detalleSeries.setVisible(True)


    # ========================================================
    # PASO 2 - GENERAR EL CRANEO 3D
    # ========================================================
    def _construirPaso2(self):
        self.paso2 = ctk.ctkCollapsibleButton()
        self.paso2.text = "Paso 2 - Generar el craneo en 3D"
        self.layout.addWidget(self.paso2)
        caja = qt.QVBoxLayout(self.paso2)

        caja.addWidget(self._ayuda(
            "El sistema separa el hueso del resto de la imagen y arma el modelo "
            "3D. Despues te muestra las piezas que encontro para que revises si "
            "alguna no corresponde (una vertebra, el chupete, la camilla)."))

        self.checkPostop = qt.QCheckBox(
            "El paciente ya fue operado (tomografia postoperatoria)")
        self.checkPostop.toolTip = (
            "En un craneo ya operado los colgajos estan separados varios "
            "milimetros. Con esta opcion tildada el sistema los acepta igual, "
            "pero entran tambien las vertebras: hay que revisar todo.")
        caja.addWidget(self.checkPostop)

        self.botonGenerar = qt.QPushButton("Generar el craneo 3D")
        self.botonGenerar.setStyleSheet(ESTILO_BOTON_ACCION)
        caja.addWidget(self.botonGenerar)

        self.estadoCraneo = self._estado("Todavia no se genero el craneo.")
        caja.addWidget(self.estadoCraneo)

        # --- panel de revision ---
        self.panelRevision = qt.QWidget()
        cajaRev = qt.QVBoxLayout(self.panelRevision)
        cajaRev.setContentsMargins(0, 6, 0, 0)
        cajaRev.addWidget(self._separador())
        cajaRev.addWidget(self._titulo("Revisá las piezas encontradas"))
        cajaRev.addWidget(self._ayuda(
            "Verde = el sistema esta seguro de que es craneo.\n"
            "Naranja = no esta seguro: miralas con 'Ver' y quitá las que no "
            "sean hueso del craneo (las vertebras del cuello son las mas "
            "frecuentes).\n"
            "Si la pieza principal esta mal elegida, marcá otra con 'Es la "
            "principal'."))

        filaBotones = qt.QHBoxLayout()
        self.botonSoloDudosas = qt.QPushButton("Ver solo las dudosas")
        self.botonVerTodas = qt.QPushButton("Ver todas")
        self.botonQuitarDudosas = qt.QPushButton("Quitar todas las dudosas")
        self.botonQuitarDudosas.setStyleSheet(ESTILO_BOTON_QUITAR)
        for b in (self.botonSoloDudosas, self.botonVerTodas, self.botonQuitarDudosas):
            filaBotones.addWidget(b)
        cajaRev.addLayout(filaBotones)

        self.contenedorPiezas = qt.QWidget()
        self.layoutPiezas = qt.QVBoxLayout(self.contenedorPiezas)
        self.layoutPiezas.setContentsMargins(0, 0, 0, 0)
        cajaRev.addWidget(self.contenedorPiezas)

        self.botonConfirmar = qt.QPushButton(
            "Confirmar el craneo y crear el modelo 3D")
        self.botonConfirmar.setStyleSheet(ESTILO_BOTON_PRINCIPAL)
        self.botonConfirmar.toolTip = (
            "Une las piezas que quedaron en un solo craneo y genera la malla 3D "
            "sobre la que se van a trazar los cortes")
        cajaRev.addWidget(self.botonConfirmar)

        caja.addWidget(self.panelRevision)
        self.panelRevision.setVisible(False)

        self.estadoConfirmacion = self._estado("")
        caja.addWidget(self.estadoConfirmacion)

        self.botonGenerar.connect("clicked(bool)", self.onGenerarCraneo)
        self.checkPostop.connect("toggled(bool)", self.onCambiarModoPostop)
        self.botonSoloDudosas.connect("clicked(bool)", self.onVerSoloDudosas)
        self.botonVerTodas.connect("clicked(bool)", self.onVerTodas)
        self.botonQuitarDudosas.connect("clicked(bool)", self.onQuitarDudosas)
        self.botonConfirmar.connect("clicked(bool)", self.onConfirmarCraneo)

    def onGenerarCraneo(self):
        volumen = self._parameterNode.estudioCargado
        if volumen is None:
            self._poner(self.estadoCraneo,
                        "Primero hay que cargar la tomografia (Paso 1).", ROJO)
            return

        self._esperando(self.estadoCraneo,
                        "Generando el craneo. Esto tarda entre 30 segundos y "
                        "dos minutos segun el tamano del estudio...")
        try:
            piezas = self.logic.generarCraneo(volumen,
                                              self.checkPostop.checked)
        finally:
            self._listo()

        if not piezas:
            self._poner(self.estadoCraneo,
                        "No se pudo generar el craneo con esta serie. Fijate de "
                        "haber cargado la serie correcta.", ROJO)
            return

        self._construirPanelPiezas(piezas)
        nDudosas = sum(1 for p in piezas if p["zona"] == "DUDOSA")
        avisos = self.logic.bloqueA.avisos

        texto = "Se encontraron %d pieza(s): %d segura(s) y %d dudosa(s)." % (
            len(piezas), len(piezas) - nDudosas, nDudosas)
        if avisos:
            texto += "\n\n" + "\n".join("- " + a for a in avisos)
        self._poner(self.estadoCraneo, texto, AMBAR if avisos else AZUL)
        self.panelRevision.setVisible(True)

    def onCambiarModoPostop(self, activado):
        """Si ya se genero el craneo, cambiar el modo rehace la revision sin
        volver a segmentar (la parte lenta ya esta hecha). Sin esto, tildar la
        casilla despues de generar no hacia nada y el medico creia que si."""
        if not self.logic.bloqueA.filas:
            return
        self._esperando(self.estadoCraneo, "Rehaciendo la revision...")
        try:
            piezas = self.logic.bloqueA.redecidir(modoPostop=bool(activado))
        finally:
            self._listo()
        if not piezas:
            return
        self._construirPanelPiezas(piezas)
        nDudosas = sum(1 for p in piezas if p["zona"] == "DUDOSA")
        self._poner(self.estadoCraneo,
                    "Revision rehecha en modo %s: %d pieza(s), %d dudosa(s).%s"
                    % ("POSTOPERATORIO" if activado else "preoperatorio",
                       len(piezas), nDudosas,
                       "\nEn postoperatorio entran tambien las vertebras: hay "
                       "que revisar todo antes de confirmar." if activado else ""),
                    AMBAR if activado else AZUL)

    def _construirPanelPiezas(self, piezas):
        while self.layoutPiezas.count():
            item = self.layoutPiezas.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._filasPiezas = []

        for p in piezas:
            fila = qt.QWidget()
            cajaFila = qt.QHBoxLayout(fila)
            cajaFila.setContentsMargins(0, 2, 0, 2)

            esDudosa = p["zona"] == "DUDOSA"
            texto = "Pieza %d   %.1f cm3" % (p["numero"], p["vol"])
            if p["referencia"]:
                texto += "   - PIEZA PRINCIPAL"
            elif p["distRef"] is not None:
                texto += "   a %.1f mm del craneo" % p["distRef"]
            if esDudosa and p["nota"]:
                texto += "\n      %s" % p["nota"]

            etiqueta = qt.QLabel(texto)
            etiqueta.setWordWrap(True)
            if p["referencia"]:
                etiqueta.setStyleSheet("font-size: 11px; font-weight: bold;")
            elif esDudosa:
                etiqueta.setStyleSheet("font-size: 11px; color: #B8860B;")
            else:
                etiqueta.setStyleSheet("font-size: 11px;")
            # OJO: en Slicer 5.10 el factor de stretch va POSICIONAL.
            # Pasarlo como kwarg (stretch=3) falla en runtime.
            cajaFila.addWidget(etiqueta, 3)

            botonVer = qt.QPushButton("Ver")
            botonVer.setFixedWidth(55)
            botonVer.toolTip = "Pinta esta pieza de rojo en la vista 3D"
            cajaFila.addWidget(botonVer)

            botonPrincipal = qt.QPushButton("Es la principal")
            botonPrincipal.setFixedWidth(105)
            botonPrincipal.toolTip = (
                "Rehace la revision tomando esta pieza como el craneo. Usalo si "
                "el sistema eligio mal.")
            botonPrincipal.setEnabled(not p["referencia"])
            cajaFila.addWidget(botonPrincipal)

            botonQuitar = qt.QPushButton("Quitar")
            botonQuitar.setFixedWidth(60)
            botonQuitar.setStyleSheet(ESTILO_BOTON_QUITAR)
            botonQuitar.setEnabled(not p["referencia"])
            cajaFila.addWidget(botonQuitar)

            self.layoutPiezas.addWidget(fila)
            self._filasPiezas.append({"numero": p["numero"], "widget": fila})

            botonVer.connect("clicked(bool)",
                             lambda _c, n=p["numero"]: self.onVerPieza(n))
            botonPrincipal.connect("clicked(bool)",
                                   lambda _c, n=p["numero"]: self.onHacerPrincipal(n))
            botonQuitar.connect("clicked(bool)",
                                lambda _c, n=p["numero"]: self.onQuitarPieza(n))

    def onVerPieza(self, numero):
        self.logic.bloqueA.resaltar(numero)

    def onVerTodas(self):
        self.logic.bloqueA.mostrarTodas()

    def onVerSoloDudosas(self):
        n = self.logic.bloqueA.soloDudosas()
        self._poner(self.estadoCraneo,
                    "Quedaron a la vista %d pieza(s) dudosa(s)." % n, AZUL)

    def onQuitarPieza(self, numero):
        if not self.logic.bloqueA.eliminar(numero):
            return
        fila = next((f for f in self._filasPiezas if f["numero"] == numero), None)
        if fila:
            fila["widget"].setVisible(False)
            self._filasPiezas = [f for f in self._filasPiezas
                                 if f["numero"] != numero]
        self._poner(self.estadoCraneo,
                    "Pieza %d quitada. Quedan %d."
                    % (numero, len(self.logic.bloqueA.piezas)), AZUL)

    def onQuitarDudosas(self):
        n = self.logic.bloqueA.eliminarDudosas()
        self._construirPanelPiezas(self.logic.bloqueA.piezas)
        self._poner(self.estadoCraneo,
                    "Se quitaron %d pieza(s) dudosa(s)." % n, AZUL)

    def onHacerPrincipal(self, numero):
        self._esperando(self.estadoCraneo, "Rehaciendo la revision...")
        try:
            piezas = self.logic.cambiarPiezaPrincipal(numero)
        finally:
            self._listo()
        if not piezas:
            self._poner(self.estadoCraneo,
                        "No se pudo usar esa pieza como principal.", ROJO)
            return
        self._construirPanelPiezas(piezas)
        nDudosas = sum(1 for p in piezas if p["zona"] == "DUDOSA")
        self._poner(self.estadoCraneo,
                    "Revision rehecha con la pieza %d como craneo: %d pieza(s), "
                    "%d dudosa(s)." % (numero, len(piezas), nDudosas), AZUL)

    def onConfirmarCraneo(self):
        if not self.logic.bloqueA.piezas:
            self._poner(self.estadoConfirmacion,
                        "No queda ninguna pieza. Volvé a generar el craneo.", ROJO)
            return

        self._esperando(self.estadoConfirmacion,
                        "Uniendo las piezas y armando el modelo 3D...")
        try:
            resultado = self.logic.confirmarCraneo()
        finally:
            self._listo()

        if resultado is None:
            self._poner(self.estadoConfirmacion,
                        "No se pudo confirmar el craneo. Mirá el detalle de "
                        "abajo.", ROJO)
            return

        texto = ("Craneo confirmado: %.1f cm3 de hueso, en %d pieza(s) "
                 "separadas.\nYa podés trazar los cortes."
                 % (resultado["volumenCM3"], resultado["piezasConexas"]))
        if resultado["piezasConexas"] > 1:
            texto += ("\nQue venga en varias piezas es normal en un craneo "
                      "pediatrico: son las suturas todavia abiertas. El corte "
                      "las tiene en cuenta y no las confunde con fragmentos.")
        self._poner(self.estadoConfirmacion, texto, VERDE)

        self.panelRevision.setVisible(False)
        self._habilitar(self.paso3, True)
        self.paso2.collapsed = True
        self._resetearVista3D()

    # ========================================================
    # PASO 3 - TRAZAR Y HACER LOS CORTES
    # ========================================================
    def _construirPaso3(self):
        self.paso3 = ctk.ctkCollapsibleButton()
        self.paso3.text = "Paso 3 - Planificar y hacer los cortes"
        self.layout.addWidget(self.paso3)
        caja = qt.QVBoxLayout(self.paso3)

        caja.addWidget(self._ayuda(
            "Trazá cada osteotomia como una linea sobre el craneo, en la vista "
            "3D: hacé click en el punto donde empieza el corte, seguí con "
            "clicks a lo largo del recorrido, y cuando llegues al final apretá "
            "'Terminar esta linea'.\n"
            "El corte entra siempre perpendicular al hueso, con la profundidad "
            "justa para atravesarlo. Podés trazar todas las lineas que "
            "necesites antes de cortar."))

        self.botonTrazar = qt.QPushButton("Trazar una linea de corte")
        self.botonTrazar.setStyleSheet(ESTILO_BOTON_ACCION)
        caja.addWidget(self.botonTrazar)

        self.botonTerminarLinea = qt.QPushButton("Terminar esta linea")
        self.botonTerminarLinea.setEnabled(False)
        caja.addWidget(self.botonTerminarLinea)

        self.estadoTrazado = self._estado("Todavia no trazaste ninguna linea.")
        caja.addWidget(self.estadoTrazado)

        self.contenedorLineas = qt.QWidget()
        self.layoutLineas = qt.QVBoxLayout(self.contenedorLineas)
        self.layoutLineas.setContentsMargins(0, 0, 0, 0)
        caja.addWidget(self.contenedorLineas)

        caja.addWidget(self._separador())

        filaGrosor = qt.QHBoxLayout()
        filaGrosor.addWidget(qt.QLabel("Ancho de la sierra (mm):"))
        self.spinGrosor = qt.QDoubleSpinBox()
        self.spinGrosor.setRange(0.3, 5.0)
        self.spinGrosor.setSingleStep(0.1)
        self.spinGrosor.setValue(BloqueF.GROSOR_CORTE_MM)
        self.spinGrosor.toolTip = (
            "Ancho real de la ranura que deja la sierra en el hueso. El valor "
            "por defecto es 1.2 mm; confirmalo con el equipo del Garrahan.")
        filaGrosor.addWidget(self.spinGrosor)
        filaGrosor.addStretch(1)
        caja.addLayout(filaGrosor)

        self.botonPrevisualizar = qt.QPushButton(
            "Previsualizar los cortes (no corta nada)")
        self.botonPrevisualizar.toolTip = (
            "Dibuja, punto por punto, hasta donde llega el corte. Azul: "
            "atraviesa el hueso de lado a lado. Rojo: se queda adentro y ahi el "
            "hueso NO se va a separar.")
        caja.addWidget(self.botonPrevisualizar)

        self.botonOcultarPrevia = qt.QPushButton("Ocultar la previsualizacion")
        self.botonOcultarPrevia.setVisible(False)
        caja.addWidget(self.botonOcultarPrevia)

        self.botonCortar = qt.QPushButton("Hacer los cortes")
        self.botonCortar.setStyleSheet(ESTILO_BOTON_PRINCIPAL)
        self.botonCortar.setEnabled(False)
        caja.addWidget(self.botonCortar)

        self.estadoCorte = self._estado("")
        caja.addWidget(self.estadoCorte)

        # El widget oficial de Slicer para colocar puntos. Se usa oculto,
        # aprovechando solo su mecanismo interno: es mas confiable que conectar
        # a mano los clicks del mouse con el nodo activo.
        self.placeWidget = slicer.qSlicerMarkupsPlaceWidget()
        self.placeWidget.setMRMLScene(slicer.mrmlScene)
        self.placeWidget.buttonsVisible = False
        caja.addWidget(self.placeWidget)
        self.placeWidget.hide()

        self.botonTrazar.connect("clicked(bool)", self.onTrazarLinea)
        self.botonTerminarLinea.connect("clicked(bool)", self.onTerminarLinea)
        self.botonPrevisualizar.connect("clicked(bool)", self.onPrevisualizar)
        self.botonOcultarPrevia.connect("clicked(bool)", self.onOcultarPrevia)
        self.botonCortar.connect("clicked(bool)", self.onCortar)

    def onTrazarLinea(self):
        self._quitarObservadorCurva()
        curva = self.logic.nuevaCurvaDeCorte(len(self._curvas) + 1)
        self._curvaEnCurso = curva
        self._observadorCurva = curva.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self._onPuntoColocado)
        self.placeWidget.setCurrentNode(curva)
        self.placeWidget.setPlaceModePersistency(True)
        self.placeWidget.setPlaceModeEnabled(True)
        self.botonTerminarLinea.setEnabled(True)
        self.botonTrazar.setEnabled(False)
        self._poner(self.estadoTrazado,
                    "Hacé click sobre el craneo en la vista 3D para ir marcando "
                    "el recorrido del corte. Puntos colocados: 0.", NARANJA)

    def _onPuntoColocado(self, caller, event):
        if self._curvaEnCurso is None:
            return
        n = self._curvaEnCurso.GetNumberOfControlPoints()
        self._poner(self.estadoTrazado,
                    "Puntos colocados: %d. Cuando llegues al final del corte, "
                    "apretá 'Terminar esta linea'." % n, NARANJA)

    def _quitarObservadorCurva(self):
        if self._observadorCurva is not None and self._curvaEnCurso is not None:
            try:
                self._curvaEnCurso.RemoveObserver(self._observadorCurva)
            except Exception:
                pass
        self._observadorCurva = None

    def onTerminarLinea(self):
        self.placeWidget.setPlaceModeEnabled(False)
        self._quitarObservadorCurva()
        curva = self._curvaEnCurso
        self.botonTerminarLinea.setEnabled(False)
        self.botonTrazar.setEnabled(True)

        n = 0 if curva is None else curva.GetNumberOfControlPoints()
        if n < 2:
            if curva is not None:
                slicer.mrmlScene.RemoveNode(curva)
            self._curvaEnCurso = None
            self._poner(self.estadoTrazado,
                        "La linea necesita al menos 2 puntos (colocaste %d). "
                        "Probá de nuevo." % n, ROJO)
            return

        self._curvas.append(curva)
        self._curvaEnCurso = None
        self._construirPanelLineas()
        self.botonCortar.setEnabled(True)
        self._poner(self.estadoTrazado,
                    "Linea guardada con %d puntos. Ya tenés %d linea(s). Podés "
                    "trazar otra o pasar a cortar."
                    % (n, len(self._curvas)), AZUL)

    def _construirPanelLineas(self):
        while self.layoutLineas.count():
            item = self.layoutLineas.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        for curva in list(self._curvas):
            fila = qt.QWidget()
            cajaFila = qt.QHBoxLayout(fila)
            cajaFila.setContentsMargins(0, 2, 0, 2)
            etiqueta = qt.QLabel("%s   -   %d puntos"
                                 % (curva.GetName(),
                                    curva.GetNumberOfControlPoints()))
            etiqueta.setStyleSheet("font-size: 11px;")
            cajaFila.addWidget(etiqueta, 3)

            botonVer = qt.QPushButton("Ver")
            botonVer.setFixedWidth(55)
            cajaFila.addWidget(botonVer)

            botonBorrar = qt.QPushButton("Borrar")
            botonBorrar.setFixedWidth(60)
            botonBorrar.setStyleSheet(ESTILO_BOTON_QUITAR)
            cajaFila.addWidget(botonBorrar)

            self.layoutLineas.addWidget(fila)
            botonVer.connect("clicked(bool)",
                             lambda _c, c=curva: self.onVerLinea(c))
            botonBorrar.connect("clicked(bool)",
                                lambda _c, c=curva: self.onBorrarLinea(c))

    def onVerLinea(self, curva):
        for c in self._curvas:
            d = c.GetDisplayNode()
            if d is None:
                continue
            # se resalta solo con color: el grosor queda el de Slicer
            if c is curva:
                d.SetSelectedColor(1.0, 1.0, 0.0)
                d.SetColor(1.0, 1.0, 0.0)
            else:
                d.SetSelectedColor(1.0, 0.2, 0.2)
                d.SetColor(1.0, 0.4, 0.4)

    def onBorrarLinea(self, curva):
        nombre = curva.GetName()
        self._curvas = [c for c in self._curvas if c is not curva]
        try:
            slicer.mrmlScene.RemoveNode(curva)
        except Exception:
            pass
        self._construirPanelLineas()
        self.botonCortar.setEnabled(bool(self._curvas))
        self._poner(self.estadoTrazado,
                    "Se borro %s. Quedan %d linea(s)."
                    % (nombre, len(self._curvas)), AZUL)

    def onPrevisualizar(self):
        if not self._curvas:
            self._poner(self.estadoCorte,
                        "Primero trazá al menos una linea de corte.", ROJO)
            return
        self._esperando(self.estadoCorte,
                        "Calculando la previsualizacion, esperá...")
        try:
            resultado = self.logic.previsualizarCortes(self.spinGrosor.value)
        finally:
            self._listo()

        if resultado is None:
            self._poner(self.estadoCorte,
                        "No se pudo calcular la previsualizacion. Mirá el "
                        "detalle de abajo.", ROJO)
            return

        lineas = ["Previsualizacion lista. En la vista 3D:",
                  "   AZUL = el corte atraviesa el hueso de lado a lado.",
                  "   ROJO = el corte se queda adentro y ahi NO se va a separar.",
                  ""]
        problemas = False
        for c in resultado["curvas"]:
            if c.get("error"):
                lineas.append("- %s: %s" % (c["nombre"], c["error"]))
                problemas = True
                continue
            pct = 100.0 * (c["nPuntos"] - c["nNoAtraviesa"]) / max(c["nPuntos"], 1)
            lineas.append("- %s: atraviesa en el %.0f%% del recorrido, hueso de "
                          "%.1f mm de espesor"
                          % (c["nombre"], pct, c["espesorMediano"]))
            if c["nNoAtraviesa"]:
                problemas = True
                lineas.append("     %d punto(s) en rojo: en esos puntos ninguna "
                              "inclinacion de la sierra logra atravesar el hueso "
                              "(tipico detras de los ojos). Si ese tramo tiene que "
                              "separarse, trazá una linea extra sobre la cara por "
                              "la que si se puede entrar." % c["nNoAtraviesa"])
            if c.get("tiradaMM", 0.0) > 0:
                lineas.append("     Tramo continuo mas largo sin cortar: %.1f mm. "
                              "Este es el numero que importa: un solo tramo sin "
                              "cortar mantiene la pieza unida, aunque el resto "
                              "este cortado." % c["tiradaMM"])
            for h in c.get("huecosUnion", []):
                problemas = True
                lineas.append("     La punta del %s queda a %.1f mm de otra linea "
                              "y no llega a tocarla: queda un puente de hueso de "
                              "%.1f mm. Si esas dos lineas tienen que unirse, "
                              "acercá las puntas."
                              % (h["punta"], h["distancia"], h["hueco"]))
        self._poner(self.estadoCorte, "\n".join(lineas),
                    AMBAR if problemas else VERDE)
        self.botonOcultarPrevia.setVisible(True)

    def onOcultarPrevia(self):
        self.logic.ocultarPrevisualizacion()
        self.botonOcultarPrevia.setVisible(False)

    def onCortar(self):
        if not self._curvas:
            self._poner(self.estadoCorte,
                        "Primero trazá al menos una linea de corte.", ROJO)
            return
        self._esperando(self.estadoCorte,
                        "Haciendo los cortes. Esto puede tardar un par de "
                        "minutos...")
        try:
            resultado = self.logic.cortar(self.spinGrosor.value)
        finally:
            self._listo()

        if resultado is None:
            self._poner(self.estadoCorte,
                        "No se pudo hacer el corte. Mirá el detalle de abajo.",
                        ROJO)
            return

        lineas = ["Cortes realizados: quedaron %d pieza(s), %d listas para "
                  "imprimir." % (resultado["total"], resultado["watertight"])]
        for p in resultado["piezas"]:
            lineas.append("   %-22s %6.1f cm3%s"
                          % (p["nombre"], p["volumenCM3"],
                             "" if p["watertight"] else "   (malla abierta)"))
        if resultado["noSepararon"]:
            lineas.append("")
            lineas.append("ATENCION: %d linea(s) quedaron marcadas pero NO "
                          "separaron el hueso: %s. La pieza que tenia que "
                          "salir sigue pegada al resto del craneo, asi que en "
                          "el Paso 4 no se va a poder mover sola. Mirá los "
                          "puntos rojos de la previsualizacion, corregí esa "
                          "linea y volvé a cortar."
                          % (len(resultado["noSepararon"]),
                             ", ".join(resultado["noSepararon"])))
            lineas.append("Recordá: una linea abierta sola nunca separa una "
                          "pieza del craneo. Para separarla, las lineas tienen "
                          "que cerrar una vuelta completa, unidas en los dos "
                          "extremos.")
        puente = resultado.get("puente")
        if puente and puente["desconocidas"]:
            lineas.append("")
            lineas.append("AVISO: %d pieza(s) quedaron sin clasificar: %s. En "
                          "el Paso 4 se decide si se mueven segun lo cerca que "
                          "esten de una linea de corte."
                          % (len(puente["desconocidas"]),
                             ", ".join(puente["desconocidas"])))
        hayAvisos = bool(resultado["noSepararon"]) or bool(
            puente and puente["desconocidas"])
        self._poner(self.estadoCorte, "\n".join(lineas),
                    AMBAR if hayAvisos else VERDE)

        self._habilitar(self.paso4, True)
        self._habilitar(self.paso5, True)
        self.paso5.collapsed = True
        self.paso3.collapsed = True
        self._resetearVista3D()

    # ========================================================
    # PASO 4 - REACOMODAR LAS PIEZAS
    # ========================================================
    def _construirPaso4(self):
        self.paso4 = ctk.ctkCollapsibleButton()
        self.paso4.text = "Paso 4 - Reacomodar las piezas"
        self.layout.addWidget(self.paso4)
        caja = qt.QVBoxLayout(self.paso4)

        caja.addWidget(self._ayuda(
            "Ahora se rearma la boveda. La base del craneo queda fija (es la "
            "referencia de toda la cirugia) y el resto de las piezas se mueven "
            "respecto de ella. Podés arrastrarlas con el mouse o moverlas con "
            "valores exactos en milimetros."))

        self.botonPrepararPiezas = qt.QPushButton("Preparar las piezas")
        self.botonPrepararPiezas.setStyleSheet(ESTILO_BOTON_ACCION)
        caja.addWidget(self.botonPrepararPiezas)

        self.estadoPiezas = self._estado("Todavia no se prepararon las piezas.")
        caja.addWidget(self.estadoPiezas)

        self.contenedorFragmentos = qt.QWidget()
        self.layoutFragmentos = qt.QVBoxLayout(self.contenedorFragmentos)
        self.layoutFragmentos.setContentsMargins(0, 0, 0, 0)
        caja.addWidget(self.contenedorFragmentos)

        # --- controles de la pieza elegida ---
        self.panelMovimiento = qt.QWidget()
        cajaMov = qt.QVBoxLayout(self.panelMovimiento)
        cajaMov.setContentsMargins(0, 6, 0, 0)
        cajaMov.addWidget(self._separador())
        self.etiquetaPiezaActiva = self._titulo("Ninguna pieza seleccionada")
        cajaMov.addWidget(self.etiquetaPiezaActiva)
        cajaMov.addWidget(self._ayuda(
            "'Mover con el mouse' enciende las flechas y los aros sobre la "
            "pieza en la vista 3D. Los botones de abajo hacen lo mismo con "
            "valores exactos."))

        filaPaso = qt.QHBoxLayout()
        filaPaso.addWidget(qt.QLabel("Mover de a (mm):"))
        self.spinPasoMM = qt.QDoubleSpinBox()
        self.spinPasoMM.setRange(0.5, 30.0)
        self.spinPasoMM.setSingleStep(0.5)
        self.spinPasoMM.setValue(2.0)
        filaPaso.addWidget(self.spinPasoMM)
        filaPaso.addWidget(qt.QLabel("   Girar de a (grados):"))
        self.spinPasoGrados = qt.QDoubleSpinBox()
        self.spinPasoGrados.setRange(1.0, 90.0)
        self.spinPasoGrados.setSingleStep(1.0)
        self.spinPasoGrados.setValue(5.0)
        filaPaso.addWidget(self.spinPasoGrados)
        filaPaso.addStretch(1)
        cajaMov.addLayout(filaPaso)

        # Traslaciones, con los nombres anatomicos escritos como los dice el
        # cirujano y no como ejes RAS. dr = derecha, da = adelante, ds = arriba.
        filaMover = qt.QHBoxLayout()
        filaMover.addWidget(qt.QLabel("Mover:"))
        for texto, eje, signo in (("Izquierda", "dr", -1), ("Derecha", "dr", 1),
                                  ("Adelante", "da", 1), ("Atras", "da", -1),
                                  ("Arriba", "ds", 1), ("Abajo", "ds", -1)):
            b = qt.QPushButton(texto)
            b.setFixedWidth(78)
            filaMover.addWidget(b)
            b.connect("clicked(bool)",
                      lambda _c, e=eje, s=signo: self.onMoverPieza(e, s))
        cajaMov.addLayout(filaMover)

        filaGiro = qt.QHBoxLayout()
        filaGiro.addWidget(qt.QLabel("Girar:"))
        for texto, eje, signo in (("Cabeceo +", "LR", 1), ("Cabeceo -", "LR", -1),
                                  ("Rolido +", "AP", 1), ("Rolido -", "AP", -1),
                                  ("Guinada +", "SI", 1), ("Guinada -", "SI", -1)):
            b = qt.QPushButton(texto)
            b.setFixedWidth(78)
            filaGiro.addWidget(b)
            b.connect("clicked(bool)",
                      lambda _c, e=eje, s=signo: self.onGirarPieza(e, s))
        cajaMov.addLayout(filaGiro)

        filaEspecial = qt.QHBoxLayout()
        self.botonVoltear = qt.QPushButton("Dar vuelta 180 (adelante <-> atras)")
        self.botonVoltear.toolTip = (
            "Gira la pieza media vuelta, como cuando se pasa una placa de "
            "occipital a frontal. No es un espejo: el hueso no se puede "
            "reflejar, se voltea.")
        filaEspecial.addWidget(self.botonVoltear)
        self.botonAjustarAlMolde = qt.QPushButton("Acercar al molde")
        self.botonAjustarAlMolde.toolTip = (
            "Apoya la pieza sobre el casquete objetivo. Es una posicion de "
            "arranque: despues hay que ajustar a mano.")
        filaEspecial.addWidget(self.botonAjustarAlMolde)
        cajaMov.addLayout(filaEspecial)

        caja.addWidget(self.panelMovimiento)
        self.panelMovimiento.setVisible(False)

        caja.addWidget(self._separador())

        filaMolde = qt.QHBoxLayout()
        filaMolde.addWidget(qt.QLabel("Indice cefalico objetivo:"))
        self.spinIC = qt.QDoubleSpinBox()
        self.spinIC.setRange(60.0, 95.0)
        self.spinIC.setSingleStep(0.5)
        self.spinIC.setValue(BloqueG.IC_OBJETIVO)
        self.spinIC.toolTip = ("Un craneo normal esta entre 76 y 81. En "
                               "escafocefalia el indice es mucho mas bajo.")
        filaMolde.addWidget(self.spinIC)
        self.botonMolde = qt.QPushButton("Mostrar el craneo objetivo")
        self.botonMolde.toolTip = (
            "Dibuja un casquete semitransparente con la forma a la que hay que "
            "llegar, como el molde blanco que usan en el quirofano.")
        filaMolde.addWidget(self.botonMolde)
        caja.addLayout(filaMolde)

        filaMedidas = qt.QHBoxLayout()
        self.botonMedidas = qt.QPushButton("Ver las medidas")
        filaMedidas.addWidget(self.botonMedidas)
        self.botonLandmarks = qt.QPushButton("Marcar la linea media")
        self.botonLandmarks.toolTip = (
            "Crea tres puntos (nasion, bregma, inion) para comprobar que la "
            "cabeza no este torcida en la tomografia, que sesgaria las medidas.")
        filaMedidas.addWidget(self.botonLandmarks)
        self.botonHuecos = qt.QPushButton("Ver separacion entre piezas")
        self.botonHuecos.toolTip = (
            "Para cada pieza, cuantos milimetros la separan de la pieza mas "
            "cercana. Sirve para decidir donde hace falta injerto o placa. "
            "0 mm = se tocan o se superponen.")
        filaMedidas.addWidget(self.botonHuecos)
        self.botonResetear = qt.QPushButton("Devolver todo a su lugar")
        filaMedidas.addWidget(self.botonResetear)
        caja.addLayout(filaMedidas)

        self.estadoMedidas = self._estado("")
        caja.addWidget(self.estadoMedidas)

        self.botonPrepararPiezas.connect("clicked(bool)", self.onPrepararPiezas)
        self.botonVoltear.connect("clicked(bool)", self.onVoltearPieza)
        self.botonAjustarAlMolde.connect("clicked(bool)", self.onAjustarAlMolde)
        self.botonMolde.connect("clicked(bool)", self.onMostrarMolde)
        self.botonMedidas.connect("clicked(bool)", self.onVerMedidas)
        self.botonLandmarks.connect("clicked(bool)", self.onMarcarLineaMedia)
        self.botonResetear.connect("clicked(bool)", self.onResetearPiezas)
        self.botonHuecos.connect("clicked(bool)", self.onVerHuecos)

    def onPrepararPiezas(self):
        self._esperando(self.estadoPiezas, "Preparando las piezas...")
        try:
            filas = self.logic.prepararPiezas()
        finally:
            self._listo()
        if not filas:
            self._poner(self.estadoPiezas,
                        "No encontre las piezas del corte. Volvé al Paso 3.", ROJO)
            return
        self._construirPanelFragmentos(filas)
        movibles = sum(1 for f in filas if f["movible"])
        self._poner(self.estadoPiezas,
                    "%d pieza(s) en total: %d se pueden mover. La base '%s' "
                    "queda fija; si esa no es la base del craneo, marcá otra con "
                    "'Es la base'." % (len(filas), movibles,
                                       self.logic.bloqueG.ancla or "-"), AZUL)
        self._habilitar(self.paso5, True)

    def _construirPanelFragmentos(self, filas):
        while self.layoutFragmentos.count():
            item = self.layoutFragmentos.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._filasFragmentos = []

        for f in filas:
            fila = qt.QWidget()
            cajaFila = qt.QHBoxLayout(fila)
            cajaFila.setContentsMargins(0, 2, 0, 2)

            texto = "%s   %.1f cm3   -   %s" % (f["nombre"], f["volumen"],
                                                f["estado"])
            if f["desplazamiento"] > 0.05:
                texto += "  (movida %.1f mm)" % f["desplazamiento"]
            etiqueta = qt.QLabel(texto)
            etiqueta.setWordWrap(True)
            if f["esAncla"]:
                etiqueta.setStyleSheet("font-size: 11px; font-weight: bold;")
            elif not f["activa"]:
                etiqueta.setStyleSheet("font-size: 11px; color: #999999;")
            else:
                etiqueta.setStyleSheet("font-size: 11px;")
            cajaFila.addWidget(etiqueta, 3)

            botonMover = qt.QPushButton("Mover con el mouse")
            botonMover.setFixedWidth(130)
            botonMover.setEnabled(f["movible"])
            cajaFila.addWidget(botonMover)

            # Corrige el criterio automatico: una pieza que no toco ningun
            # corte queda fija, y una que paso cerca de una linea queda movible.
            if f["preexistente"]:
                botonFija = qt.QPushButton("Liberar")
                botonFija.toolTip = ("Esta pieza quedo fija porque ninguna "
                                     "linea de corte pasa cerca. Liberala "
                                     "si hay que moverla.")
            else:
                botonFija = qt.QPushButton("Dejar fija")
                botonFija.toolTip = ("Para piezas que no se tienen que mover "
                                     "(mandibula, vertebras, tubos) y quedaron "
                                     "movibles por estar cerca de una linea.")
            botonFija.setFixedWidth(80)
            botonFija.setEnabled(not f["esAncla"] and f["activa"])
            cajaFila.addWidget(botonFija)

            botonBase = qt.QPushButton("Es la base")
            botonBase.setFixedWidth(80)
            botonBase.setEnabled(not f["esAncla"] and f["activa"])
            cajaFila.addWidget(botonBase)

            if f["activa"]:
                botonQuitar = qt.QPushButton("Sacar")
                botonQuitar.setStyleSheet(ESTILO_BOTON_QUITAR)
            else:
                botonQuitar = qt.QPushButton("Devolver")
            botonQuitar.setFixedWidth(70)
            botonQuitar.setEnabled(not f["esAncla"])
            cajaFila.addWidget(botonQuitar)

            self.layoutFragmentos.addWidget(fila)
            self._filasFragmentos.append(f)

            nombre = f["nombre"]
            botonMover.connect("clicked(bool)",
                               lambda _c, n=nombre: self.onSeleccionarPieza(n))
            botonBase.connect("clicked(bool)",
                              lambda _c, n=nombre: self.onFijarBase(n))
            preexistente = f["preexistente"]
            botonFija.connect(
                "clicked(bool)",
                lambda _c, n=nombre, pre=preexistente: self.onFijarOLiberar(n, pre))
            activa = f["activa"]
            botonQuitar.connect(
                "clicked(bool)",
                lambda _c, n=nombre, a=activa: self.onSacarODevolver(n, a))

    def _refrescarFragmentos(self):
        self._construirPanelFragmentos(self.logic.bloqueG.tablaDePiezas())

    def onSeleccionarPieza(self, nombre):
        if not self.logic.bloqueG.manipular(nombre, True):
            return
        self._piezaSeleccionada = nombre
        self.etiquetaPiezaActiva.setText("Pieza seleccionada: %s" % nombre)
        self.panelMovimiento.setVisible(True)
        self._poner(self.estadoPiezas,
                    "Arrastrá las flechas para mover '%s' y los aros para "
                    "girarla, o usá los botones de abajo." % nombre, AZUL)

    def onFijarBase(self, nombre):
        self.logic.bloqueG.fijar(nombre)
        self._refrescarFragmentos()
        self._poner(self.estadoPiezas,
                    "'%s' quedo como base fija del armado." % nombre, AZUL)

    def onFijarOLiberar(self, nombre, estabaFija):
        if estabaFija:
            self.logic.bloqueG.liberar_pieza(nombre)
            texto = "'%s' ahora se puede mover." % nombre
        else:
            self.logic.bloqueG.fijar_pieza(nombre)
            texto = ("'%s' queda fija y no se cuenta en las medidas."
                     % nombre)
            if self._piezaSeleccionada == nombre:
                self._piezaSeleccionada = None
                self.panelMovimiento.setVisible(False)
        self._refrescarFragmentos()
        self._poner(self.estadoPiezas, texto, AZUL)

    def onSacarODevolver(self, nombre, estabaActiva):
        if estabaActiva:
            self.logic.bloqueG.descartar(nombre)
        else:
            self.logic.bloqueG.restaurar(nombre)
        if self._piezaSeleccionada == nombre and estabaActiva:
            self._piezaSeleccionada = None
            self.panelMovimiento.setVisible(False)
        self._refrescarFragmentos()

    def _piezaActiva(self):
        if not self._piezaSeleccionada:
            self._poner(self.estadoPiezas,
                        "Primero elegí una pieza con 'Mover con el mouse'.", ROJO)
            return None
        return self._piezaSeleccionada

    def onMoverPieza(self, eje, signo):
        nombre = self._piezaActiva()
        if not nombre:
            return
        paso = self.spinPasoMM.value * signo
        kwargs = {"dr": 0.0, "da": 0.0, "ds": 0.0}
        kwargs[eje] = paso
        self.logic.bloqueG.mover(nombre, **kwargs)
        self._refrescarFragmentos()

    def onGirarPieza(self, eje, signo):
        nombre = self._piezaActiva()
        if not nombre:
            return
        self.logic.bloqueG.rotar(nombre, eje, self.spinPasoGrados.value * signo)
        self._refrescarFragmentos()

    def onVoltearPieza(self):
        nombre = self._piezaActiva()
        if not nombre:
            return
        self.logic.bloqueG.voltear(nombre, "SI")
        self._refrescarFragmentos()

    def onAjustarAlMolde(self):
        nombre = self._piezaActiva()
        if not nombre:
            return
        if not self.logic.bloqueG.ajustar(nombre):
            self._poner(self.estadoMedidas,
                        "Primero generá el craneo objetivo.", ROJO)
            return
        self._refrescarFragmentos()

    def onMostrarMolde(self):
        info = self.logic.bloqueG.molde(self.spinIC.value)
        if info is None:
            self._poner(self.estadoMedidas,
                        "No hay piezas activas para calcular el objetivo.", ROJO)
            return
        self._poner(self.estadoMedidas,
                    "Craneo actual: %.0f mm de largo por %.0f mm de ancho "
                    "(indice cefalico %.1f).\n"
                    "Objetivo: %.0f x %.0f mm (indice %.1f).\n"
                    "Moviendo las piezas se puede acortar de adelante hacia "
                    "atras; el ancho no se gana moviendo placas enteras."
                    % (info["largoActual"], info["anchoActual"], info["icActual"],
                       info["largoObjetivo"], info["anchoObjetivo"],
                       info["icObjetivo"]),
                    AZUL)

    def onVerMedidas(self):
        m = self.logic.bloqueG.metricas()
        if m is None:
            self._poner(self.estadoMedidas, "No hay piezas preparadas.", ROJO)
            return
        lineas = [
            "Largo (adelante-atras): %.0f mm" % m["largo"],
            "Ancho (lado a lado):    %.0f mm" % m["ancho"],
            "Alto:                   %.0f mm" % m["alto"],
            "Indice cefalico: %.1f   (objetivo %.1f)" % (m["ic"], m["icObjetivo"]),
        ]
        if m["ic"] < m["icObjetivo"] - 1.0:
            lineas.append("Faltan %.1f puntos: seguí acercando las piezas de "
                          "adelante y de atras." % (m["icObjetivo"] - m["ic"]))
        elif m["ic"] > m["icObjetivo"] + 1.0:
            lineas.append("Te pasaste %.1f puntos del objetivo."
                          % (m["ic"] - m["icObjetivo"]))
        else:
            lineas.append("El indice cefalico esta dentro del objetivo.")
        if m["descartadas"]:
            lineas.append("Hueso resecado: %.1f cm3 en %d pieza(s)."
                          % (m["volumenDescartado"], len(m["descartadas"])))
        if m["nivelS"] is None:
            lineas.append("ATENCION: el indice se midio sobre todo el craneo, "
                          "incluida la base y la cara, asi que puede no ser "
                          "exacto (suele salir mas bajo que el real).")
        if m["excluidas"]:
            lineas.append("No se cuentan en las medidas (no las toco ningun "
                          "corte): %s." % ", ".join(m["excluidas"]))
        if m["avisoInclinacion"]:
            lineas.append("ATENCION: " + m["avisoInclinacion"])
        self._poner(self.estadoMedidas, "\n".join(lineas),
                    AMBAR if (m["avisoInclinacion"] or m["nivelS"] is None)
                    else AZUL)

    def onMarcarLineaMedia(self):
        self.logic.bloqueG.landmarks()
        self._poner(self.estadoMedidas,
                    "Se crearon los tres puntos de linea media. Colocalos sobre "
                    "el craneo desde el modulo Markups (nasion en la raiz de la "
                    "nariz, bregma arriba, inion atras) y despues volvé a "
                    "'Ver las medidas'.", AZUL)

    def onVerHuecos(self):
        self._esperando(self.estadoMedidas,
                        "Midiendo la separacion entre piezas...")
        inicio = time.time()
        try:
            filas = self.logic.bloqueG.huecos()
        finally:
            self._listo()
        segundos = time.time() - inicio
        self.logic.log("Separacion entre piezas calculada en %.1f s." % segundos)
        if not filas:
            self._poner(self.estadoMedidas,
                        "Hacen falta al menos dos piezas en el armado.", ROJO)
            return
        lineas = ["Separacion de cada pieza con la mas cercana:"]
        for h in filas:
            nombre = h["pieza"] + (" (%s)" % h["etiqueta"] if h["etiqueta"] else "")
            lineas.append("   %-24s -> %-22s %5.1f mm"
                          % (nombre, h["vecina"] or "-", h["mm"]))
        lineas.append("0 mm = las piezas se tocan o se superponen (el modulo "
                      "no detecta superposiciones). Es una medida aproximada.")
        lineas.append("(calculado en %.1f s)" % segundos)
        self._poner(self.estadoMedidas, "\n".join(lineas), AZUL)

    def onResetearPiezas(self):
        n = self.logic.bloqueG.resetear()
        self._refrescarFragmentos()
        self._poner(self.estadoPiezas,
                    "%d pieza(s) volvieron a su posicion original." % n, AZUL)

    # ========================================================
    # PASO 5 - EXPORTAR
    # ========================================================
    def _construirPaso5(self):
        self.paso5 = ctk.ctkCollapsibleButton()
        self.paso5.text = "Paso 5 - Exportar el resultado"
        self.layout.addWidget(self.paso5)
        caja = qt.QVBoxLayout(self.paso5)

        caja.addWidget(self._ayuda(
            "Guarda cada pieza como archivo STL, ya con los movimientos "
            "aplicados. Son los archivos que van a la impresora 3D para hacer "
            "el molde y las guias de corte."))

        self.botonExportarArmado = qt.QPushButton(
            "Guardar el armado final (STL)...")
        self.botonExportarArmado.setStyleSheet(ESTILO_BOTON_PRINCIPAL)
        caja.addWidget(self.botonExportarArmado)

        self.botonExportarCrudo = qt.QPushButton(
            "Guardar las piezas sin mover (STL)...")
        self.botonExportarCrudo.toolTip = (
            "Las piezas tal como quedaron despues del corte, en su posicion "
            "original. Sirve para imprimir el craneo del paciente.")
        caja.addWidget(self.botonExportarCrudo)

        self.estadoExportar = self._estado("")
        caja.addWidget(self.estadoExportar)

        caja.addWidget(self._ayuda(
            "Pendiente: reporte prequirurgico en PDF con las capturas, las "
            "medidas y los movimientos de cada pieza."))

        self.botonExportarArmado.connect("clicked(bool)", self.onExportarArmado)
        self.botonExportarCrudo.connect("clicked(bool)", self.onExportarCrudo)

    def _pedirCarpeta(self):
        return qt.QFileDialog.getExistingDirectory(
            self.parent, "Elegi la carpeta donde guardar los STL")

    def onExportarArmado(self):
        carpeta = self._pedirCarpeta()
        if not carpeta:
            return
        self._esperando(self.estadoExportar, "Guardando...")
        try:
            n = self.logic.bloqueG.exportarSTL(carpeta)
        finally:
            self._listo()
        self._poner(self.estadoExportar,
                    "Se guardaron %d archivo(s) STL en:\n%s" % (n, carpeta),
                    VERDE if n else ROJO)

    def onExportarCrudo(self):
        carpeta = self._pedirCarpeta()
        if not carpeta:
            return
        self._esperando(self.estadoExportar, "Guardando...")
        try:
            n = BloqueF.exportar_stl(carpeta)
        finally:
            self._listo()
        self._poner(self.estadoExportar,
                    "Se guardaron %d archivo(s) STL en:\n%s" % (n, carpeta),
                    VERDE if n else ROJO)

    # ========================================================
    # CONSOLA DE DETALLE
    # ========================================================
    def _construirConsola(self):
        self.consola = ctk.ctkCollapsibleButton()
        self.consola.text = "Detalle tecnico"
        self.consola.collapsed = True
        self.layout.addWidget(self.consola)
        caja = qt.QVBoxLayout(self.consola)
        caja.addWidget(self._ayuda(
            "Todo lo que fue haciendo el sistema, con los numeros. Sirve para "
            "revisar un caso raro o para pegarlo en un reporte de error."))
        self.textoConsola = qt.QTextEdit()
        self.textoConsola.setReadOnly(True)
        self.textoConsola.setMinimumHeight(220)
        self.textoConsola.setStyleSheet("font-family: monospace; font-size: 11px;")
        caja.addWidget(self.textoConsola)
        botonLimpiar = qt.QPushButton("Limpiar")
        caja.addWidget(botonLimpiar)
        botonLimpiar.connect("clicked(bool)",
                             lambda _c: self.textoConsola.clear())

    def _log(self, mensaje):
        """Todo lo que imprimen los bloques pasa por aca: va a la consola de
        Python (para nosotros) y al panel de detalle (para el medico, que no
        abre la consola de Python)."""
        texto = str(mensaje)
        print(texto)
        try:
            self.textoConsola.append(texto)
            barra = self.textoConsola.verticalScrollBar()
            barra.setValue(barra.maximum)
        except Exception:
            pass

    # ========================================================
    # CICLO DE VIDA
    # ========================================================
    def _resetearVista3D(self):
        try:
            slicer.app.layoutManager().threeDWidget(0).threeDView().resetFocalPoint()
        except Exception:
            pass

    def cleanup(self) -> None:
        self.removeObservers()
        if self.logic is not None:
            self.logic.limpiarTemporal()

    def enter(self) -> None:
        self.initializeParameterNode()

    def exit(self) -> None:
        if self._parameterNode and self._parameterNodeGuiTag:
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
        if self._parameterNode and self._parameterNodeGuiTag:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
            self._parameterNodeGuiTag = None
        self._parameterNode = inputParameterNode


# ============================================================
# LOGICA - orquesta los bloques, sin un solo widget adentro
# ============================================================
class CranioPlanLogic(ScriptedLoadableModuleLogic):
    """
    No implementa algoritmos: los llama en orden y traduce sus resultados.

    Esa es la razon de que sea corta. Si algo de aca empieza a crecer, casi
    seguro pertenece a un bloque de CranioPlanLib.
    """

    def __init__(self, log=None) -> None:
        ScriptedLoadableModuleLogic.__init__(self)
        self.log = log if log is not None else print
        self.bloqueA = BloqueA.BloqueA(log=self.log)
        self.bloqueG = BloqueG.BloqueG(log=self.log)
        # Carpeta del estudio en curso. Si el medico eligio un .zip, es la
        # carpeta temporal donde se descomprimio, y hay que borrarla al final.
        self._carpetaEstudio = None
        self._carpetaTemporal = None
        # El Bloque F trabaja con constantes de modulo (asi sigue siendo
        # pegable en la consola tal cual). Se le enchufa el mismo log y el
        # prefijo de curvas para que no agarre curvas viejas de la escena.
        BloqueF.configurar(log=self.log, PREFIJO_CURVAS=Comun.PREFIJO_CURVA_CORTE)
        PuenteFG.configurar(log=self.log)

    def getParameterNode(self):
        return CranioPlanParameterNode(super().getParameterNode())

    # -------- Paso 1 --------
    def analizarEstudio(self, ruta):
        """
        Acepta una carpeta o un .zip. Devuelve el ranking de series.

        La carpeta del estudio queda guardada aca y no en el widget porque es
        la Logic la que despues tiene que cargar la serie: si el widget la
        perdiera (por ejemplo al cerrarse la escena), la carga fallaria con un
        error que no dice nada.
        """
        self.limpiarTemporal()
        try:
            carpeta, temporal = BloqueC.prepararCarpeta(ruta)
        except Exception as e:
            self.log("CranioPlan: no pude abrir el .zip: %s" % e)
            return None
        self._carpetaEstudio = carpeta
        self._carpetaTemporal = temporal
        return BloqueC.analizarEstudio(carpeta, log=self.log)

    def cargarSerie(self, seriesUID):
        if not self._carpetaEstudio:
            self.log("CranioPlan: no hay ningun estudio analizado todavia.")
            return None
        return BloqueC.cargarSerie(self._carpetaEstudio, seriesUID, log=self.log)

    def limpiarTemporal(self):
        """Borra la carpeta temporal del .zip, si habia."""
        if self._carpetaTemporal:
            BloqueC.borrarTemporal(self._carpetaTemporal)
            self._carpetaTemporal = None

    # -------- Paso 2 --------
    def generarCraneo(self, volumeNode, modoPostop=False):
        return self.bloqueA.generar(volumeNode, modoPostop=modoPostop)

    def cambiarPiezaPrincipal(self, numero):
        return self.bloqueA.redecidir(referenciaForzada=int(numero))

    def confirmarCraneo(self):
        resultado = self.bloqueA.confirmarCraneo()
        if resultado is None:
            return None
        modelo = self.bloqueA.enviarAPlanner()
        if modelo is None:
            return None
        ok, lineas = self.bloqueA.verificarCompatibilidad()
        for l in lineas:
            self.log("CranioPlan: " + l)
        if not ok:
            return None
        return resultado

    # -------- Paso 3 --------
    def nuevaCurvaDeCorte(self, indice):
        """
        Crea una curva ABIERTA para una osteotomia.

        Solo abiertas: el lazo cerrado del Bloque F sigue soportado en el
        codigo, pero no se ofrece en la interfaz porque en la practica todas
        las osteotomias se trazan como lineas y tener las dos opciones era una
        forma mas de equivocarse.
        """
        nombre = slicer.mrmlScene.GenerateUniqueName(
            "%s%d" % (Comun.PREFIJO_CURVA_CORTE, indice))
        curva = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsCurveNode", nombre)
        curva.CreateDefaultDisplayNodes()
        d = curva.GetDisplayNode()
        if d is not None:
            # Tamano de puntos y grosor de linea: los de Slicer por defecto
            # (no se tocan), igual que al trazar desde la consola. Antes se
            # forzaban mas grandes y la linea se veia como un tubo grueso.
            d.SetSelectedColor(1.0, 0.2, 0.2)
            d.SetColor(1.0, 0.4, 0.4)
            d.SetPropertiesLabelVisibility(False)
            d.SetPointLabelsVisibility(False)
        return curva

    def previsualizarCortes(self, grosorMM):
        BloqueF.configurar(GROSOR_CORTE_MM=float(grosorMM))
        return BloqueF.diagnosticar()

    def ocultarPrevisualizacion(self):
        BloqueF.ocultarPicos(False)

    def cortar(self, grosorMM):
        BloqueF.configurar(GROSOR_CORTE_MM=float(grosorMM))
        resultado = BloqueF.cortar()
        if resultado is not None:
            # Puente F -> G: deja los nombres Tapa_/Resto_/Hueso_ y los
            # atributos CranioPlan.* en cada pieza, sin que el medico apriete
            # nada. Con el Bloque F arreglado solo escribe los atributos.
            resultado["puente"] = PuenteFG.puente()
        return resultado

    # -------- Paso 4 --------
    def prepararPiezas(self):
        return self.bloqueG.preparar()


# ============================================================
# TEST
# ============================================================
class CranioPlanTest(ScriptedLoadableModuleTest):

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.delayDisplay("Sin tests automatizados por ahora.")
