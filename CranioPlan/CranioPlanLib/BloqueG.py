# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - BLOQUE G v1.1  (version libreria, para el modulo)
#              REACOMODAMIENTO DE PIEZAS POST-CORTE
# ============================================================
# Entra DESPUES del Bloque F. Toma los modelos de la carpeta
# 'CranioPlan_Fragmentos' y permite descartar piezas, moverlas y rotarlas para
# rearmar la boveda craneal (remodelacion tipo Melbourne).
#
# Es el script desarrollo/valentino/Bloque_G_V1.1.py convertido en libreria.
# La logica es la de V1.1. Los cambios son de envoltorio:
#
#   - el estado global del script (_ESTADO, _ANCLA, _NIVEL_BOVEDA) vive
#     adentro de la clase BloqueG (self.estado, self.ancla, self.nivelBoveda),
#     para que el modulo pueda tirar el estado de un paciente y empezar otro
#     sin reiniciar Slicer;
#   - print() -> self.log(), que escribe en la consola Y en la interfaz;
#   - los metodos DEVUELVEN diccionarios/listas para que la interfaz muestre
#     un resumen (tablaDePiezas(), metricas(), molde(), huecos());
#   - IC_OBJETIVO del script -> self.icObjetivo, que elige el medico en la
#     interfaz (arranca en IC_OBJETIVO);
#   - nombres de nodos y carpetas desde Comun.py.
#
# ############################################################
# CAMBIOS v1.1 RESPECTO DE v1.0
# ############################################################
# [F1] LA MOVILIDAD YA NO SE DECIDE POR EL NOMBRE. Una pieza es movible si
#      su superficie pasa a menos de UMBRAL_CONTACTO_MM de alguna curva de
#      corte de la escena ("lo que toco el corte se mueve"). Si no hay curvas
#      en la escena, cae al criterio por nombre (Hueso_ = fija) y lo avisa.
#      v1.0 decidia solo por el prefijo Hueso_; con el bug de nombres del
#      Bloque F (piezas 'Segment_N') todo quedaba movible, incluida la
#      mandibula. v1.0 del modulo no tenia ningun parche para eso.
# [F2] EL INDICE CEFALICO SE MIDE SOLO SOBRE LA BOVEDA: por encima de un
#      nivel S tomado de la curva CERRADA de la escena, y solo sobre piezas
#      craneales (se excluyen las preexistentes: mandibula, vertebras,
#      tubuladuras).
# [F3] NUEVO: huecos(). Distancia minima de cada pieza a su vecina mas
#      cercana, en mm.
#
# ------------------------------------------------------------
# DECISIONES DE ARQUITECTURA (vigentes desde v1.0)
# ------------------------------------------------------------
# [G1] LA GEOMETRIA NUNCA SE TOCA. Cada pieza movible recibe su propio
#      vtkMRMLLinearTransformNode. Todo el movimiento vive en esa matriz:
#      reset instantaneo, la matriz ES el dato del reporte prequirurgico, y
#      se recalculan metricas sin recalcular geometria. El hardening ocurre
#      solo al exportar STL, sobre copias.
# [G3] HAY UNA PIEZA ANCLA (la base). Tocar el corte no alcanza para
#      moverse: la base craneal y el macizo facial tambien tocan el corte y
#      son el marco de referencia de toda la cirugia. Se detecta como la
#      pieza grande con el centroide mas inferior en S; se corrige con
#      fijar().
#         FIJA (gris claro):  piezas lejos de toda curva + el ANCLA
#         MOVIBLE (color):    piezas en contacto con alguna curva, menos el ancla
# [G4] EL CENTRO DE TRANSFORMACION VA EN EL CENTROIDE. Sin esto el gizmo rota
#      alrededor del origen RAS y la pieza sale volando en el primer arrastre.
# [G5] EL MOLDE ES OPCIONAL Y DE BAJO VALOR: un elipsoide liso derivado del
#      bbox del paciente, sin anatomia real. Referencia visual.
# [G6] NO EXISTE 'ESPEJAR' UNA PIEZA: el hueso no se refleja, se VOLTEA 180
#      grados (rotacion propia). Por eso el comando es voltear().
# [G7] DESCARTAR ES UN ESTADO, NO UN BORRADO: la pieza queda gris translucida
#      para el reporte (volumen resecado) y porque el cirujano cambia de idea.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (del script V1.1)
# ------------------------------------------------------------
# - El IC se calcula sobre los ejes RAS del estudio. Si la cabeza esta
#   rotada en el CT, ancho y largo salen sesgados. metricas() avisa si los
#   landmarks de linea media muestran desviacion.
# - Con transformaciones rigidas NO se puede ensanchar una placa. Ganar ancho
#   biparietal requiere barrel staves o partir el parietal en tiras.
# - ajustar() (ICP) contra un elipsoide liso desliza tangencialmente: pose
#   inicial, no resultado final.
# - huecos() mide distancia superficie-a-superficie por muestreo de puntos:
#   es una cota, no una medida exacta.
# - No hay deteccion de colisiones: dos piezas pueden interpenetrarse y el
#   modulo no avisa. huecos() devuelve 0.0 en ese caso.
#
# ------------------------------------------------------------
# CASOS CONOCIDOS EN QUE FALLA EN EL MODULO
# (identificados por revision de codigo, no observados todavia en casos reales)
# ------------------------------------------------------------
# - SIN NIVEL DE BOVEDA: el nivel sale de una curva CERRADA y el modulo solo
#   usa lineas abiertas. Entonces nivelBoveda queda en None y el IC se mide
#   sobre todas las piezas craneales, incluida la base y la cara: puede salir
#   mas bajo que el real. metricas() lo avisa. Solucion prevista en V1.2
#   (puntos anatomicos), primero en el script de consola.
# - Una pieza craneal a mas de UMBRAL_CONTACTO_MM de toda linea queda fija
#   aunque haya que moverla; una mandibula, vertebra o tubuladura a menos de
#   UMBRAL_CONTACTO_MM de una linea queda movible. Se corrigen a mano con
#   liberar_pieza() / fijar_pieza().
# - Se toman TODAS las curvas de la escena, no solo las 'Corte_': una curva
#   ajena cerca de una pieza la vuelve movible.
# - Si se corta dos veces en la misma escena, los modelos del primer corte
#   siguen en la carpeta y aparecen como piezas duplicadas.
# - huecos() y _distancia_a_curvas() recorren puntos con bucles Python: con
#   mallas grandes huecos() puede tardar. Se vectoriza primero en el script.
# ============================================================

import os

import numpy as np
import vtk
import slicer

from .Comun import NOMBRE_CARPETA_SH, NOMBRE_MOLDE, NOMBRE_LANDMARKS

BLOQUE_G_VERSION = "v1.1"
CRANIOPLAN_G_VERSION = "Bloque G v1.1"

# ==================== CONFIGURACION ====================
UMBRAL_CONTACTO_MM     = 5.0     # (F1) distancia pieza-curva para considerar
#   que la pieza nacio del corte. El kerf es ~1.2mm y la curva corre sobre la
#   tabla externa, asi que las piezas cortadas quedan a 1-3mm. La mandibula y
#   las vertebras estan a decenas de mm. 5.0 separa con holgura.
IC_OBJETIVO            = 78.0    # indice cefalico objetivo (normocefalia ~76-81)
PREFIJOS_FIJOS         = ("Hueso_",)   # solo se usa si NO hay curvas (fallback)
VOL_MIN_ANCLA_CM3      = 5.0     # el ancla tiene que ser una pieza grande
MUESTREO_HUECOS        = 4       # 1 de cada N puntos al medir huecos

COLOR_FIJA             = (0.85, 0.85, 0.85)
COLOR_DESCARTADA       = (0.55, 0.55, 0.55)
OPACIDAD_DESCARTADA    = 0.25
OPACIDAD_MOLDE         = 0.20

PALETA = [
    (0.90, 0.30, 0.30), (0.30, 0.65, 0.95), (0.35, 0.80, 0.45),
    (0.95, 0.75, 0.25), (0.75, 0.45, 0.90), (0.30, 0.85, 0.85),
    (0.95, 0.55, 0.30), (0.60, 0.80, 0.30),
]

EJES = {"LR": (1.0, 0.0, 0.0),
        "AP": (0.0, 1.0, 0.0),
        "SI": (0.0, 0.0, 1.0)}


# ============================================================
# UTILIDADES (no dependen del estado; identicas a V1.1)
# ============================================================
def _sh():
    return slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)


def _carpeta_fragmentos():
    shNode = _sh()
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(shNode.GetSceneItemID(), hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if shNode.GetItemName(itemId) == NOMBRE_CARPETA_SH:
            return shNode, itemId
    return shNode, None


def _modelos_de_carpeta(shNode, itemId):
    encontrados = []

    def recorrer(parent):
        hijos = vtk.vtkIdList()
        shNode.GetItemChildren(parent, hijos)
        for i in range(hijos.GetNumberOfIds()):
            hid = hijos.GetId(i)
            nodo = shNode.GetItemDataNode(hid)
            if nodo is not None and nodo.IsA("vtkMRMLModelNode"):
                encontrados.append(nodo)
            elif nodo is None:
                recorrer(hid)

    recorrer(itemId)
    return encontrados


def _centroide_y_volumen(modelNode):
    poly = modelNode.GetPolyData()
    if poly is None or poly.GetNumberOfPoints() == 0:
        return np.zeros(3), 0.0
    cdm = vtk.vtkCenterOfMass()
    cdm.SetInputData(poly)
    cdm.SetUseScalarsAsWeights(False)
    cdm.Update()
    centro = np.array(cdm.GetCenter(), dtype=float)
    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(poly)
    tri.Update()
    masa = vtk.vtkMassProperties()
    masa.SetInputConnection(tri.GetOutputPort())
    masa.Update()
    return centro, masa.GetVolume() / 1000.0


def _curvas_de_corte():
    """Todas las curvas de corte presentes en la escena."""
    cerradas = list(slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode'))
    abiertas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')]
    return cerradas + abiertas


def _puntos_de_curvas(curvas):
    """Array (N,3) con los puntos densos de todas las curvas juntas."""
    todos = []
    for c in curvas:
        pts = c.GetCurvePointsWorld()
        if pts is None:
            continue
        for i in range(pts.GetNumberOfPoints()):
            todos.append(pts.GetPoint(i))
    return np.asarray(todos, dtype=float) if todos else None


def _distancia_a_curvas(modelNode, puntosCurva):
    """
    (F1) Distancia minima entre la superficie de la pieza y el conjunto de
    puntos de las curvas de corte. Se construye un locator sobre la MALLA y
    se consulta cada punto de curva: las curvas tienen cientos de puntos,
    las mallas decenas de miles, asi que este sentido es el barato.
    """
    poly = modelNode.GetPolyData()
    if poly is None or poly.GetNumberOfPoints() == 0 or puntosCurva is None:
        return float('inf')
    loc = vtk.vtkPointLocator()
    loc.SetDataSet(poly)
    loc.BuildLocator()
    dmin = float('inf')
    for p in puntosCurva:
        idc = loc.FindClosestPoint(p.tolist())
        if idc < 0:
            continue
        q = np.array(poly.GetPoint(idc), dtype=float)
        d = float(np.linalg.norm(q - p))
        if d < dmin:
            dmin = d
            if dmin < 0.05:
                break
    return dmin


def _poly_mundo(modelNode):
    poly = modelNode.GetPolyData()
    if poly is None:
        return None
    tg = vtk.vtkGeneralTransform()
    slicer.vtkMRMLTransformNode.GetTransformBetweenNodes(
        modelNode.GetParentTransformNode(), None, tg)
    f = vtk.vtkTransformPolyDataFilter()
    f.SetInputData(poly)
    f.SetTransform(tg)
    f.Update()
    salida = vtk.vtkPolyData()
    salida.DeepCopy(f.GetOutput())
    return salida


# ============================================================
# LA CLASE QUE USA EL MODULO
# ============================================================
class BloqueG:
    """Estado del rearmado de UN paciente."""

    def __init__(self, log=print):
        self.log = log
        self.estado = {}          # nombre -> dict de la pieza   (_ESTADO)
        self.ancla = None         # nombre de la pieza ancla     (_ANCLA)
        self.nivelBoveda = None   # (F2) nivel S del corte basal (_NIVEL_BOVEDA)
        self.icObjetivo = IC_OBJETIVO

    # --------------------------------------------------------
    def _centroide_actual(self, info):
        m = vtk.vtkMatrix4x4()
        info['tNode'].GetMatrixTransformToWorld(m)
        c = info['centroide0']
        p = m.MultiplyPoint([c[0], c[1], c[2], 1.0])
        return np.array(p[:3], dtype=float)

    def _componer(self, tNode, matrizMundo):
        actual = vtk.vtkMatrix4x4()
        tNode.GetMatrixTransformToParent(actual)
        t = vtk.vtkTransform()
        t.PostMultiply()
        t.SetMatrix(actual)
        t.Concatenate(matrizMundo)
        tNode.SetMatrixTransformToParent(t.GetMatrix())

    def _buscar(self, nombre):
        if nombre in self.estado:
            return self.estado[nombre]
        for n, info in self.estado.items():
            if info['etiqueta'] and info['etiqueta'].upper() == str(nombre).upper():
                return info
        self.log(f"  no encuentro la pieza '{nombre}'. Corre piezas() para ver la lista.")
        return None

    def _nombre_de(self, info):
        for n, i in self.estado.items():
            if i is info:
                return n
        return "?"

    def _pintar(self, info):
        dn = info['node'].GetDisplayNode()
        if dn is None:
            info['node'].CreateDefaultDisplayNodes()
            dn = info['node'].GetDisplayNode()
        if dn is None:
            return
        if not info['activa']:
            dn.SetColor(*COLOR_DESCARTADA)
            dn.SetOpacity(OPACIDAD_DESCARTADA)
        elif not info['movible']:
            dn.SetColor(*COLOR_FIJA)
            dn.SetOpacity(1.0)
        else:
            dn.SetColor(*info['color'])
            dn.SetOpacity(1.0)
        dn.SetVisibility(True)

    # ========================================================
    # PREPARACION
    # ========================================================
    def preparar(self, umbral_contacto_mm=UMBRAL_CONTACTO_MM):
        """
        Inventaria los fragmentos, decide movilidad POR CONTACTO CON LAS CURVAS
        (F1), crea un transform por pieza movible (G1) y detecta el ancla (G3).
        Se puede volver a correr: resetea el estado completo.

        IMPORTANTE: las curvas de corte tienen que seguir en la escena. Si se
        borraron, esto cae al criterio por nombre y lo avisa.

        Devuelve la lista de piezas (tablaDePiezas) o None.
        """
        shNode, itemId = _carpeta_fragmentos()
        if itemId is None:
            self.log(f"ERROR: no encuentro la carpeta '{NOMBRE_CARPETA_SH}'. Corre cortar().")
            return None

        modelos = _modelos_de_carpeta(shNode, itemId)
        if not modelos:
            self.log("ERROR: la carpeta existe pero no tiene modelos adentro.")
            return None

        for info in self.estado.values():
            t = info.get('tNode')
            if t is not None and slicer.mrmlScene.IsNodePresent(t):
                info['node'].SetAndObserveTransformNodeID(None)
                slicer.mrmlScene.RemoveNode(t)
        self.estado = {}
        self.ancla = None
        self.nivelBoveda = None

        self.log(f"--- {CRANIOPLAN_G_VERSION} : preparar() ---")

        curvas = _curvas_de_corte()
        puntosCurva = _puntos_de_curvas(curvas)
        usaGeometria = puntosCurva is not None and len(puntosCurva) > 0

        if usaGeometria:
            self.log(f"Curvas de corte en escena: {[c.GetName() for c in curvas]}")
            self.log(f"Movilidad por contacto geometrico (umbral {umbral_contacto_mm:.1f} mm).")
            # (F2) nivel de boveda = S medio de la curva CERRADA (craneotomia basal)
            cerradas = [c for c in curvas if c.IsA('vtkMRMLMarkupsClosedCurveNode')]
            if cerradas:
                p = _puntos_de_curvas([cerradas[0]])
                if p is not None and len(p):
                    self.nivelBoveda = float(np.mean(p[:, 2]))
                    self.log(f"Nivel de boveda tomado de '{cerradas[0].GetName()}': "
                             f"S = {self.nivelBoveda:.1f} mm")
        else:
            self.log("AVISO: no hay curvas de corte en la escena.")
            self.log("  Caigo al criterio por NOMBRE, que falla si las piezas se llaman")
            self.log("  'Segment_N'. Volve a correr el corte sin borrar las curvas, o")
            self.log("  marca las fijas a mano con fijar_pieza('Nombre').")

        idxColor = 0
        for node in modelos:
            nombre = node.GetName()
            centro, vol = _centroide_y_volumen(node)

            if usaGeometria:
                d = _distancia_a_curvas(node, puntosCurva)
                movible = (d <= umbral_contacto_mm)
            else:
                d = float('nan')
                movible = not any(nombre.startswith(p) for p in PREFIJOS_FIJOS)

            tNode = None
            if movible:
                tNode = slicer.mrmlScene.AddNewNodeByClass(
                    'vtkMRMLLinearTransformNode', f"T_{nombre}")
                tNode.SetCenterOfTransformation(centro[0], centro[1], centro[2])
                node.SetAndObserveTransformNodeID(tNode.GetID())

            self.estado[nombre] = {
                'node': node,
                'tNode': tNode,
                'centroide0': centro,
                'volumen_cm3': vol,
                'movible': movible,
                'preexistente': (not movible),   # lejos de todo corte
                'activa': True,
                'etiqueta': None,
                'dist_curva': d,
                'color': PALETA[idxColor % len(PALETA)],
            }
            if movible:
                idxColor += 1

        candidatas = [(n, i) for n, i in self.estado.items()
                      if i['movible'] and i['volumen_cm3'] >= VOL_MIN_ANCLA_CM3]
        if not candidatas:
            candidatas = [(n, i) for n, i in self.estado.items() if i['movible']]
        if candidatas:
            nombreAncla = min(candidatas, key=lambda kv: kv[1]['centroide0'][2])[0]
            self._fijar_interno(nombreAncla, avisar=False)
            self.log(f"Ancla detectada: '{nombreAncla}' (centroide mas inferior en S).")
            self.log("Si no es la base craneal + macizo facial, corregilo con fijar('Nombre').")

        for info in self.estado.values():
            self._pintar(info)

        return self.piezas()

    def _fijar_interno(self, nombre, avisar=True):
        info = self._buscar(nombre)
        if info is None:
            return
        anterior = self.ancla
        if anterior is not None and anterior in self.estado and anterior != nombre:
            ant = self.estado[anterior]
            ant['movible'] = not ant['preexistente']
            self._pintar(ant)

        info['movible'] = False
        if info['tNode'] is not None:
            info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
            dn = info['tNode'].GetDisplayNode()
            if dn is not None:
                dn.SetEditorVisibility(False)
        self.ancla = self._nombre_de(info)
        self._pintar(info)
        if avisar:
            self.log(f"Ancla: '{self.ancla}'. Queda fija y gris; el resto se mueve respecto a ella.")

    def fijar(self, nombre):
        """Define manualmente la pieza ancla (base craneal + macizo facial)."""
        self._fijar_interno(nombre, avisar=True)
        return self.tablaDePiezas()

    def fijar_pieza(self, nombre):
        """Marca una pieza como FIJA preexistente (mandibula, vertebras, tubuladura).
        Util si el criterio automatico la dejo movible."""
        info = self._buscar(nombre)
        if info is None:
            return
        info['movible'] = False
        info['preexistente'] = True
        if info['tNode'] is not None:
            info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
            info['node'].SetAndObserveTransformNodeID(None)
            slicer.mrmlScene.RemoveNode(info['tNode'])
            info['tNode'] = None
        self._pintar(info)
        self.log(f"'{nombre}' marcada como fija preexistente (excluida de las metricas).")

    def liberar_pieza(self, nombre):
        """Vuelve movible una pieza marcada como fija por error."""
        info = self._buscar(nombre)
        if info is None:
            return
        if info['tNode'] is None:
            t = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLLinearTransformNode', f"T_{self._nombre_de(info)}")
            c = info['centroide0']
            t.SetCenterOfTransformation(c[0], c[1], c[2])
            info['node'].SetAndObserveTransformNodeID(t.GetID())
            info['tNode'] = t
        info['movible'] = True
        info['preexistente'] = False
        self._pintar(info)
        self.log(f"'{nombre}' liberada: ahora es movible.")

    def piezas(self):
        """Tabla de estado de todas las piezas (al log) + la lista para la
        interfaz."""
        if not self.estado:
            self.log("No hay estado. Corre preparar().")
            return None
        self.log("")
        self.log(f"{'PIEZA':<22}{'ET':<4}{'VOL cm3':>9}{'d.CURVA':>9}  {'ESTADO':<16}{'DESPL mm':>10}")
        self.log("-" * 74)
        for nombre, info in sorted(self.estado.items()):
            if not info['activa']:
                estado = "DESCARTADA"
            elif nombre == self.ancla:
                estado = "ANCLA (fija)"
            elif info['preexistente']:
                estado = "fija (preexist)"
            elif not info['movible']:
                estado = "fija"
            else:
                estado = "movible"
            d = 0.0
            if info['tNode'] is not None:
                d = float(np.linalg.norm(self._centroide_actual(info) - info['centroide0']))
            dc = info['dist_curva']
            dcs = "  n/a" if (dc != dc or dc == float('inf')) else f"{dc:>9.1f}"
            et = info['etiqueta'] or "-"
            self.log(f"{nombre:<22}{et:<4}{info['volumen_cm3']:>9.2f}{dcs}  {estado:<16}{d:>10.1f}")
        activas = [i for i in self.estado.values() if i['activa']]
        desc = [i for i in self.estado.values() if not i['activa']]
        volDesc = sum(i['volumen_cm3'] for i in desc)
        self.log("-" * 74)
        self.log(f"{len(activas)} pieza(s) en el armado, {len(desc)} descartada(s) "
                 f"({volDesc:.2f} cm3 de hueso resecado).")
        self.log("La columna d.CURVA es la distancia al corte: chica = nacio del corte.")
        self.log("")
        return self.tablaDePiezas()

    def tablaDePiezas(self):
        """Lista de dicts lista para la tabla del widget (sin imprimir)."""
        filas = []
        for nombre, info in sorted(self.estado.items()):
            if not info['activa']:
                estado = "descartada"
            elif nombre == self.ancla:
                estado = "base (fija)"
            elif info['preexistente']:
                estado = "fija (no la toco ningun corte)"
            elif not info['movible']:
                estado = "fija"
            else:
                estado = "movible"
            d = 0.0
            if info['tNode'] is not None:
                d = float(np.linalg.norm(
                    self._centroide_actual(info) - info['centroide0']))
            filas.append({
                "nombre": nombre,
                "etiqueta": info['etiqueta'] or "",
                "volumen": info['volumen_cm3'],
                "estado": estado,
                "movible": bool(info['movible'] and info['activa']),
                "activa": bool(info['activa']),
                "esAncla": nombre == self.ancla,
                "preexistente": bool(info['preexistente']),
                "distCurva": info['dist_curva'],
                "desplazamiento": d,
                "color": info['color'],
            })
        return filas

    def etiquetar(self, nombre, letra):
        """Renombra con el vocabulario del quirofano (A, B, C, R, L, F)."""
        info = self._buscar(nombre)
        if info is None:
            return
        info['etiqueta'] = str(letra).upper()
        self.log(f"'{self._nombre_de(info)}' etiquetada como '{info['etiqueta']}'.")

    # ========================================================
    # DESCARTE
    # ========================================================
    def descartar(self, nombre):
        """(G7) Saca la pieza del armado sin borrarla: gris translucida."""
        info = self._buscar(nombre)
        if info is None:
            return False
        if self._nombre_de(info) == self.ancla:
            self.log("No podes descartar el ancla. Primero asigna otra con fijar().")
            return False
        info['activa'] = False
        if info['tNode'] is not None:
            dn = info['tNode'].GetDisplayNode()
            if dn is not None:
                dn.SetEditorVisibility(False)
        self._pintar(info)
        self.log(f"'{self._nombre_de(info)}' descartada ({info['volumen_cm3']:.2f} cm3).")
        return True

    def restaurar(self, nombre):
        """Devuelve al armado una pieza descartada."""
        info = self._buscar(nombre)
        if info is None:
            return False
        info['activa'] = True
        self._pintar(info)
        self.log(f"'{self._nombre_de(info)}' restaurada al armado.")
        return True

    # ========================================================
    # MANIPULACION
    # ========================================================
    def manipular(self, nombre, encender=True):
        """Enciende el gizmo interactivo. Rota alrededor del centroide (G4)."""
        info = self._buscar(nombre)
        if info is None:
            return False
        if info['tNode'] is None or not info['movible']:
            self.log(f"'{self._nombre_de(info)}' es una pieza FIJA. Si tiene que moverse, "
                     f"usa liberar_pieza() o reasigna el ancla con fijar().")
            return False
        if not info['activa']:
            self.log(f"'{self._nombre_de(info)}' esta descartada. Corre restaurar() primero.")
            return False

        info['tNode'].CreateDefaultDisplayNodes()
        dn = info['tNode'].GetDisplayNode()
        dn.SetEditorVisibility(bool(encender))
        dn.SetEditorSliceIntersectionVisibility(False)
        dn.SetEditorTranslationEnabled(True)
        dn.SetEditorRotationEnabled(True)
        dn.SetEditorScalingEnabled(False)
        try:
            dn.UpdateEditorBounds()
        except Exception:
            pass

        if encender:
            for i in self.estado.values():
                if i is not info and i['tNode'] is not None:
                    d = i['tNode'].GetDisplayNode()
                    if d is not None:
                        d.SetEditorVisibility(False)
            self.log(f"Gizmo activo en '{self._nombre_de(info)}'. Arrastra flechas (traslacion) "
                     f"o aros (rotacion) en la vista 3D.")
        else:
            self.log(f"Gizmo apagado en '{self._nombre_de(info)}'.")
        return True

    def apagarGizmos(self):
        """Apaga todas las agarraderas (solo del modulo)."""
        for i in self.estado.values():
            if i['tNode'] is not None:
                d = i['tNode'].GetDisplayNode()
                if d is not None:
                    d.SetEditorVisibility(False)

    def mover(self, nombre, dr=0.0, da=0.0, ds=0.0):
        """
        Traslacion en mm sobre ejes anatomicos RAS.
          dr > 0 -> derecha del paciente | da > 0 -> anterior | ds > 0 -> superior
        """
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        t = vtk.vtkTransform()
        t.Translate(float(dr), float(da), float(ds))
        self._componer(info['tNode'], t.GetMatrix())
        d = float(np.linalg.norm(self._centroide_actual(info) - info['centroide0']))
        self.log(f"'{self._nombre_de(info)}' movida ({dr:+.1f}, {da:+.1f}, {ds:+.1f}) mm. "
                 f"Desplazamiento total: {d:.1f} mm.")
        return True

    def rotar(self, nombre, eje, grados):
        """Rotacion en grados alrededor del CENTROIDE. eje: 'LR', 'AP' o 'SI'."""
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        key = str(eje).upper()
        if key not in EJES:
            self.log(f"Eje '{eje}' invalido. Usa 'LR', 'AP' o 'SI'.")
            return False
        v = EJES[key]
        c = self._centroide_actual(info)
        t = vtk.vtkTransform()
        t.PostMultiply()
        t.Translate(-c[0], -c[1], -c[2])
        t.RotateWXYZ(float(grados), v[0], v[1], v[2])
        t.Translate(c[0], c[1], c[2])
        self._componer(info['tNode'], t.GetMatrix())
        self.log(f"'{self._nombre_de(info)}' rotada {grados:+.1f} grados sobre el eje {key}.")
        return True

    def voltear(self, nombre, eje="SI"):
        """
        (G6) Voltea 180 grados. NO es un espejo: una reflexion tiene determinante
        negativo y el hueso no se puede reflejar. Es la rotacion propia que hace
        el cirujano al pasar una placa de un lado al otro.
          'SI' -> occipital pasa a frontal (swap clasico de escafocefalia)
          'AP' -> cruza izquierda-derecha invirtiendo arriba-abajo
        """
        return self.rotar(nombre, eje, 180.0)

    def resetear(self, nombre=None):
        """Devuelve una pieza (o todas) a su posicion original post-corte."""
        if nombre is None:
            objetivos = list(self.estado.values())
        else:
            info = self._buscar(nombre)
            objetivos = [info] if info else []
        n = 0
        for info in objetivos:
            if info['tNode'] is not None:
                info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
                n += 1
        self.log(f"{n} pieza(s) devuelta(s) a su posicion original.")
        return n

    # ========================================================
    # MEDICION
    # ========================================================
    def _piezas_de_medicion(self):
        """(F2) Piezas que cuentan para el IC: activas y craneales.
        Se excluyen las preexistentes (mandibula, vertebras, tubuladuras)."""
        return [i for i in self.estado.values() if i['activa'] and not i['preexistente']]

    def _bbox_boveda(self, nivel_S):
        """Bounding box de la boveda: solo la porcion por ENCIMA de nivel_S."""
        bounds = None
        for info in self._piezas_de_medicion():
            poly = _poly_mundo(info['node'])
            if poly is None or poly.GetNumberOfPoints() == 0:
                continue
            if nivel_S is not None:
                plano = vtk.vtkPlane()
                plano.SetOrigin(0.0, 0.0, float(nivel_S))
                plano.SetNormal(0.0, 0.0, 1.0)
                clip = vtk.vtkClipPolyData()
                clip.SetInputData(poly)
                clip.SetClipFunction(plano)
                clip.Update()
                poly = clip.GetOutput()
                if poly.GetNumberOfPoints() < 5:
                    continue
            b = list(poly.GetBounds())
            if bounds is None:
                bounds = b
            else:
                bounds = [min(bounds[0], b[0]), max(bounds[1], b[1]),
                          min(bounds[2], b[2]), max(bounds[3], b[3]),
                          min(bounds[4], b[4]), max(bounds[5], b[5])]
        return bounds

    def metricas(self, nivel_S=None):
        """
        (F2) Indice cefalico medido SOLO sobre la boveda, por encima del plano
        de la craneotomia basal, y solo sobre piezas craneales.

        nivel_S : plano de referencia en RAS. None = el detectado en preparar()
                  a partir de la curva cerrada. Pasalo a mano para forzarlo.

        Devuelve un dict para la interfaz, o None.
        """
        if not self.estado:
            self.log("No hay estado. Corre preparar().")
            return None

        if nivel_S is None:
            nivel_S = self.nivelBoveda

        b = self._bbox_boveda(nivel_S)
        if b is None:
            self.log("No hay piezas craneales activas por encima del nivel indicado.")
            return None
        L = b[3] - b[2]
        W = b[1] - b[0]
        H = b[5] - b[4]
        ic = 100.0 * W / L if L > 1e-6 else 0.0

        self.log("")
        self.log(f"--- METRICAS ({CRANIOPLAN_G_VERSION}) ---")
        if nivel_S is None:
            self.log("AVISO: sin nivel de boveda. Estoy midiendo TODO, incluida la cara")
            self.log("  y la base: el IC va a salir falsamente bajo. Pasa nivel_S a mano.")
        else:
            self.log(f"Medido por encima de S = {nivel_S:.1f} mm (plano de craneotomia).")
        excl = [self._nombre_de(i) for i in self.estado.values() if i['preexistente']]
        if excl:
            self.log(f"Excluidas de la medicion: {excl}")

        self.log(f"Largo A-P : {L:.1f} mm")
        self.log(f"Ancho L-R : {W:.1f} mm")
        self.log(f"Alto  S-I : {H:.1f} mm")
        self.log(f"Indice cefalico : {ic:.1f}   (objetivo {self.icObjetivo:.1f})")
        if ic < self.icObjetivo - 1.0:
            falta_L = L - (W * 100.0 / self.icObjetivo)
            self.log(f"  faltan {self.icObjetivo - ic:.1f} puntos de IC.")
            self.log(f"  a ancho constante habria que acortar {falta_L:.1f} mm en A-P,")
            self.log(f"  o ensanchar {(L * self.icObjetivo / 100.0) - W:.1f} mm en transverso")
            self.log("  (el ancho NO se gana con transformaciones rigidas: requiere")
            self.log("   barrel staves o partir la placa en tiras).")
        elif ic > self.icObjetivo + 1.0:
            self.log(f"  te pasaste {ic - self.icObjetivo:.1f} puntos de IC.")
        else:
            self.log("  IC dentro del objetivo.")

        avisoInclinacion = None
        lm = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_LANDMARKS)
        if lm is not None and lm.GetNumberOfControlPoints() >= 2:
            pts = []
            for i in range(lm.GetNumberOfControlPoints()):
                if lm.GetNthControlPointPositionStatus(i) == lm.PositionDefined:
                    p = [0.0, 0.0, 0.0]
                    lm.GetNthControlPointPosition(i, p)
                    pts.append(np.array(p))
            if len(pts) >= 2:
                disp = float(np.std([p[0] for p in pts]))
                if disp > 3.0:
                    self.log(f"  AVISO: la linea media se desvia {disp:.1f} mm en R. "
                             f"La cabeza esta rotada en el CT y el IC sale sesgado.")
                    avisoInclinacion = (
                        "La linea media se desvia %.1f mm. La cabeza esta "
                        "torcida en la tomografia y el indice cefalico sale "
                        "sesgado." % disp)

        self.log("")
        self.log("Transformaciones aplicadas:")
        self.log(f"{'PIEZA':<22}{'ET':<4}{'dR':>7}{'dA':>7}{'dS':>7}{'ROT':>8}  EJE")
        self.log("-" * 74)
        hubo = False
        transformaciones = []
        for nombre, info in sorted(self.estado.items()):
            if not info['activa'] or info['tNode'] is None:
                continue
            m = vtk.vtkMatrix4x4()
            info['tNode'].GetMatrixTransformToParent(m)
            d = self._centroide_actual(info) - info['centroide0']
            R = np.array([[m.GetElement(i, j) for j in range(3)] for i in range(3)])
            traza = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
            ang = float(np.degrees(np.arccos(traza)))
            eje = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
            nn = np.linalg.norm(eje)
            eje = eje / nn if nn > 1e-9 else np.array([0.0, 0.0, 1.0])
            if np.linalg.norm(d) < 0.05 and ang < 0.05:
                continue
            hubo = True
            et = info['etiqueta'] or "-"
            self.log(f"{nombre:<22}{et:<4}{d[0]:>7.1f}{d[1]:>7.1f}{d[2]:>7.1f}{ang:>8.1f}  "
                     f"({eje[0]:.2f}, {eje[1]:.2f}, {eje[2]:.2f})")
            transformaciones.append({
                "nombre": nombre, "etiqueta": info['etiqueta'] or "",
                "dR": float(d[0]), "dA": float(d[1]), "dS": float(d[2]),
                "rotacion": ang, "eje": [float(x) for x in eje],
            })
        if not hubo:
            self.log("  (ninguna pieza movida todavia)")

        desc = [(n, i) for n, i in self.estado.items() if not i['activa']]
        vol = 0.0
        if desc:
            vol = sum(i['volumen_cm3'] for _, i in desc)
            self.log("")
            self.log(f"Hueso resecado: {vol:.2f} cm3 en {len(desc)} pieza(s) "
                     f"({', '.join(n for n, _ in desc)}).")
        self.log("")

        return {"largo": L, "ancho": W, "alto": H, "ic": ic,
                "icObjetivo": self.icObjetivo,
                "nivelS": nivel_S,
                "excluidas": excl,
                "transformaciones": transformaciones,
                "descartadas": [n for n, _ in desc],
                "volumenDescartado": vol,
                "avisoInclinacion": avisoInclinacion}

    def huecos(self, muestreo=MUESTREO_HUECOS):
        """
        (F3) Distancia minima de cada pieza activa a su vecina mas cercana.

        Es el dato que decide conducta en quirofano: que separacion osifica sola,
        cual necesita injerto, donde va una placa. Un 0.0 significa contacto o
        interpenetracion (el modulo no distingue: no hay deteccion de colisiones).

        LIMITE: se mide por muestreo de puntos de superficie, asi que es una cota
        razonable pero no una medida exacta de area de defecto.

        Devuelve una lista [{"pieza", "etiqueta", "vecina", "mm"}], o None.
        """
        activas = [(n, i) for n, i in self.estado.items() if i['activa']]
        if len(activas) < 2:
            self.log("Hacen falta al menos dos piezas activas.")
            return None

        polys = {}
        for n, i in activas:
            p = _poly_mundo(i['node'])
            if p is not None and p.GetNumberOfPoints() > 0:
                polys[n] = p

        locs = {}
        for n, p in polys.items():
            l = vtk.vtkPointLocator()  # noqa: E741 (mismo nombre que en V1.1)
            l.SetDataSet(p)
            l.BuildLocator()
            locs[n] = l

        self.log("")
        self.log(f"--- HUECOS ENTRE PIEZAS ({CRANIOPLAN_G_VERSION}) ---")
        self.log(f"{'PIEZA':<22}{'VECINA MAS CERCANA':<24}{'HUECO mm':>10}")
        self.log("-" * 58)

        resultado = []
        nombres = list(polys.keys())
        for n in nombres:
            pA = polys[n]
            nptsA = pA.GetNumberOfPoints()
            idxs = range(0, nptsA, max(1, int(muestreo)))
            mejor, mejorD = None, float('inf')
            for m in nombres:
                if m == n:
                    continue
                loc = locs[m]
                pB = polys[m]
                dmin = float('inf')
                for k in idxs:
                    q = pA.GetPoint(k)
                    idc = loc.FindClosestPoint(q)
                    if idc < 0:
                        continue
                    r = pB.GetPoint(idc)
                    d = ((q[0] - r[0]) ** 2 + (q[1] - r[1]) ** 2 + (q[2] - r[2]) ** 2) ** 0.5
                    if d < dmin:
                        dmin = d
                        if dmin < 0.01:
                            break
                if dmin < mejorD:
                    mejorD, mejor = dmin, m
            et = self.estado[n]['etiqueta']
            etq = f" ({et})" if et else ""
            self.log(f"{n + etq:<22}{mejor or '-':<24}{mejorD:>10.1f}")
            resultado.append({"pieza": n, "etiqueta": et or "",
                              "vecina": mejor, "mm": float(mejorD)})

        self.log("-" * 58)
        self.log("0.0 = contacto o interpenetracion (no hay deteccion de colisiones).")
        self.log("")
        return resultado

    # ========================================================
    # MOLDE (opcional, ver G5) Y LINEA MEDIA
    # ========================================================
    def molde(self, ic_objetivo=None, nivel_S=None):
        """
        (G5) Casquete elipsoidal de referencia, semitransparente. OPCIONAL y de
        valor limitado: es un elipsoide liso derivado del bbox del paciente, no
        corresponde a ninguna anatomia real y tiende a tapar mas que orientar.
        Preserva el producto L*W (el hueso no se estira) y reparte con el IC
        pedido:  L = sqrt(L0*W0*100/IC) ,  W = L*IC/100

        Devuelve un dict para la interfaz, o None.
        """
        ic_objetivo = float(ic_objetivo if ic_objetivo is not None
                            else self.icObjetivo)
        self.icObjetivo = ic_objetivo
        if nivel_S is None:
            nivel_S = self.nivelBoveda
        b = self._bbox_boveda(nivel_S)
        if b is None:
            self.log("No hay boveda medible. Corre preparar().")
            return None

        L0 = b[3] - b[2]
        W0 = b[1] - b[0]
        H0 = b[5] - b[4]
        L = float(np.sqrt(L0 * W0 * 100.0 / float(ic_objetivo)))
        W = L * float(ic_objetivo) / 100.0
        cx = 0.5 * (b[0] + b[1])
        cy = 0.5 * (b[2] + b[3])
        cz = 0.5 * (b[4] + b[5])

        fuente = vtk.vtkSphereSource()
        fuente.SetRadius(0.5)
        fuente.SetThetaResolution(64)
        fuente.SetPhiResolution(64)
        fuente.Update()

        t = vtk.vtkTransform()
        t.PostMultiply()
        t.Scale(W, L, H0)
        t.Translate(cx, cy, cz)
        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetInputConnection(fuente.GetOutputPort())
        tf.SetTransform(t)
        tf.Update()

        plano = vtk.vtkPlane()
        plano.SetOrigin(cx, cy, cz - 0.15 * H0)
        plano.SetNormal(0.0, 0.0, 1.0)
        clip = vtk.vtkClipPolyData()
        clip.SetInputConnection(tf.GetOutputPort())
        clip.SetClipFunction(plano)
        clip.Update()

        nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
        if nodo is None:
            nodo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', NOMBRE_MOLDE)
            nodo.CreateDefaultDisplayNodes()
        nodo.SetAndObservePolyData(clip.GetOutput())
        dn = nodo.GetDisplayNode()
        dn.SetColor(0.95, 0.95, 0.90)
        dn.SetOpacity(OPACIDAD_MOLDE)
        dn.SetVisibility(True)
        dn.SetBackfaceCulling(False)

        self.log(f"Molde de referencia (IC {ic_objetivo:.1f}).")
        self.log(f"  boveda actual : largo {L0:.1f}  ancho {W0:.1f}  -> IC {100.0*W0/L0:.1f}")
        self.log(f"  objetivo      : largo {L:.1f}  ancho {W:.1f}")

        return {"nodo": nodo, "largoActual": L0, "anchoActual": W0,
                "altoActual": H0, "icActual": 100.0 * W0 / L0,
                "largoObjetivo": L, "anchoObjetivo": W,
                "icObjetivo": ic_objetivo, "nivelS": nivel_S,
                "acortar": L0 - L, "ensanchar": W - W0}

    def verMolde(self, visible=True):
        """Prende o apaga el molde (ocultar_molde() del script = visible False)."""
        nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
        if nodo is not None and nodo.GetDisplayNode():
            nodo.GetDisplayNode().SetVisibility(bool(visible))

    def landmarks(self):
        """Crea los 3 puntos de linea media (nasion, bregma, inion) para
        verificar la inclinacion de la cabeza en el CT."""
        nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_LANDMARKS)
        if nodo is None:
            nodo = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLMarkupsFiducialNode', NOMBRE_LANDMARKS)
            nodo.CreateDefaultDisplayNodes()
            for etiqueta in ("Nasion", "Bregma", "Inion"):
                idx = nodo.AddControlPoint([0.0, 0.0, 0.0])
                nodo.SetNthControlPointLabel(idx, etiqueta)
                nodo.UnsetNthControlPointPosition(idx)
            self.log(f"Nodo '{NOMBRE_LANDMARKS}' creado con 3 puntos sin colocar.")
        self.log("Colocalos sobre la linea media y corre metricas().")
        return nodo

    def ajustar(self, nombre, iteraciones=60):
        """ICP rigido contra el molde. Pose inicial, no resultado final: el
        elipsoide es liso y el ICP desliza tangencialmente."""
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        moldeNode = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
        if moldeNode is None:
            self.log("No hay molde. Corre molde() primero.")
            return False
        origen = _poly_mundo(info['node'])
        destino = moldeNode.GetPolyData()
        if origen is None or destino is None:
            return False
        icp = vtk.vtkIterativeClosestPointTransform()
        icp.SetSource(origen)
        icp.SetTarget(destino)
        icp.GetLandmarkTransform().SetModeToRigidBody()
        icp.SetMaximumNumberOfIterations(int(iteraciones))
        icp.SetCheckMeanDistance(1)
        icp.SetMaximumMeanDistance(0.01)
        icp.StartByMatchingCentroidsOff()
        icp.Modified()
        icp.Update()
        m = vtk.vtkMatrix4x4()
        m.DeepCopy(icp.GetMatrix())
        self._componer(info['tNode'], m)
        self.log(f"'{self._nombre_de(info)}' ajustada al molde por ICP. Revisa a ojo.")
        return True

    # ========================================================
    # EXPORTACION
    # ========================================================
    def exportarSTL(self, carpetaDestino):
        """Exporta el armado final: cada pieza activa con su transformacion
        aplicada (G1: el hardening ocurre aca, sobre copias).
        (exportar_stl_G() del script.)"""
        if not self.estado:
            self.log("No hay estado. Corre preparar().")
            return 0
        if not os.path.isdir(carpetaDestino):
            os.makedirs(carpetaDestino)
        n = 0
        for nombre, info in sorted(self.estado.items()):
            if not info['activa']:
                continue
            poly = _poly_mundo(info['node'])
            if poly is None or poly.GetNumberOfPoints() == 0:
                self.log(f"  {nombre}: malla vacia, no se exporta.")
                continue
            etiqueta = f"_{info['etiqueta']}" if info['etiqueta'] else ""
            ruta = os.path.join(carpetaDestino, f"{nombre}{etiqueta}.stl")
            w = vtk.vtkSTLWriter()
            w.SetFileName(ruta)
            w.SetInputData(poly)
            w.SetFileTypeToBinary()
            w.Write()
            self.log(f"  {nombre} -> {ruta}")
            n += 1
        self.log(f"{n} STL exportado(s) con las transformaciones aplicadas.")
        return n

    def limpiar(self):
        """Borra transforms, molde y landmarks. Los fragmentos NO se tocan.
        (limpiar_G() del script.)"""
        for info in self.estado.values():
            info['node'].SetAndObserveTransformNodeID(None)
            if info['tNode'] is not None and slicer.mrmlScene.IsNodePresent(info['tNode']):
                slicer.mrmlScene.RemoveNode(info['tNode'])
            dn = info['node'].GetDisplayNode()
            if dn is not None:
                dn.SetOpacity(1.0)
        for nombre in (NOMBRE_MOLDE, NOMBRE_LANDMARKS):
            nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
            if nodo is not None:
                slicer.mrmlScene.RemoveNode(nodo)
        self.estado = {}
        self.ancla = None
        self.nivelBoveda = None
        self.log("Bloque G limpiado. Los fragmentos del Bloque F siguen intactos.")
