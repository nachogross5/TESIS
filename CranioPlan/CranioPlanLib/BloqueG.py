# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - BLOQUE G v1  (version libreria, para el modulo)
#              REACOMODAMIENTO DE PIEZAS POST-CORTE
# ============================================================
# Entra DESPUES del Bloque F. Toma los modelos de la carpeta
# 'CranioPlan_Fragmentos' y permite descartar piezas, moverlas y rotarlas para
# rearmar la boveda craneal (remodelacion tipo Melbourne).
#
# Es el script BloqueG_v1.py convertido en libreria: mismo algoritmo, con el
# estado global _ESTADO/_ANCLA metido adentro de una clase para que el modulo
# pueda tirar el estado de un paciente y empezar otro sin reiniciar Slicer.
#
# ------------------------------------------------------------
# DECISIONES DE ARQUITECTURA (leer antes de tocar nada)
# ------------------------------------------------------------
# [G1] LA GEOMETRIA NUNCA SE TOCA. Cada pieza movible recibe su propio
#      vtkMRMLLinearTransformNode. Todo el movimiento vive en esa matriz.
#      Consecuencias: reset instantaneo por pieza, la matriz ES el dato del
#      reporte prequirurgico (traslacion en mm + angulo/eje de rotacion), y se
#      pueden recalcular metricas sin recalcular geometria. El hardening solo
#      ocurre al exportar STL, y sobre una copia.
#
# [G2] MOVIBLE = NACIO DEL CORTE. El Bloque F ya resuelve esto y lo guarda en
#      el nombre:
#         Hueso_N                        -> pieza PREEXISTENTE (el corte no la
#                                           toco). FIJA.
#         Tapa_*, Resto_N, Fragmento_N   -> nacieron del corte. MOVIBLES.
#      No hace falta dilatar el kerf y medir solapamiento: la procedencia ya
#      esta determinada por conectividad, que es mas confiable.
#
# [G3] TOCAR UN CORTE NO ALCANZA PARA SER MOVIBLE. La base craneal y el macizo
#      facial tambien tocan el corte inferior y sin embargo NO se mueven: son
#      el marco de referencia de toda la cirugia. Si la base se moviera, no
#      habria sistema de coordenadas y ninguna metrica significaria nada. Por
#      eso hay un ANCLA: se detecta como la pieza grande con el centroide mas
#      inferior en S. La heuristica no siempre acierta -> fijar('Nombre').
#
#         FIJA (gris claro):  piezas Hueso_N  +  el ANCLA
#         MOVIBLE (color):    todo lo que nacio del corte y no es el ancla
#
# [G4] EL CENTRO DE TRANSFORMACION VA EN EL CENTROIDE DE LA PIEZA. Si no se
#      setea, el gizmo rota alrededor del origen RAS y la pieza sale volando
#      fuera de pantalla en el primer arrastre.
#
# [G5] EL OBJETIVO ES UN MOLDE ELIPSOIDAL POR INDICE CEFALICO, NO LA SIMETRIA.
#      El espejo medio-sagital sirve para deformidades ASIMETRICAS
#      (plagiocefalia, coronal unilateral). En escafocefalia los dos lados
#      estan igual de mal: el espejo no informa nada. Lo que define el objetivo
#      es el IC. El molde se genera preservando el producto L*W (el hueso no se
#      estira) y sirve de casquete de referencia semitransparente, como el
#      casquete blanco de las fotos del Garrahan.
#
# [G6] NO EXISTE 'ESPEJAR' UNA PIEZA. Una reflexion tiene determinante negativo
#      y el hueso no se puede reflejar. Cuando el cirujano pasa una placa de un
#      lado al otro la VOLTEA 180 grados (la tabla interna queda mirando
#      afuera), que es una rotacion propia. Por eso el comando es voltear() y
#      pide el eje.
#
# [G7] DESCARTAR ES UN ESTADO, NO UN BORRADO. La pieza descartada se apaga y
#      queda gris translucida pero sigue en la escena: se necesita para el
#      reporte (volumen de hueso resecado) y porque el cirujano cambia de idea.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (documentados a proposito, no ignorar)
# ------------------------------------------------------------
# - El IC se calcula sobre los ejes RAS del estudio. Si la cabeza esta rotada
#   en el CT, el ancho y el largo salen sesgados. Mitigacion actual: colocar
#   los 3 landmarks de linea media y revisar el aviso de inclinacion que
#   imprime metricas(). Alineacion automatica: pendiente.
# - ajustar() (ICP) contra un elipsoide liso puede deslizar tangencialmente: la
#   superficie no tiene rasgos que anclen la solucion. Sirve como pose inicial,
#   no como resultado final.
# - No hay deteccion de colisiones entre piezas ni de huecos residuales. Se ven
#   a ojo en el 3D. Cuantificacion: pendiente para el reporte.
# ============================================================

import os

import numpy as np
import vtk
import slicer

from .Comun import NOMBRE_CARPETA_SH, NOMBRE_MOLDE, NOMBRE_LANDMARKS

BLOQUE_G_VERSION = "v1.0"

# ==================== CONFIGURACION ====================
IC_OBJETIVO            = 78.0    # indice cefalico objetivo (normocefalia ~76-81)
PREFIJOS_FIJOS         = ("Hueso_",)   # (G2) piezas preexistentes
VOL_MIN_ANCLA_CM3      = 5.0     # el ancla tiene que ser una pieza grande

COLOR_FIJA             = (0.85, 0.85, 0.85)
COLOR_DESCARTADA       = (0.55, 0.55, 0.55)
OPACIDAD_DESCARTADA    = 0.25
OPACIDAD_MOLDE         = 0.22

PALETA = [
    (0.90, 0.30, 0.30), (0.30, 0.65, 0.95), (0.35, 0.80, 0.45),
    (0.95, 0.75, 0.25), (0.75, 0.45, 0.90), (0.30, 0.85, 0.85),
    (0.95, 0.55, 0.30), (0.60, 0.80, 0.30),
]

EJES = {"LR": (1.0, 0.0, 0.0),   # izquierda-derecha (eje R)
        "AP": (0.0, 1.0, 0.0),   # antero-posterior  (eje A)
        "SI": (0.0, 0.0, 1.0)}   # supero-inferior   (eje S)


# ------------------------------------------------------------
# Utilidades de escena (no dependen del estado)
# ------------------------------------------------------------
def _sh():
    return slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)


def _carpeta_fragmentos():
    """Devuelve (shNode, itemId) de la carpeta del Bloque F, o (shNode, None)."""
    shNode = _sh()
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(shNode.GetSceneItemID(), hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if shNode.GetItemName(itemId) == NOMBRE_CARPETA_SH:
            return shNode, itemId
    return shNode, None


def _modelos_de_carpeta(shNode, itemId):
    """Todos los modelos colgando de la carpeta, incluyendo subcarpetas (el
    Bloque F agrupa las tapas en 'Tapas_<curva>')."""
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
                recorrer(hid)   # subcarpeta

    recorrer(itemId)
    return encontrados


def _centroide_y_volumen(modelNode):
    """Centroide (por centro de masa) y volumen en cm3 de la malla original."""
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


def _es_movible_por_nombre(nombre):
    """(G2) Nacio del corte? Todo lo que no sea pieza preexistente."""
    return not any(nombre.startswith(p) for p in PREFIJOS_FIJOS)


def _poly_mundo(modelNode):
    """Malla de la pieza CON su transformacion aplicada (coordenadas mundo).
    No modifica el nodo: devuelve una copia transformada."""
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


def _angulo_y_eje(matriz):
    """Extrae angulo (grados) y eje de la parte rotacional de una 4x4."""
    R = np.array([[matriz.GetElement(i, j) for j in range(3)] for i in range(3)])
    traza = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    ang = np.degrees(np.arccos(traza))
    if ang < 1e-3:
        return 0.0, np.array([0.0, 0.0, 1.0])
    eje = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(eje)
    eje = eje / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])
    return float(ang), eje


# ============================================================
# LA CLASE QUE USA EL MODULO
# ============================================================
class BloqueG:
    """Estado del rearmado de UN paciente."""

    def __init__(self, log=print):
        self.log = log
        self.estado = {}      # nombre -> dict de la pieza
        self.ancla = None     # nombre de la pieza ancla
        self.icObjetivo = IC_OBJETIVO

    # --------------------------------------------------------
    def _buscar(self, nombre):
        """Busca una pieza por nombre o por etiqueta (A, B, R, L, F...)."""
        if nombre in self.estado:
            return self.estado[nombre]
        for n, info in self.estado.items():
            if info['etiqueta'] and info['etiqueta'].upper() == str(nombre).upper():
                return info
        self.log("Bloque G: no encuentro la pieza '%s'." % nombre)
        return None

    def _nombreDe(self, info):
        for n, i in self.estado.items():
            if i is info:
                return n
        return None

    def _pintar(self, info):
        """Aplica color/opacidad segun el estado actual de la pieza."""
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

    def _centroide_actual(self, info):
        """Centroide de la pieza en su posicion ACTUAL (mundo)."""
        if info['tNode'] is None:
            return np.array(info['centroide0'], dtype=float)
        m = vtk.vtkMatrix4x4()
        info['tNode'].GetMatrixTransformToWorld(m)
        c = info['centroide0']
        p = m.MultiplyPoint([c[0], c[1], c[2], 1.0])
        return np.array(p[:3], dtype=float)

    def _componer(self, tNode, matrizMundo):
        """Compone una matriz de mundo POR IZQUIERDA sobre la transformacion
        actual de la pieza: primero lo que ya tenia, despues esto."""
        actual = vtk.vtkMatrix4x4()
        tNode.GetMatrixTransformToParent(actual)
        t = vtk.vtkTransform()
        t.PostMultiply()
        t.SetMatrix(actual)
        t.Concatenate(matrizMundo)
        tNode.SetMatrixTransformToParent(t.GetMatrix())

    # ========================================================
    # PREPARACION
    # ========================================================
    def preparar(self):
        """
        Inventaria los fragmentos del Bloque F, crea un transform por pieza
        movible (G1), detecta el ancla (G3) y colorea todo segun su estado.
        Se puede volver a correr: resetea el estado completo.

        Devuelve la lista de piezas o None.
        """
        shNode, itemId = _carpeta_fragmentos()
        if itemId is None:
            self.log("Bloque G: no encuentro la carpeta '%s'. Falta hacer los "
                     "cortes." % NOMBRE_CARPETA_SH)
            return None
        modelos = _modelos_de_carpeta(shNode, itemId)
        if not modelos:
            self.log("Bloque G: la carpeta existe pero no tiene piezas adentro.")
            return None

        # limpiar transforms previos de una corrida anterior
        for info in self.estado.values():
            t = info.get('tNode')
            if t is not None and slicer.mrmlScene.IsNodePresent(t):
                info['node'].SetAndObserveTransformNodeID(None)
                slicer.mrmlScene.RemoveNode(t)
        self.estado = {}
        self.ancla = None

        self.log("Bloque G %s: preparando %d pieza(s)."
                 % (BLOQUE_G_VERSION, len(modelos)))

        idxColor = 0
        for node in modelos:
            nombre = node.GetName()
            centro, vol = _centroide_y_volumen(node)
            movible = _es_movible_por_nombre(nombre)
            tNode = None
            if movible:
                tNode = slicer.mrmlScene.AddNewNodeByClass(
                    'vtkMRMLLinearTransformNode', "T_%s" % nombre)
                # (G4) el gizmo tiene que girar la pieza sobre si misma
                tNode.SetCenterOfTransformation(centro[0], centro[1], centro[2])
                node.SetAndObserveTransformNodeID(tNode.GetID())
            self.estado[nombre] = {
                'node': node,
                'tNode': tNode,
                'centroide0': centro,
                'volumen_cm3': vol,
                'movible': movible,
                'activa': True,
                'etiqueta': None,
                'color': PALETA[idxColor % len(PALETA)],
            }
            if movible:
                idxColor += 1

        # (G3) Ancla: pieza grande con el centroide mas inferior en S
        candidatas = [(n, i) for n, i in self.estado.items()
                      if i['movible'] and i['volumen_cm3'] >= VOL_MIN_ANCLA_CM3]
        if not candidatas:
            candidatas = [(n, i) for n, i in self.estado.items() if i['movible']]
        if candidatas:
            nombreAncla = min(candidatas, key=lambda kv: kv[1]['centroide0'][2])[0]
            self._fijar_interno(nombreAncla, avisar=False)
            self.log("Bloque G: base detectada automaticamente -> '%s' (la pieza "
                     "grande mas baja). Si no es la base del craneo, cambiala."
                     % nombreAncla)

        for info in self.estado.values():
            self._pintar(info)
        return self.tablaDePiezas()

    def _fijar_interno(self, nombre, avisar=True):
        info = self._buscar(nombre)
        if info is None:
            return
        anterior = self.ancla
        if anterior is not None and anterior in self.estado and anterior != nombre:
            ant = self.estado[anterior]
            ant['movible'] = _es_movible_por_nombre(anterior)
            self._pintar(ant)
        info['movible'] = False
        if info['tNode'] is not None:
            info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
            dn = info['tNode'].GetDisplayNode()
            if dn is not None:
                dn.SetEditorVisibility(False)
        self.ancla = self._nombreDe(info)
        self._pintar(info)
        if avisar:
            self.log("Bloque G: base -> '%s'. Queda fija y en gris; el resto se "
                     "mueve respecto de ella." % self.ancla)

    def fijar(self, nombre):
        """Define manualmente cual es la pieza base (base craneal + macizo
        facial)."""
        self._fijar_interno(nombre, avisar=True)
        return self.tablaDePiezas()

    # --------------------------------------------------------
    def tablaDePiezas(self):
        """Lista de dicts lista para la tabla del widget."""
        filas = []
        for nombre, info in sorted(self.estado.items()):
            if not info['activa']:
                estado = "descartada"
            elif nombre == self.ancla:
                estado = "base (fija)"
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
                "desplazamiento": d,
                "color": info['color'],
            })
        return filas

    def etiquetar(self, nombre, letra):
        """Renombra la pieza con el vocabulario del quirofano (A, B, C, R, L)."""
        info = self._buscar(nombre)
        if info is None:
            return
        info['etiqueta'] = str(letra).upper()
        self.log("Bloque G: '%s' etiquetada como '%s'." % (nombre, info['etiqueta']))

    # ========================================================
    # DESCARTE
    # ========================================================
    def descartar(self, nombre):
        """(G7) Saca la pieza del armado sin borrarla: gris translucida."""
        info = self._buscar(nombre)
        if info is None:
            return False
        if self._nombreDe(info) == self.ancla:
            self.log("Bloque G: no se puede descartar la base. Primero asignar "
                     "otra pieza como base.")
            return False
        info['activa'] = False
        if info['tNode'] is not None:
            dn = info['tNode'].GetDisplayNode()
            if dn is not None:
                dn.SetEditorVisibility(False)
        self._pintar(info)
        self.log("Bloque G: '%s' descartada (%.2f cm3 de hueso resecado). Sigue "
                 "en la escena para el reporte." % (nombre, info['volumen_cm3']))
        return True

    def restaurar(self, nombre):
        """Devuelve al armado una pieza descartada."""
        info = self._buscar(nombre)
        if info is None:
            return False
        info['activa'] = True
        self._pintar(info)
        self.log("Bloque G: '%s' restaurada al armado." % nombre)
        return True

    # ========================================================
    # MANIPULACION
    # ========================================================
    def manipular(self, nombre, encender=True):
        """Enciende el gizmo interactivo (traslacion + rotacion) sobre la
        pieza. El gizmo gira alrededor del centroide (G4)."""
        info = self._buscar(nombre)
        if info is None:
            return False
        if info['tNode'] is None or not info['movible']:
            self.log("Bloque G: '%s' es una pieza FIJA. Si tiene que moverse, "
                     "asignar la base a otra pieza." % nombre)
            return False
        if not info['activa']:
            self.log("Bloque G: '%s' esta descartada. Restaurala primero." % nombre)
            return False

        info['tNode'].CreateDefaultDisplayNodes()
        dn = info['tNode'].GetDisplayNode()
        dn.SetEditorVisibility(bool(encender))
        dn.SetEditorSliceIntersectionVisibility(False)
        dn.SetEditorTranslationEnabled(True)
        dn.SetEditorRotationEnabled(True)
        dn.SetEditorScalingEnabled(False)   # el hueso no se escala
        try:
            dn.UpdateEditorBounds()
        except Exception:
            pass

        if encender:
            # apagar los demas gizmos para no llenar la pantalla
            for n, i in self.estado.items():
                if i is not info and i['tNode'] is not None:
                    d = i['tNode'].GetDisplayNode()
                    if d is not None:
                        d.SetEditorVisibility(False)
            self.log("Bloque G: agarradera activa en '%s'. Arrastrar las flechas "
                     "(mover) o los aros (girar) en la vista 3D." % nombre)
        return True

    def apagarGizmos(self):
        for i in self.estado.values():
            if i['tNode'] is not None:
                d = i['tNode'].GetDisplayNode()
                if d is not None:
                    d.SetEditorVisibility(False)

    def mover(self, nombre, dr=0.0, da=0.0, ds=0.0):
        """
        Traslacion numerica en mm sobre los ejes anatomicos RAS.
          dr > 0 -> hacia la derecha del paciente
          da > 0 -> hacia adelante
          ds > 0 -> hacia arriba
        """
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        t = vtk.vtkTransform()
        t.Translate(float(dr), float(da), float(ds))
        self._componer(info['tNode'], t.GetMatrix())
        d = float(np.linalg.norm(self._centroide_actual(info) - info['centroide0']))
        self.log("Bloque G: '%s' movida (%+.1f, %+.1f, %+.1f) mm. Desplazamiento "
                 "total: %.1f mm." % (nombre, dr, da, ds, d))
        return True

    def rotar(self, nombre, eje, grados):
        """Rotacion en grados alrededor del CENTROIDE de la pieza.
        eje: 'LR' (izq-der), 'AP' (antero-posterior), 'SI' (arriba-abajo)."""
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        key = str(eje).upper()
        if key not in EJES:
            self.log("Bloque G: eje '%s' invalido. Usar 'LR', 'AP' o 'SI'." % eje)
            return False
        v = EJES[key]
        c = self._centroide_actual(info)
        t = vtk.vtkTransform()
        t.PostMultiply()
        t.Translate(-c[0], -c[1], -c[2])
        t.RotateWXYZ(float(grados), v[0], v[1], v[2])
        t.Translate(c[0], c[1], c[2])
        self._componer(info['tNode'], t.GetMatrix())
        self.log("Bloque G: '%s' rotada %+.1f grados sobre el eje %s."
                 % (nombre, grados, key))
        return True

    def voltear(self, nombre, eje="SI"):
        """
        (G6) Voltea la pieza 180 grados. NO es un espejo: una reflexion tiene
        determinante negativo y el hueso no se puede reflejar. Esto es la
        rotacion propia que hace el cirujano cuando pasa una placa de un lado
        al otro.
          eje='SI' -> occipital pasa a frontal (el swap clasico de escafocefalia)
          eje='AP' -> cruza izquierda-derecha invirtiendo arriba-abajo
          eje='LR' -> invierte adelante-atras y arriba-abajo
        """
        ok = self.rotar(nombre, eje, 180.0)
        if ok:
            self.log("Bloque G:   (voltear = rotacion de 180 grados; la tabla "
                     "interna queda hacia afuera, como en la cirugia real)")
        return ok

    def resetear(self, nombre=None):
        """Devuelve una pieza (o todas) a su posicion original post-corte."""
        if nombre is None:
            objetivos = list(self.estado.values())
        else:
            info = self._buscar(nombre)
            objetivos = [info] if info is not None else []
        n = 0
        for info in objetivos:
            if info['tNode'] is not None:
                info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
                n += 1
        self.log("Bloque G: %d pieza(s) devuelta(s) a su posicion original." % n)
        return n

    # ========================================================
    # LINEA MEDIA Y MOLDE OBJETIVO
    # ========================================================
    def landmarks(self):
        """
        Crea (o muestra) el nodo de 3 puntos de linea media: nasion, bregma,
        inion. Sirven para (a) verificar la inclinacion de la cabeza en el CT y
        (b) definir el plano medio-sagital de referencia.
        """
        nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_LANDMARKS)
        if nodo is None:
            nodo = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLMarkupsFiducialNode', NOMBRE_LANDMARKS)
            nodo.CreateDefaultDisplayNodes()
            for etiqueta in ("Nasion", "Bregma", "Inion"):
                idx = nodo.AddControlPoint([0.0, 0.0, 0.0])
                nodo.SetNthControlPointLabel(idx, etiqueta)
                nodo.UnsetNthControlPointPosition(idx)
            self.log("Bloque G: nodo '%s' creado con 3 puntos sin colocar."
                     % NOMBRE_LANDMARKS)
        self.log("Bloque G: colocar Nasion, Bregma e Inion sobre la linea media "
                 "del craneo. Despues, 'Ver medidas' avisa si la cabeza esta "
                 "torcida en la tomografia.")
        return nodo

    def _bbox_activas(self):
        """Caja envolvente RAS de todas las piezas ACTIVAS en su posicion
        actual."""
        bounds = None
        for info in self.estado.values():
            if not info['activa']:
                continue
            poly = _poly_mundo(info['node'])
            if poly is None or poly.GetNumberOfPoints() == 0:
                continue
            b = list(poly.GetBounds())
            if bounds is None:
                bounds = b
            else:
                bounds = [min(bounds[0], b[0]), max(bounds[1], b[1]),
                          min(bounds[2], b[2]), max(bounds[3], b[3]),
                          min(bounds[4], b[4]), max(bounds[5], b[5])]
        return bounds

    def molde(self, ic_objetivo=None):
        """
        (G5) Genera el casquete elipsoidal objetivo, semitransparente.

        Preserva el producto L*W actual (el hueso no se estira) y reparte esa
        superficie con el indice cefalico pedido:

            IC = W/L*100    y    L*W = L0*W0
            =>  L = sqrt(L0*W0*100/IC)  ,  W = L*IC/100

        Es el equivalente digital del casquete blanco de las fotos del
        Garrahan.
        """
        ic_objetivo = float(ic_objetivo if ic_objetivo is not None
                            else self.icObjetivo)
        self.icObjetivo = ic_objetivo
        b = self._bbox_activas()
        if b is None:
            self.log("Bloque G: no hay piezas activas.")
            return None

        L0 = b[3] - b[2]      # largo antero-posterior (eje A)
        W0 = b[1] - b[0]      # ancho izq-der (eje R)
        H0 = b[5] - b[4]      # alto (eje S)
        ic0 = 100.0 * W0 / L0 if L0 > 1e-6 else 0.0

        L = float(np.sqrt(L0 * W0 * 100.0 / ic_objetivo))
        W = L * ic_objetivo / 100.0

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

        # recortar por debajo del nivel de la base para que quede casquete
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

        self.log("Bloque G: molde objetivo generado (IC %.1f)." % ic_objetivo)
        self.log("  actual  : largo %.1f mm  ancho %.1f mm  -> IC %.1f"
                 % (L0, W0, ic0))
        self.log("  objetivo: largo %.1f mm  ancho %.1f mm  -> IC %.1f"
                 % (L, W, ic_objetivo))
        self.log("  hay que ACORTAR %.1f mm de adelante hacia atras y ENSANCHAR "
                 "%.1f mm a lo ancho." % (L0 - L, W - W0))

        return {"nodo": nodo, "largoActual": L0, "anchoActual": W0,
                "altoActual": H0, "icActual": ic0, "largoObjetivo": L,
                "anchoObjetivo": W, "icObjetivo": ic_objetivo,
                "acortar": L0 - L, "ensanchar": W - W0}

    def verMolde(self, visible=True):
        nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
        if nodo is not None and nodo.GetDisplayNode() is not None:
            nodo.GetDisplayNode().SetVisibility(bool(visible))

    def ajustar(self, nombre, iteraciones=60):
        """
        ICP rigido de la pieza contra el molde: la 'pega' a la superficie
        objetivo.

        LIMITE: el elipsoide es liso y no tiene rasgos que anclen la solucion,
        asi que el ICP puede deslizar tangencialmente. Sirve como pose inicial;
        el ajuste fino lo hace el cirujano con la agarradera.
        """
        info = self._buscar(nombre)
        if info is None or info['tNode'] is None:
            return False
        moldeNode = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
        if moldeNode is None:
            self.log("Bloque G: no hay molde. Generalo primero.")
            return False
        origen = _poly_mundo(info['node'])
        destino = moldeNode.GetPolyData()
        if origen is None or destino is None:
            return False

        icp = vtk.vtkIterativeClosestPointTransform()
        icp.SetSource(origen)
        icp.SetTarget(destino)
        icp.GetLandmarkTransform().SetModeToRigidBody()   # sin escala ni reflexion
        icp.SetMaximumNumberOfIterations(int(iteraciones))
        icp.SetCheckMeanDistance(1)
        icp.SetMaximumMeanDistance(0.01)
        icp.StartByMatchingCentroidsOff()
        icp.Modified()
        icp.Update()

        m = vtk.vtkMatrix4x4()
        m.DeepCopy(icp.GetMatrix())
        self._componer(info['tNode'], m)
        self.log("Bloque G: '%s' acercada al molde. Revisar a ojo: sobre una "
                 "superficie lisa el ajuste automatico puede deslizar." % nombre)
        return True

    # ========================================================
    # METRICAS
    # ========================================================
    def metricas(self):
        """
        Indice cefalico actual vs objetivo + transformacion aplicada a cada
        pieza. Es la salida que despues alimenta el reporte prequirurgico.
        """
        if not self.estado:
            self.log("Bloque G: no hay piezas preparadas.")
            return None
        b = self._bbox_activas()
        if b is None:
            self.log("Bloque G: no hay piezas activas.")
            return None

        L = b[3] - b[2]
        W = b[1] - b[0]
        H = b[5] - b[4]
        ic = 100.0 * W / L if L > 1e-6 else 0.0

        self.log("")
        self.log("--- MEDIDAS ---")
        self.log("Largo (adelante-atras) : %.1f mm" % L)
        self.log("Ancho (lado a lado)    : %.1f mm" % W)
        self.log("Alto                   : %.1f mm" % H)
        self.log("Indice cefalico actual : %.1f   (objetivo %.1f)"
                 % (ic, self.icObjetivo))

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
                    avisoInclinacion = (
                        "La linea media se desvia %.1f mm. La cabeza esta "
                        "torcida en la tomografia y el indice cefalico sale "
                        "sesgado." % disp)
                    self.log("  AVISO: " + avisoInclinacion)

        transformaciones = []
        self.log("")
        self.log("Movimientos aplicados:")
        for nombre, info in sorted(self.estado.items()):
            if not info['activa'] or info['tNode'] is None:
                continue
            m = vtk.vtkMatrix4x4()
            info['tNode'].GetMatrixTransformToParent(m)
            d = self._centroide_actual(info) - info['centroide0']
            ang, eje = _angulo_y_eje(m)
            if np.linalg.norm(d) < 0.05 and ang < 0.05:
                continue
            transformaciones.append({
                "nombre": nombre, "etiqueta": info['etiqueta'] or "",
                "dR": float(d[0]), "dA": float(d[1]), "dS": float(d[2]),
                "rotacion": ang, "eje": [float(x) for x in eje],
            })
            self.log("  %-22s  dR %+6.1f  dA %+6.1f  dS %+6.1f   giro %5.1f deg"
                     % (nombre, d[0], d[1], d[2], ang))
        if not transformaciones:
            self.log("  (ninguna pieza se movio todavia)")

        descartadas = [(n, i) for n, i in self.estado.items() if not i['activa']]
        volDescartado = sum(i['volumen_cm3'] for _, i in descartadas)
        if descartadas:
            self.log("")
            self.log("Hueso resecado: %.2f cm3 en %d pieza(s): %s"
                     % (volDescartado, len(descartadas),
                        ", ".join(n for n, _ in descartadas)))

        return {"largo": L, "ancho": W, "alto": H, "ic": ic,
                "icObjetivo": self.icObjetivo,
                "transformaciones": transformaciones,
                "descartadas": [n for n, _ in descartadas],
                "volumenDescartado": volDescartado,
                "avisoInclinacion": avisoInclinacion}

    # ========================================================
    # EXPORTACION
    # ========================================================
    def exportarSTL(self, carpetaDestino):
        """
        Exporta el armado FINAL: cada pieza activa con su transformacion
        aplicada (G1: el hardening ocurre aca, sobre una copia; los nodos de la
        escena quedan intactos y se pueden seguir moviendo).
        """
        if not self.estado:
            self.log("Bloque G: no hay piezas preparadas.")
            return 0
        if not os.path.isdir(carpetaDestino):
            os.makedirs(carpetaDestino)
        n = 0
        for nombre, info in sorted(self.estado.items()):
            if not info['activa']:
                continue
            poly = _poly_mundo(info['node'])
            if poly is None or poly.GetNumberOfPoints() == 0:
                self.log("Bloque G:   %s: malla vacia, no se exporta." % nombre)
                continue
            etiqueta = "_%s" % info['etiqueta'] if info['etiqueta'] else ""
            ruta = os.path.join(carpetaDestino, "%s%s.stl" % (nombre, etiqueta))
            w = vtk.vtkSTLWriter()
            w.SetFileName(ruta)
            w.SetInputData(poly)
            w.SetFileTypeToBinary()
            w.Write()
            self.log("Bloque G:   %s -> %s" % (nombre, ruta))
            n += 1
        self.log("Bloque G: %d STL exportado(s) con los movimientos aplicados." % n)
        return n

    def limpiar(self):
        """Borra los transforms, el molde y los landmarks. Los fragmentos del
        Bloque F NO se tocan."""
        for info in self.estado.values():
            try:
                info['node'].SetAndObserveTransformNodeID(None)
                if (info['tNode'] is not None
                        and slicer.mrmlScene.IsNodePresent(info['tNode'])):
                    slicer.mrmlScene.RemoveNode(info['tNode'])
                dn = info['node'].GetDisplayNode()
                if dn is not None:
                    dn.SetOpacity(1.0)
            except Exception:
                pass
        for nombre in (NOMBRE_MOLDE, NOMBRE_LANDMARKS):
            nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
            if nodo is not None:
                slicer.mrmlScene.RemoveNode(nodo)
        self.estado = {}
        self.ancla = None
        self.log("Bloque G: limpiado. Los fragmentos siguen intactos.")
