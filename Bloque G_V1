# ============================================================
# CRANIOPLAN - BLOQUE G v1 (REACOMODAMIENTO DE PIEZAS POST-CORTE)
#
# Entra DESPUES del Bloque F (cortar()). Toma los modelos de la carpeta
# 'CranioPlan_Fragmentos' y permite descartar piezas, moverlas y rotarlas
# para rearmar la boveda craneal (remodelacion tipo Melbourne).
#
# ------------------------------------------------------------
# DECISIONES DE ARQUITECTURA (leer antes de tocar nada)
# ------------------------------------------------------------
# [G1] LA GEOMETRIA NUNCA SE TOCA. Cada pieza movible recibe su propio
#      vtkMRMLLinearTransformNode. Todo el movimiento vive en esa matriz.
#      Consecuencias: reset instantaneo por pieza, la matriz ES el dato del
#      reporte prequirurgico (traslacion en mm + angulo/eje de rotacion), y
#      se pueden recalcular metricas sin recalcular geometria. El hardening
#      solo ocurre al exportar STL, y sobre una copia.
#
# [G2] MOVIBLE = NACIO DEL CORTE. El Bloque F ya resuelve esto y lo guarda en
#      el nombre, solo que nadie lo estaba usando:
#         Hueso_N                        -> isla PREEXISTENTE (creadaPorCorte=False).
#                                           No toco ninguna linea trazada. FIJA.
#         Tapa_*, Resto_N, Fragmento_N   -> nacieron del corte. MOVIBLES.
#      No hace falta dilatar el kerf y medir solapamiento: la procedencia ya
#      esta determinada por conectividad, que es mas confiable.
#
# [G3] TOCAR UN CORTE NO ALCANZA PARA SER MOVIBLE. La base craneal y el macizo
#      facial tambien tocan el corte inferior y sin embargo NO se mueven: son
#      el marco de referencia de toda la cirugia. Si la base se moviera, no
#      habria sistema de coordenadas y ninguna metrica significaria nada.
#      Por eso hay un ANCLA: se detecta como la pieza grande con el centroide
#      mas inferior en S. La heuristica no siempre acierta -> fijar('Nombre').
#
#      Regla final:
#         FIJA (gris claro):  islas Hueso_N  +  el ANCLA
#         MOVIBLE (color):    todo lo que nacio del corte y no es el ancla
#
# [G4] EL CENTRO DE TRANSFORMACION VA EN EL CENTROIDE DE LA PIEZA. Si no se
#      setea, el gizmo rota alrededor del origen RAS y la pieza sale volando
#      fuera de pantalla en el primer arrastre. Con el centroide, el gizmo
#      gira la pieza sobre si misma, que es lo que el cirujano espera.
#
# [G5] EL OBJETIVO ES UN MOLDE ELIPSOIDAL POR INDICE CEFALICO, NO LA SIMETRIA.
#      El espejo medio-sagital sirve para deformidades ASIMETRICAS
#      (plagiocefalia, coronal unilateral). En escafocefalia (caso Fleitas,
#      IC 61.2) los dos lados estan igual de mal: el espejo no informa nada.
#      Lo que define el objetivo es el IC: llevarlo de 61.2 a ~78. El molde se
#      genera preservando el producto L*W (el hueso no se estira) y sirve de
#      casquete de referencia semitransparente, como el casquete blanco de las
#      fotos del Garrahan.
#
# [G6] NO EXISTE 'ESPEJAR' UNA PIEZA. Una reflexion tiene determinante negativo
#      y el hueso no se puede reflejar. Cuando el cirujano pasa una placa de un
#      lado al otro la VOLTEA 180 grados (la tabla interna queda mirando
#      afuera), que es una rotacion propia. Por eso el comando es voltear() y
#      pide el eje: 180 sobre SI lleva occipital a frontal (el swap clasico de
#      escafocefalia); 180 sobre AP cruza izq-der invirtiendo arriba-abajo.
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
#   los 3 landmarks de linea media con landmarks() y revisar el aviso de
#   inclinacion que imprime metricas(). Alineacion automatica: pendiente.
# - ajustar() (ICP) contra un elipsoide liso puede deslizar tangencialmente:
#   la superficie no tiene rasgos que anclen la solucion. Sirve como pose
#   inicial, no como resultado final. El cirujano ajusta despues a mano.
# - No hay deteccion de colisiones entre piezas ni de huecos residuales.
#   Se ven a ojo en el 3D. Cuantificacion: pendiente para el reporte (Paso 5).
#
# FLUJO:
#   Bloque A+B -> confirmar_craneo() -> enviar_a_planner()
#   -> curvas -> Bloque F cortar()
#   -> cargar este archivo -> preparar() -> landmarks() -> molde()
#   -> descartar(...) / manipular(...) / mover(...) / rotar(...) / voltear(...)
#   -> metricas() -> exportar_stl_G('C:/salida')
#
# REQUISITO: haber corrido cortar() del Bloque F en esta misma escena.
# ============================================================

import os
import numpy as np
import vtk

CRANIOPLAN_G_VERSION = "Bloque G v1.0"

# ==================== CONFIGURACION ====================

NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"   # carpeta que deja el Bloque F
NOMBRE_MOLDE           = "CranioPlan_Molde"
NOMBRE_LANDMARKS       = "CranioPlan_LineaMedia"

IC_OBJETIVO            = 78.0    # indice cefalico objetivo (normocefalia ~76-81)
PREFIJOS_FIJOS         = ("Hueso_",)   # (G2) islas preexistentes: nunca movibles
VOL_MIN_ANCLA_CM3      = 5.0     # el ancla tiene que ser una pieza grande

COLOR_FIJA             = (0.85, 0.85, 0.85)
COLOR_DESCARTADA       = (0.55, 0.55, 0.55)
OPACIDAD_DESCARTADA    = 0.25
OPACIDAD_MOLDE         = 0.22

# paleta para piezas movibles (se recorre ciclicamente)
PALETA = [
    (0.90, 0.30, 0.30), (0.30, 0.65, 0.95), (0.35, 0.80, 0.45),
    (0.95, 0.75, 0.25), (0.75, 0.45, 0.90), (0.30, 0.85, 0.85),
    (0.95, 0.55, 0.30), (0.60, 0.80, 0.30),
]

EJES = {"LR": (1.0, 0.0, 0.0),   # izquierda-derecha (eje R)
        "AP": (0.0, 1.0, 0.0),   # antero-posterior  (eje A)
        "SI": (0.0, 0.0, 1.0)}   # supero-inferior   (eje S)

# Estado global del modulo. nombre -> dict con:
#   node, tNode, centroide0 (np3), volumen_cm3, movible, activa, etiqueta, color
_ESTADO = {}
_ANCLA = [None]      # lista de 1 para poder reasignar desde funciones


# ============================================================
# PARTE 1 - UTILIDADES
# ============================================================

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
    """Todos los vtkMRMLModelNode colgando de la carpeta, incluyendo
    subcarpetas (el Bloque F agrupa las tapas en 'Tapas_<curva>')."""
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
    """(G2) Nacio del corte? Todo lo que no sea isla preexistente."""
    return not any(nombre.startswith(p) for p in PREFIJOS_FIJOS)


def _poly_mundo(modelNode):
    """Polydata de la pieza CON su transformacion aplicada (coordenadas mundo).
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


def _centroide_actual(info):
    """Centroide de la pieza en su posicion ACTUAL (mundo)."""
    m = vtk.vtkMatrix4x4()
    info['tNode'].GetMatrixTransformToWorld(m)
    c = info['centroide0']
    p = m.MultiplyPoint([c[0], c[1], c[2], 1.0])
    return np.array(p[:3], dtype=float)


def _componer(tNode, matrizMundo):
    """Compone una matriz de mundo POR IZQUIERDA sobre la transformacion
    actual de la pieza. Es decir: primero lo que ya tenia, despues esto."""
    actual = vtk.vtkMatrix4x4()
    tNode.GetMatrixTransformToParent(actual)
    t = vtk.vtkTransform()
    t.PostMultiply()
    t.SetMatrix(actual)
    t.Concatenate(matrizMundo)
    tNode.SetMatrixTransformToParent(t.GetMatrix())


def _buscar(nombre):
    """Busca una pieza por nombre o por etiqueta (A, B, R, L, F...)."""
    if nombre in _ESTADO:
        return _ESTADO[nombre]
    for n, info in _ESTADO.items():
        if info['etiqueta'] and info['etiqueta'].upper() == str(nombre).upper():
            return info
    print(f"  no encuentro la pieza '{nombre}'. Corre piezas() para ver la lista.")
    return None


def _pintar(info):
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


# ============================================================
# PARTE 2 - PREPARACION
# ============================================================

def preparar():
    """
    Inventaria los fragmentos del Bloque F, crea un transform por pieza movible
    (G1), detecta el ancla (G3) y colorea todo segun su estado.
    Se puede volver a correr: resetea el estado completo.
    """
    global _ESTADO
    shNode, itemId = _carpeta_fragmentos()
    if itemId is None:
        print(f"ERROR: no encuentro la carpeta '{NOMBRE_CARPETA_SH}'. Corre cortar() del Bloque F.")
        return

    modelos = _modelos_de_carpeta(shNode, itemId)
    if not modelos:
        print("ERROR: la carpeta existe pero no tiene modelos adentro.")
        return

    # limpiar transforms previos de una corrida anterior
    for info in _ESTADO.values():
        t = info.get('tNode')
        if t is not None and slicer.mrmlScene.IsNodePresent(t):
            info['node'].SetAndObserveTransformNodeID(None)
            slicer.mrmlScene.RemoveNode(t)
    _ESTADO = {}
    _ANCLA[0] = None

    print(f"--- {CRANIOPLAN_G_VERSION} : preparar() ---")
    idxColor = 0
    for node in modelos:
        nombre = node.GetName()
        centro, vol = _centroide_y_volumen(node)
        movible = _es_movible_por_nombre(nombre)

        tNode = None
        if movible:
            tNode = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLLinearTransformNode', f"T_{nombre}")
            # (G4) el gizmo tiene que girar la pieza sobre si misma
            tNode.SetCenterOfTransformation(centro[0], centro[1], centro[2])
            node.SetAndObserveTransformNodeID(tNode.GetID())

        _ESTADO[nombre] = {
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
    candidatas = [(n, i) for n, i in _ESTADO.items()
                  if i['movible'] and i['volumen_cm3'] >= VOL_MIN_ANCLA_CM3]
    if not candidatas:
        candidatas = [(n, i) for n, i in _ESTADO.items() if i['movible']]
    if candidatas:
        nombreAncla = min(candidatas, key=lambda kv: kv[1]['centroide0'][2])[0]
        _fijar_interno(nombreAncla, avisar=False)
        print(f"Ancla detectada automaticamente: '{nombreAncla}' "
              f"(centroide mas inferior en S). Si no es la base craneal, "
              f"corregilo con fijar('Nombre').")

    for info in _ESTADO.values():
        _pintar(info)

    piezas()


def _fijar_interno(nombre, avisar=True):
    info = _buscar(nombre)
    if info is None:
        return
    anterior = _ANCLA[0]
    if anterior is not None and anterior in _ESTADO and anterior != nombre:
        ant = _ESTADO[anterior]
        ant['movible'] = _es_movible_por_nombre(anterior)
        _pintar(ant)

    info['movible'] = False
    if info['tNode'] is not None:
        info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
        dn = info['tNode'].GetDisplayNode()
        if dn is not None:
            dn.SetEditorVisibility(False)
    _ANCLA[0] = [n for n, i in _ESTADO.items() if i is info][0]
    _pintar(info)
    if avisar:
        print(f"Ancla: '{_ANCLA[0]}'. Queda fija y en gris; el resto se mueve respecto a ella.")


def fijar(nombre):
    """Define manualmente cual es la pieza ancla (base craneal + macizo facial)."""
    _fijar_interno(nombre, avisar=True)


def piezas():
    """Tabla de estado de todas las piezas."""
    if not _ESTADO:
        print("No hay estado. Corre preparar().")
        return
    print("")
    print(f"{'PIEZA':<24}{'ET':<5}{'VOL cm3':>10}  {'ESTADO':<14}{'DESPLAZ mm':>12}")
    print("-" * 70)
    for nombre, info in sorted(_ESTADO.items()):
        if not info['activa']:
            estado = "DESCARTADA"
        elif nombre == _ANCLA[0]:
            estado = "ANCLA (fija)"
        elif not info['movible']:
            estado = "fija"
        else:
            estado = "movible"
        d = 0.0
        if info['tNode'] is not None:
            d = float(np.linalg.norm(_centroide_actual(info) - info['centroide0']))
        et = info['etiqueta'] or "-"
        print(f"{nombre:<24}{et:<5}{info['volumen_cm3']:>10.2f}  {estado:<14}{d:>12.1f}")
    activas = [i for i in _ESTADO.values() if i['activa']]
    desc = [i for i in _ESTADO.values() if not i['activa']]
    volDesc = sum(i['volumen_cm3'] for i in desc)
    print("-" * 70)
    print(f"{len(activas)} pieza(s) en el armado, {len(desc)} descartada(s) "
          f"({volDesc:.2f} cm3 de hueso resecado).")
    print("")


def etiquetar(nombre, letra):
    """Renombra la pieza con el vocabulario del quirofano (A, B, C, R, L, F)."""
    info = _buscar(nombre)
    if info is None:
        return
    info['etiqueta'] = str(letra).upper()
    print(f"'{nombre}' etiquetada como '{info['etiqueta']}'. "
          f"Ahora la podes referenciar por la letra.")


# ============================================================
# PARTE 3 - DESCARTE
# ============================================================

def descartar(nombre):
    """(G7) Saca la pieza del armado sin borrarla: gris translucida."""
    info = _buscar(nombre)
    if info is None:
        return
    if nombre == _ANCLA[0] or info is _ESTADO.get(_ANCLA[0]):
        print("No podes descartar el ancla. Primero asigna otra con fijar().")
        return
    info['activa'] = False
    if info['tNode'] is not None:
        dn = info['tNode'].GetDisplayNode()
        if dn is not None:
            dn.SetEditorVisibility(False)
    _pintar(info)
    print(f"'{nombre}' descartada ({info['volumen_cm3']:.2f} cm3). "
          f"Sigue en la escena para el reporte; restaurar('{nombre}') la devuelve.")


def restaurar(nombre):
    """Devuelve al armado una pieza descartada."""
    info = _buscar(nombre)
    if info is None:
        return
    info['activa'] = True
    _pintar(info)
    print(f"'{nombre}' restaurada al armado.")


# ============================================================
# PARTE 4 - MANIPULACION
# ============================================================

def manipular(nombre, encender=True):
    """Enciende el gizmo interactivo (traslacion + rotacion) sobre la pieza.
    El gizmo gira alrededor del centroide (G4), no del origen RAS."""
    info = _buscar(nombre)
    if info is None:
        return
    if info['tNode'] is None or not info['movible']:
        print(f"'{nombre}' es una pieza FIJA. Si tiene que moverse, "
              f"reasigna el ancla con fijar() sobre otra pieza.")
        return
    if not info['activa']:
        print(f"'{nombre}' esta descartada. Corre restaurar('{nombre}') primero.")
        return

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

    # apagar los demas gizmos para no llenar la pantalla
    if encender:
        for n, i in _ESTADO.items():
            if i is not info and i['tNode'] is not None:
                d = i['tNode'].GetDisplayNode()
                if d is not None:
                    d.SetEditorVisibility(False)
        print(f"Gizmo activo en '{nombre}'. Arrastra las flechas (traslacion) "
              f"o los aros (rotacion) en la vista 3D.")
    else:
        print(f"Gizmo apagado en '{nombre}'.")


def mover(nombre, dr=0.0, da=0.0, ds=0.0):
    """
    Traslacion numerica en mm sobre los ejes anatomicos RAS.
      dr > 0 -> hacia la derecha del paciente
      da > 0 -> hacia anterior
      ds > 0 -> hacia superior
    Ej: acercar la tapa occipital al segmento medio -> mover('C', da=8)
    """
    info = _buscar(nombre)
    if info is None or info['tNode'] is None:
        return
    t = vtk.vtkTransform()
    t.Translate(float(dr), float(da), float(ds))
    _componer(info['tNode'], t.GetMatrix())
    d = float(np.linalg.norm(_centroide_actual(info) - info['centroide0']))
    print(f"'{nombre}' movida ({dr:+.1f}, {da:+.1f}, {ds:+.1f}) mm. "
          f"Desplazamiento total desde el origen: {d:.1f} mm.")


def rotar(nombre, eje, grados):
    """
    Rotacion en grados alrededor del CENTROIDE de la pieza.
      eje: 'LR' (izq-der), 'AP' (antero-posterior), 'SI' (supero-inferior)
    """
    info = _buscar(nombre)
    if info is None or info['tNode'] is None:
        return
    key = str(eje).upper()
    if key not in EJES:
        print(f"Eje '{eje}' invalido. Usa 'LR', 'AP' o 'SI'.")
        return
    v = EJES[key]
    c = _centroide_actual(info)

    t = vtk.vtkTransform()
    t.PostMultiply()
    t.Translate(-c[0], -c[1], -c[2])
    t.RotateWXYZ(float(grados), v[0], v[1], v[2])
    t.Translate(c[0], c[1], c[2])
    _componer(info['tNode'], t.GetMatrix())
    print(f"'{nombre}' rotada {grados:+.1f} grados sobre el eje {key} "
          f"(alrededor de su centroide).")


def voltear(nombre, eje="SI"):
    """
    (G6) Voltea la pieza 180 grados. NO es un espejo: una reflexion tiene
    determinante negativo y el hueso no se puede reflejar. Esto es la rotacion
    propia que hace el cirujano cuando pasa una placa de un lado al otro.
      eje='SI' -> occipital pasa a frontal (swap clasico de escafocefalia)
      eje='AP' -> cruza izquierda-derecha invirtiendo arriba-abajo
      eje='LR' -> invierte anterior-posterior y arriba-abajo
    La tabla interna queda mirando hacia afuera, igual que en el quirofano.
    """
    rotar(nombre, eje, 180.0)
    print("  (voltear = rotacion propia de 180 grados; la tabla interna "
          "queda hacia afuera, como en la cirugia real)")


def resetear(nombre=None):
    """Devuelve una pieza (o todas) a su posicion original post-corte."""
    objetivos = _ESTADO.items() if nombre is None else \
        [(n, i) for n, i in _ESTADO.items() if i is _buscar(nombre)]
    n = 0
    for _, info in objetivos:
        if info['tNode'] is not None:
            info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
            n += 1
    print(f"{n} pieza(s) devuelta(s) a su posicion original.")


# ============================================================
# PARTE 5 - LINEA MEDIA Y MOLDE OBJETIVO
# ============================================================

def landmarks():
    """
    Crea (o muestra) el nodo de 3 puntos de linea media: nasion, bregma, inion.
    Sirven para (a) verificar la inclinacion de la cabeza en el CT y
    (b) definir el plano medio-sagital de referencia.
    Colocalos en la vista 3D sobre el ancla / macizo facial.
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
        print(f"Nodo '{NOMBRE_LANDMARKS}' creado con 3 puntos sin colocar.")
    print("Coloca Nasion, Bregma e Inion sobre la linea media del craneo "
          "(modulo Markups -> click en cada punto -> ubicar en la vista 3D).")
    print("Despues corre metricas(): avisa si la cabeza esta inclinada en el CT.")
    return nodo


def _bbox_activas():
    """Bounding box RAS de todas las piezas ACTIVAS en su posicion actual."""
    bounds = None
    for info in _ESTADO.values():
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


def molde(ic_objetivo=IC_OBJETIVO):
    """
    (G5) Genera el casquete elipsoidal objetivo, semitransparente.
    Preserva el producto L*W actual (el hueso no se estira) y reparte esa
    superficie con el indice cefalico pedido:
        IC = W/L*100    y    L*W = L0*W0
        =>  L = sqrt(L0*W0*100/IC)  ,  W = L*IC/100
    Es el equivalente digital del casquete blanco de las fotos del Garrahan.
    """
    b = _bbox_activas()
    if b is None:
        print("No hay piezas activas. Corre preparar().")
        return

    L0 = b[3] - b[2]      # largo antero-posterior (eje A)
    W0 = b[1] - b[0]      # ancho izq-der (eje R)
    H0 = b[5] - b[4]      # alto (eje S)
    ic0 = 100.0 * W0 / L0 if L0 > 1e-6 else 0.0

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

    # recortar por debajo del nivel de la base para que quede casquete, no huevo
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

    print(f"Molde objetivo generado (IC {ic_objetivo:.1f}).")
    print(f"  actual : largo {L0:.1f} mm  ancho {W0:.1f} mm  -> IC {ic0:.1f}")
    print(f"  objetivo: largo {L:.1f} mm  ancho {W:.1f} mm  -> IC {ic_objetivo:.1f}")
    print(f"  hay que ACORTAR {L0 - L:.1f} mm en antero-posterior y "
          f"ENSANCHAR {W - W0:.1f} mm en transverso.")
    print("  Acerca las piezas anteriores y posteriores entre si: cada mm que "
          "acortas en A-P es IC que ganas.")
    return nodo


def ajustar(nombre, iteraciones=60):
    """
    ICP rigido de la pieza contra el molde: la 'pega' a la superficie objetivo.
    LIMITE: el elipsoide es liso y no tiene rasgos que anclen la solucion, asi
    que el ICP puede deslizar tangencialmente. Sirve como pose inicial; el
    ajuste fino lo hace el cirujano con el gizmo.
    """
    info = _buscar(nombre)
    if info is None or info['tNode'] is None:
        return
    moldeNode = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
    if moldeNode is None:
        print("No hay molde. Corre molde() primero.")
        return

    origen = _poly_mundo(info['node'])
    destino = moldeNode.GetPolyData()
    if origen is None or destino is None:
        return

    icp = vtk.vtkIterativeClosestPointTransform()
    icp.SetSource(origen)
    icp.SetTarget(destino)
    icp.GetLandmarkTransform().SetModeToRigidBody()   # sin escala, sin reflexion
    icp.SetMaximumNumberOfIterations(int(iteraciones))
    icp.SetCheckMeanDistance(1)
    icp.SetMaximumMeanDistance(0.01)
    icp.StartByMatchingCentroidsOff()
    icp.Modified()
    icp.Update()

    m = vtk.vtkMatrix4x4()
    m.DeepCopy(icp.GetMatrix())
    _componer(info['tNode'], m)
    print(f"'{nombre}' ajustada al molde por ICP ({icp.GetNumberOfIterations()} it.). "
          f"Revisa a ojo: el ICP puede deslizar sobre una superficie lisa.")


# ============================================================
# PARTE 6 - METRICAS
# ============================================================

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


def metricas():
    """
    Indice cefalico actual vs objetivo + transformacion aplicada a cada pieza.
    Esta es la salida que despues alimenta el reporte prequirurgico (Paso 5).
    """
    if not _ESTADO:
        print("No hay estado. Corre preparar().")
        return

    b = _bbox_activas()
    if b is None:
        print("No hay piezas activas.")
        return
    L = b[3] - b[2]
    W = b[1] - b[0]
    H = b[5] - b[4]
    ic = 100.0 * W / L if L > 1e-6 else 0.0

    print("")
    print(f"--- METRICAS ({CRANIOPLAN_G_VERSION}) ---")
    print(f"Largo A-P : {L:.1f} mm")
    print(f"Ancho L-R : {W:.1f} mm")
    print(f"Alto  S-I : {H:.1f} mm")
    print(f"Indice cefalico actual : {ic:.1f}   (objetivo {IC_OBJETIVO:.1f})")
    if ic < IC_OBJETIVO - 1.0:
        print(f"  faltan {IC_OBJETIVO - ic:.1f} puntos de IC: sigue acortando en A-P.")
    elif ic > IC_OBJETIVO + 1.0:
        print(f"  te pasaste {ic - IC_OBJETIVO:.1f} puntos de IC.")
    else:
        print("  IC dentro del objetivo.")

    # aviso de inclinacion usando los landmarks de linea media
    lm = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_LANDMARKS)
    if lm is not None and lm.GetNumberOfControlPoints() >= 2:
        pts = []
        for i in range(lm.GetNumberOfControlPoints()):
            if lm.GetNthControlPointPositionStatus(i) == lm.PositionDefined:
                p = [0.0, 0.0, 0.0]
                lm.GetNthControlPointPosition(i, p)
                pts.append(np.array(p))
        if len(pts) >= 2:
            desvio = float(np.max([abs(p[0]) for p in pts]) -
                           np.min([abs(p[0]) for p in pts]))
            disp = float(np.std([p[0] for p in pts]))
            if disp > 3.0:
                print(f"  AVISO: la linea media se desvia {disp:.1f} mm en R. "
                      f"La cabeza esta rotada en el CT y el IC sale sesgado.")

    print("")
    print("Transformaciones aplicadas:")
    print(f"{'PIEZA':<24}{'ET':<5}{'dR':>7}{'dA':>7}{'dS':>7}{'ROT':>8}  EJE")
    print("-" * 76)
    for nombre, info in sorted(_ESTADO.items()):
        if not info['activa'] or info['tNode'] is None:
            continue
        m = vtk.vtkMatrix4x4()
        info['tNode'].GetMatrixTransformToParent(m)
        d = _centroide_actual(info) - info['centroide0']
        ang, eje = _angulo_y_eje(m)
        if np.linalg.norm(d) < 0.05 and ang < 0.05:
            continue
        et = info['etiqueta'] or "-"
        print(f"{nombre:<24}{et:<5}{d[0]:>7.1f}{d[1]:>7.1f}{d[2]:>7.1f}{ang:>8.1f}  "
              f"({eje[0]:.2f}, {eje[1]:.2f}, {eje[2]:.2f})")

    desc = [(n, i) for n, i in _ESTADO.items() if not i['activa']]
    if desc:
        vol = sum(i['volumen_cm3'] for _, i in desc)
        print("")
        print(f"Hueso resecado: {vol:.2f} cm3 en {len(desc)} pieza(s) "
              f"({', '.join(n for n, _ in desc)}).")
    print("")


# ============================================================
# PARTE 7 - EXPORTACION
# ============================================================

def exportar_stl_G(carpeta_destino):
    """
    Exporta el armado FINAL: cada pieza activa con su transformacion aplicada
    (G1: el hardening ocurre aca, sobre una copia; los nodos de la escena
    quedan intactos y se pueden seguir moviendo).
    """
    if not _ESTADO:
        print("No hay estado. Corre preparar().")
        return
    if not os.path.isdir(carpeta_destino):
        os.makedirs(carpeta_destino)

    n = 0
    for nombre, info in sorted(_ESTADO.items()):
        if not info['activa']:
            continue
        poly = _poly_mundo(info['node'])
        if poly is None or poly.GetNumberOfPoints() == 0:
            print(f"  {nombre}: malla vacia, no se exporta.")
            continue
        etiqueta = f"_{info['etiqueta']}" if info['etiqueta'] else ""
        ruta = os.path.join(carpeta_destino, f"{nombre}{etiqueta}.stl")
        w = vtk.vtkSTLWriter()
        w.SetFileName(ruta)
        w.SetInputData(poly)
        w.SetFileTypeToBinary()
        w.Write()
        print(f"  {nombre} -> {ruta}")
        n += 1
    print(f"{n} STL exportado(s) con las transformaciones aplicadas.")


def limpiar_G():
    """Borra los transforms, el molde y los landmarks del Bloque G.
    Los fragmentos del Bloque F NO se tocan."""
    global _ESTADO
    for info in _ESTADO.values():
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
    _ESTADO = {}
    _ANCLA[0] = None
    print("Bloque G limpiado. Los fragmentos del Bloque F siguen intactos.")


# ============================================================

print(f"{CRANIOPLAN_G_VERSION} cargado.")
print(f"  IC objetivo por defecto: {IC_OBJETIVO}   carpeta: '{NOMBRE_CARPETA_SH}'")
print("Comandos:")
print("  preparar()                     -> inventario + transform por pieza + ancla")
print("  piezas()                       -> tabla de estado")
print("  landmarks()                    -> 3 puntos de linea media")
print("  molde(78)                      -> casquete objetivo por indice cefalico")
print("  manipular('Tapa_C1')           -> gizmo interactivo en esa pieza")
print("  mover('C', da=8)               -> traslacion en mm (dr/da/ds)")
print("  rotar('C', 'LR', 15)           -> rotacion en grados sobre eje anatomico")
print("  voltear('Resto_2', 'SI')       -> 180 grados (el flip real del cirujano)")
print("  ajustar('Tapa_C1')             -> ICP contra el molde")
print("  descartar('Fragmento_3')       -> saca la pieza (gris translucida)")
print("  restaurar('Fragmento_3')       -> la devuelve al armado")
print("  etiquetar('Resto_1', 'A')      -> vocabulario del quirofano")
print("  fijar('Hueso_1')               -> cambia el ancla")
print("  metricas()                     -> IC actual vs objetivo + transformaciones")
print("  resetear('C')  /  resetear()   -> vuelve al post-corte")
print("  exportar_stl_G('C:/salida')    -> STL con transformaciones aplicadas")
print("  limpiar_G()                    -> borra transforms, molde y landmarks")
