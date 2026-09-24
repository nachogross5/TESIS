# ============================================================
# CRANIOPLAN - BLOQUE G v1.1 (REACOMODAMIENTO DE PIEZAS POST-CORTE)
#
# Entra DESPUES del Bloque F. Toma los modelos de 'CranioPlan_Fragmentos',
# permite descartar piezas, moverlas y rotarlas, y mide el resultado.
#
# ############################################################
# CAMBIOS v1.1 RESPECTO DE v1.0 (los tres son correcciones de bugs
# encontrados probando sobre un caso real, no mejoras cosmeticas)
# ############################################################
#
# [F1] LA MOVILIDAD YA NO SE DECIDE POR EL NOMBRE.
#      v1.0 asumia la convencion del Bloque F (Hueso_N = isla preexistente,
#      fija; Tapa_/Resto_/Fragmento_ = nacida del corte, movible). En la
#      practica los modelos exportados llegaron llamandose 'Segment_5',
#      'Segment_17', etc., asi que TODO quedaba movible: la mandibula, las
#      vertebras y hasta una tubuladura se podian arrastrar.
#      Ahora la movilidad se decide GEOMETRICAMENTE: una pieza es movible si
#      su superficie pasa a menos de UMBRAL_CONTACTO_MM de alguna curva de
#      corte presente en la escena. Es el criterio literal que pedia el
#      cirujano ("lo que toco el corte se mueve") y no depende de ningun
#      nombre. Si no hay curvas en la escena, cae al criterio por nombre y
#      lo avisa.
#
# [F2] EL INDICE CEFALICO SE MIDE SOLO SOBRE LA BOVEDA.
#      v1.0 media el bounding box de TODAS las piezas activas. Como el ancla
#      trae la cara, la mandibula y el occipital por debajo del corte basal,
#      el bbox estaba dominado por piezas que no se mueven: al desplazar el
#      occipital 25mm el largo A-P bajaba 0.7mm, y un craneo escafocefalico
#      daba IC 77.7 (normocefalia). El numero era ruido.
#      Ahora se mide solo por ENCIMA de un nivel S de referencia (el plano de
#      la craneotomia basal, tomado de la curva cerrada de la escena) y solo
#      sobre piezas craneales, excluyendo las preexistentes (mandibula,
#      vertebras, tubuladuras). El nivel usado se imprime siempre.
#
# [F3] NUEVO: huecos(). Distancia minima de cada pieza a su vecina mas
#      cercana, en mm. Es el dato que el cirujano necesita para decidir que
#      separacion osifica sola, cual necesita injerto y donde va una placa.
#      El IC es un numero de diagnostico y de resultado; el hueco es una
#      decision intraoperatoria.
#
# ------------------------------------------------------------
# DECISIONES DE ARQUITECTURA (siguen vigentes desde v1.0)
# ------------------------------------------------------------
# [G1] LA GEOMETRIA NUNCA SE TOCA. Cada pieza movible recibe su propio
#      vtkMRMLLinearTransformNode. Todo el movimiento vive en esa matriz:
#      reset instantaneo, la matriz ES el dato del reporte prequirurgico, y
#      se recalculan metricas sin recalcular geometria. El hardening ocurre
#      solo al exportar STL, sobre copias.
#
# [G3] HAY UNA PIEZA ANCLA. Tocar el corte no alcanza para moverse: la base
#      craneal y el macizo facial tambien tocan la craneotomia basal y sin
#      embargo son el marco de referencia de toda la cirugia. Se detecta
#      como la pieza grande con el centroide mas inferior en S, y se corrige
#      con fijar().
#
#      Regla final de v1.1:
#         FIJA (gris claro):  piezas lejos de toda curva (mandibula,
#                             vertebras, tubuladuras)  +  el ANCLA
#         MOVIBLE (color):    piezas en contacto con alguna curva de corte,
#                             menos el ancla
#
# [G4] EL CENTRO DE TRANSFORMACION VA EN EL CENTROIDE. Sin esto el gizmo rota
#      alrededor del origen RAS y la pieza sale volando en el primer arrastre.
#
# [G5] EL MOLDE ES OPCIONAL Y DE BAJO VALOR. Un elipsoide liso derivado del
#      bbox del paciente. No corresponde a ninguna anatomia real y en la
#      practica tapa mas de lo que orienta. Queda como referencia visual
#      apagada por defecto: molde() hay que pedirlo explicitamente.
#
# [G6] NO EXISTE 'ESPEJAR' UNA PIEZA. Una reflexion tiene determinante
#      negativo y el hueso no se puede reflejar. Cuando el cirujano pasa una
#      placa de un lado al otro la VOLTEA 180 grados (rotacion propia, la
#      tabla interna queda mirando afuera). Por eso el comando es voltear().
#
# [G7] DESCARTAR ES UN ESTADO, NO UN BORRADO. La pieza descartada queda gris
#      translucida: se necesita para el reporte (volumen resecado) y porque
#      el cirujano cambia de idea.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS
# ------------------------------------------------------------
# - El IC se calcula sobre los ejes RAS del estudio. Si la cabeza esta
#   rotada en el CT, ancho y largo salen sesgados. metricas() avisa si los
#   landmarks de linea media muestran desviacion. Alineacion automatica:
#   pendiente.
# - Con transformaciones rigidas NO se puede ensanchar una placa. Ganar
#   ancho biparietal requiere barrel staves (doblado en tallo verde) o
#   partir el parietal en tiras. El modulo solo puede acortar el A-P y
#   trasladar placas enteras. Es una limitacion real, no un bug.
# - ajustar() (ICP) contra un elipsoide liso desliza tangencialmente: la
#   superficie no tiene rasgos que anclen la solucion. Pose inicial, no
#   resultado final.
# - huecos() mide distancia superficie-a-superficie por muestreo de puntos.
#   Es una cota, no una medida exacta de area de defecto.
# - No hay deteccion de colisiones: dos piezas pueden interpenetrarse y el
#   modulo no avisa. huecos() devuelve 0.0 en ese caso.
#
# FLUJO:
#   Bloque F cortar()  (dejar las curvas en la escena, NO borrarlas:
#                       v1.1 las necesita para decidir la movilidad)
#   -> preparar() -> piezas() -> etiquetar(...)
#   -> descartar(...) / manipular(...) / mover(...) / rotar(...) / voltear(...)
#   -> metricas() / huecos() -> exportar_stl_G('C:/salida')
# ============================================================

import os
import numpy as np
import vtk

CRANIOPLAN_G_VERSION = "Bloque G v1.1"

# ==================== CONFIGURACION ====================

NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"
NOMBRE_MOLDE           = "CranioPlan_Molde"
NOMBRE_LANDMARKS       = "CranioPlan_LineaMedia"

UMBRAL_CONTACTO_MM     = 5.0     # (F1) distancia pieza-curva para considerar
#   que la pieza nacio del corte. El kerf es ~1.2mm y la curva corre sobre la
#   tabla externa, asi que las piezas cortadas quedan a 1-3mm. La mandibula y
#   las vertebras estan a decenas de mm. 5.0 separa con holgura.
IC_OBJETIVO            = 78.0
PREFIJOS_FIJOS         = ("Hueso_",)   # solo se usa si NO hay curvas (fallback)
VOL_MIN_ANCLA_CM3      = 5.0
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

_ESTADO = {}
_ANCLA = [None]
_NIVEL_BOVEDA = [None]     # (F2) nivel S del corte basal


# ============================================================
# PARTE 1 - UTILIDADES
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


def _centroide_actual(info):
    m = vtk.vtkMatrix4x4()
    info['tNode'].GetMatrixTransformToWorld(m)
    c = info['centroide0']
    p = m.MultiplyPoint([c[0], c[1], c[2], 1.0])
    return np.array(p[:3], dtype=float)


def _componer(tNode, matrizMundo):
    actual = vtk.vtkMatrix4x4()
    tNode.GetMatrixTransformToParent(actual)
    t = vtk.vtkTransform()
    t.PostMultiply()
    t.SetMatrix(actual)
    t.Concatenate(matrizMundo)
    tNode.SetMatrixTransformToParent(t.GetMatrix())


def _buscar(nombre):
    if nombre in _ESTADO:
        return _ESTADO[nombre]
    for n, info in _ESTADO.items():
        if info['etiqueta'] and info['etiqueta'].upper() == str(nombre).upper():
            return info
    print(f"  no encuentro la pieza '{nombre}'. Corre piezas() para ver la lista.")
    return None


def _nombre_de(info):
    for n, i in _ESTADO.items():
        if i is info:
            return n
    return "?"


def _pintar(info):
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

def preparar(umbral_contacto_mm=UMBRAL_CONTACTO_MM):
    """
    Inventaria los fragmentos, decide movilidad POR CONTACTO CON LAS CURVAS
    (F1), crea un transform por pieza movible (G1) y detecta el ancla (G3).
    Se puede volver a correr: resetea el estado completo.

    IMPORTANTE: las curvas de corte tienen que seguir en la escena. Si las
    borraste, esto cae al criterio por nombre y lo avisa.
    """
    global _ESTADO
    shNode, itemId = _carpeta_fragmentos()
    if itemId is None:
        print(f"ERROR: no encuentro la carpeta '{NOMBRE_CARPETA_SH}'. Corre cortar().")
        return

    modelos = _modelos_de_carpeta(shNode, itemId)
    if not modelos:
        print("ERROR: la carpeta existe pero no tiene modelos adentro.")
        return

    for info in _ESTADO.values():
        t = info.get('tNode')
        if t is not None and slicer.mrmlScene.IsNodePresent(t):
            info['node'].SetAndObserveTransformNodeID(None)
            slicer.mrmlScene.RemoveNode(t)
    _ESTADO = {}
    _ANCLA[0] = None
    _NIVEL_BOVEDA[0] = None

    print(f"--- {CRANIOPLAN_G_VERSION} : preparar() ---")

    curvas = _curvas_de_corte()
    puntosCurva = _puntos_de_curvas(curvas)
    usaGeometria = puntosCurva is not None and len(puntosCurva) > 0

    if usaGeometria:
        print(f"Curvas de corte en escena: {[c.GetName() for c in curvas]}")
        print(f"Movilidad por contacto geometrico (umbral {umbral_contacto_mm:.1f} mm).")
        # (F2) nivel de boveda = S medio de la curva CERRADA (craneotomia basal)
        cerradas = [c for c in curvas if c.IsA('vtkMRMLMarkupsClosedCurveNode')]
        if cerradas:
            p = _puntos_de_curvas([cerradas[0]])
            if p is not None and len(p):
                _NIVEL_BOVEDA[0] = float(np.mean(p[:, 2]))
                print(f"Nivel de boveda tomado de '{cerradas[0].GetName()}': "
                      f"S = {_NIVEL_BOVEDA[0]:.1f} mm")
    else:
        print("AVISO: no hay curvas de corte en la escena.")
        print("  Caigo al criterio por NOMBRE, que falla si las piezas se llaman")
        print("  'Segment_N'. Volve a correr el corte sin borrar las curvas, o")
        print("  marca las fijas a mano con fijar_pieza('Nombre').")

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

        _ESTADO[nombre] = {
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

    candidatas = [(n, i) for n, i in _ESTADO.items()
                  if i['movible'] and i['volumen_cm3'] >= VOL_MIN_ANCLA_CM3]
    if not candidatas:
        candidatas = [(n, i) for n, i in _ESTADO.items() if i['movible']]
    if candidatas:
        nombreAncla = min(candidatas, key=lambda kv: kv[1]['centroide0'][2])[0]
        _fijar_interno(nombreAncla, avisar=False)
        print(f"Ancla detectada: '{nombreAncla}' (centroide mas inferior en S).")
        print("Si no es la base craneal + macizo facial, corregilo con fijar('Nombre').")

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
        ant['movible'] = not ant['preexistente']
        _pintar(ant)

    info['movible'] = False
    if info['tNode'] is not None:
        info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
        dn = info['tNode'].GetDisplayNode()
        if dn is not None:
            dn.SetEditorVisibility(False)
    _ANCLA[0] = _nombre_de(info)
    _pintar(info)
    if avisar:
        print(f"Ancla: '{_ANCLA[0]}'. Queda fija y gris; el resto se mueve respecto a ella.")


def fijar(nombre):
    """Define manualmente la pieza ancla (base craneal + macizo facial)."""
    _fijar_interno(nombre, avisar=True)


def fijar_pieza(nombre):
    """Marca una pieza como FIJA preexistente (mandibula, vertebras, tubuladura).
    Util si el criterio automatico la dejo movible."""
    info = _buscar(nombre)
    if info is None:
        return
    info['movible'] = False
    info['preexistente'] = True
    if info['tNode'] is not None:
        info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
        info['node'].SetAndObserveTransformNodeID(None)
        slicer.mrmlScene.RemoveNode(info['tNode'])
        info['tNode'] = None
    _pintar(info)
    print(f"'{nombre}' marcada como fija preexistente (excluida de las metricas).")


def liberar_pieza(nombre):
    """Vuelve movible una pieza marcada como fija por error."""
    info = _buscar(nombre)
    if info is None:
        return
    if info['tNode'] is None:
        t = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLLinearTransformNode', f"T_{_nombre_de(info)}")
        c = info['centroide0']
        t.SetCenterOfTransformation(c[0], c[1], c[2])
        info['node'].SetAndObserveTransformNodeID(t.GetID())
        info['tNode'] = t
    info['movible'] = True
    info['preexistente'] = False
    _pintar(info)
    print(f"'{nombre}' liberada: ahora es movible.")


def piezas():
    """Tabla de estado de todas las piezas."""
    if not _ESTADO:
        print("No hay estado. Corre preparar().")
        return
    print("")
    print(f"{'PIEZA':<22}{'ET':<4}{'VOL cm3':>9}{'d.CURVA':>9}  {'ESTADO':<16}{'DESPL mm':>10}")
    print("-" * 74)
    for nombre, info in sorted(_ESTADO.items()):
        if not info['activa']:
            estado = "DESCARTADA"
        elif nombre == _ANCLA[0]:
            estado = "ANCLA (fija)"
        elif info['preexistente']:
            estado = "fija (preexist)"
        elif not info['movible']:
            estado = "fija"
        else:
            estado = "movible"
        d = 0.0
        if info['tNode'] is not None:
            d = float(np.linalg.norm(_centroide_actual(info) - info['centroide0']))
        dc = info['dist_curva']
        dcs = "  n/a" if (dc != dc or dc == float('inf')) else f"{dc:>9.1f}"
        et = info['etiqueta'] or "-"
        print(f"{nombre:<22}{et:<4}{info['volumen_cm3']:>9.2f}{dcs}  {estado:<16}{d:>10.1f}")
    activas = [i for i in _ESTADO.values() if i['activa']]
    desc = [i for i in _ESTADO.values() if not i['activa']]
    volDesc = sum(i['volumen_cm3'] for i in desc)
    print("-" * 74)
    print(f"{len(activas)} pieza(s) en el armado, {len(desc)} descartada(s) "
          f"({volDesc:.2f} cm3 de hueso resecado).")
    print("La columna d.CURVA es la distancia al corte: chica = nacio del corte.")
    print("")


def etiquetar(nombre, letra):
    """Renombra con el vocabulario del quirofano (A, B, C, R, L, F)."""
    info = _buscar(nombre)
    if info is None:
        return
    info['etiqueta'] = str(letra).upper()
    print(f"'{_nombre_de(info)}' etiquetada como '{info['etiqueta']}'.")


# ============================================================
# PARTE 3 - DESCARTE
# ============================================================

def descartar(nombre):
    """(G7) Saca la pieza del armado sin borrarla: gris translucida."""
    info = _buscar(nombre)
    if info is None:
        return
    if _nombre_de(info) == _ANCLA[0]:
        print("No podes descartar el ancla. Primero asigna otra con fijar().")
        return
    info['activa'] = False
    if info['tNode'] is not None:
        dn = info['tNode'].GetDisplayNode()
        if dn is not None:
            dn.SetEditorVisibility(False)
    _pintar(info)
    print(f"'{_nombre_de(info)}' descartada ({info['volumen_cm3']:.2f} cm3).")


def restaurar(nombre):
    """Devuelve al armado una pieza descartada."""
    info = _buscar(nombre)
    if info is None:
        return
    info['activa'] = True
    _pintar(info)
    print(f"'{_nombre_de(info)}' restaurada al armado.")


# ============================================================
# PARTE 4 - MANIPULACION
# ============================================================

def manipular(nombre, encender=True):
    """Enciende el gizmo interactivo. Rota alrededor del centroide (G4)."""
    info = _buscar(nombre)
    if info is None:
        return
    if info['tNode'] is None or not info['movible']:
        print(f"'{_nombre_de(info)}' es una pieza FIJA. Si tiene que moverse, "
              f"usa liberar_pieza() o reasigna el ancla con fijar().")
        return
    if not info['activa']:
        print(f"'{_nombre_de(info)}' esta descartada. Corre restaurar() primero.")
        return

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
        for i in _ESTADO.values():
            if i is not info and i['tNode'] is not None:
                d = i['tNode'].GetDisplayNode()
                if d is not None:
                    d.SetEditorVisibility(False)
        print(f"Gizmo activo en '{_nombre_de(info)}'. Arrastra flechas (traslacion) "
              f"o aros (rotacion) en la vista 3D.")
    else:
        print(f"Gizmo apagado en '{_nombre_de(info)}'.")


def mover(nombre, dr=0.0, da=0.0, ds=0.0):
    """
    Traslacion en mm sobre ejes anatomicos RAS.
      dr > 0 -> derecha del paciente | da > 0 -> anterior | ds > 0 -> superior
    """
    info = _buscar(nombre)
    if info is None or info['tNode'] is None:
        return
    t = vtk.vtkTransform()
    t.Translate(float(dr), float(da), float(ds))
    _componer(info['tNode'], t.GetMatrix())
    d = float(np.linalg.norm(_centroide_actual(info) - info['centroide0']))
    print(f"'{_nombre_de(info)}' movida ({dr:+.1f}, {da:+.1f}, {ds:+.1f}) mm. "
          f"Desplazamiento total: {d:.1f} mm.")


def rotar(nombre, eje, grados):
    """Rotacion en grados alrededor del CENTROIDE. eje: 'LR', 'AP' o 'SI'."""
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
    print(f"'{_nombre_de(info)}' rotada {grados:+.1f} grados sobre el eje {key}.")


def voltear(nombre, eje="SI"):
    """
    (G6) Voltea 180 grados. NO es un espejo: una reflexion tiene determinante
    negativo y el hueso no se puede reflejar. Es la rotacion propia que hace
    el cirujano al pasar una placa de un lado al otro; la tabla interna queda
    mirando hacia afuera, igual que en el quirofano.
      'SI' -> occipital pasa a frontal (swap clasico de escafocefalia)
      'AP' -> cruza izquierda-derecha invirtiendo arriba-abajo
    """
    rotar(nombre, eje, 180.0)


def resetear(nombre=None):
    """Devuelve una pieza (o todas) a su posicion original post-corte."""
    if nombre is None:
        objetivos = list(_ESTADO.values())
    else:
        info = _buscar(nombre)
        objetivos = [info] if info else []
    n = 0
    for info in objetivos:
        if info['tNode'] is not None:
            info['tNode'].SetMatrixTransformToParent(vtk.vtkMatrix4x4())
            n += 1
    print(f"{n} pieza(s) devuelta(s) a su posicion original.")


# ============================================================
# PARTE 5 - MEDICION
# ============================================================

def _piezas_de_medicion():
    """(F2) Piezas que cuentan para el IC: activas y craneales.
    Se excluyen las preexistentes (mandibula, vertebras, tubuladuras)."""
    return [i for i in _ESTADO.values() if i['activa'] and not i['preexistente']]


def _bbox_boveda(nivel_S):
    """Bounding box de la boveda: solo la porcion por ENCIMA de nivel_S."""
    bounds = None
    for info in _piezas_de_medicion():
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


def metricas(nivel_S=None):
    """
    (F2) Indice cefalico medido SOLO sobre la boveda, por encima del plano de
    la craneotomia basal, y solo sobre piezas craneales.

    nivel_S : plano de referencia en RAS. None = el detectado en preparar()
              a partir de la curva cerrada. Pasalo a mano para forzarlo.
    """
    if not _ESTADO:
        print("No hay estado. Corre preparar().")
        return

    if nivel_S is None:
        nivel_S = _NIVEL_BOVEDA[0]

    b = _bbox_boveda(nivel_S)
    if b is None:
        print("No hay piezas craneales activas por encima del nivel indicado.")
        return
    L = b[3] - b[2]
    W = b[1] - b[0]
    H = b[5] - b[4]
    ic = 100.0 * W / L if L > 1e-6 else 0.0

    print("")
    print(f"--- METRICAS ({CRANIOPLAN_G_VERSION}) ---")
    if nivel_S is None:
        print("AVISO: sin nivel de boveda. Estoy midiendo TODO, incluida la cara")
        print("  y la base: el IC va a salir falsamente bajo. Pasa nivel_S a mano.")
    else:
        print(f"Medido por encima de S = {nivel_S:.1f} mm (plano de craneotomia).")
    excl = [_nombre_de(i) for i in _ESTADO.values() if i['preexistente']]
    if excl:
        print(f"Excluidas de la medicion: {excl}")

    print(f"Largo A-P : {L:.1f} mm")
    print(f"Ancho L-R : {W:.1f} mm")
    print(f"Alto  S-I : {H:.1f} mm")
    print(f"Indice cefalico : {ic:.1f}   (objetivo {IC_OBJETIVO:.1f})")
    if ic < IC_OBJETIVO - 1.0:
        falta_L = L - (W * 100.0 / IC_OBJETIVO)
        print(f"  faltan {IC_OBJETIVO - ic:.1f} puntos de IC.")
        print(f"  a ancho constante habria que acortar {falta_L:.1f} mm en A-P,")
        print(f"  o ensanchar {(L * IC_OBJETIVO / 100.0) - W:.1f} mm en transverso")
        print("  (el ancho NO se gana con transformaciones rigidas: requiere")
        print("   barrel staves o partir la placa en tiras).")
    elif ic > IC_OBJETIVO + 1.0:
        print(f"  te pasaste {ic - IC_OBJETIVO:.1f} puntos de IC.")
    else:
        print("  IC dentro del objetivo.")

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
                print(f"  AVISO: la linea media se desvia {disp:.1f} mm en R. "
                      f"La cabeza esta rotada en el CT y el IC sale sesgado.")

    print("")
    print("Transformaciones aplicadas:")
    print(f"{'PIEZA':<22}{'ET':<4}{'dR':>7}{'dA':>7}{'dS':>7}{'ROT':>8}  EJE")
    print("-" * 74)
    hubo = False
    for nombre, info in sorted(_ESTADO.items()):
        if not info['activa'] or info['tNode'] is None:
            continue
        m = vtk.vtkMatrix4x4()
        info['tNode'].GetMatrixTransformToParent(m)
        d = _centroide_actual(info) - info['centroide0']
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
        print(f"{nombre:<22}{et:<4}{d[0]:>7.1f}{d[1]:>7.1f}{d[2]:>7.1f}{ang:>8.1f}  "
              f"({eje[0]:.2f}, {eje[1]:.2f}, {eje[2]:.2f})")
    if not hubo:
        print("  (ninguna pieza movida todavia)")

    desc = [(n, i) for n, i in _ESTADO.items() if not i['activa']]
    if desc:
        vol = sum(i['volumen_cm3'] for _, i in desc)
        print("")
        print(f"Hueso resecado: {vol:.2f} cm3 en {len(desc)} pieza(s) "
              f"({', '.join(n for n, _ in desc)}).")
    print("")


def huecos(muestreo=MUESTREO_HUECOS):
    """
    (F3) Distancia minima de cada pieza activa a su vecina mas cercana.

    Es el dato que decide conducta en quirofano: que separacion osifica sola,
    cual necesita injerto, donde va una placa. Un 0.0 significa contacto o
    interpenetracion (el modulo no distingue: no hay deteccion de colisiones).

    LIMITE: se mide por muestreo de puntos de superficie, asi que es una cota
    razonable pero no una medida exacta de area de defecto.
    """
    activas = [(n, i) for n, i in _ESTADO.items() if i['activa']]
    if len(activas) < 2:
        print("Hacen falta al menos dos piezas activas.")
        return

    polys = {}
    for n, i in activas:
        p = _poly_mundo(i['node'])
        if p is not None and p.GetNumberOfPoints() > 0:
            polys[n] = p

    locs = {}
    for n, p in polys.items():
        l = vtk.vtkPointLocator()
        l.SetDataSet(p)
        l.BuildLocator()
        locs[n] = l

    print("")
    print(f"--- HUECOS ENTRE PIEZAS ({CRANIOPLAN_G_VERSION}) ---")
    print(f"{'PIEZA':<22}{'VECINA MAS CERCANA':<24}{'HUECO mm':>10}")
    print("-" * 58)

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
        et = _ESTADO[n]['etiqueta']
        etq = f" ({et})" if et else ""
        print(f"{n + etq:<22}{mejor or '-':<24}{mejorD:>10.1f}")

    print("-" * 58)
    print("0.0 = contacto o interpenetracion (no hay deteccion de colisiones).")
    print("")


# ============================================================
# PARTE 6 - MOLDE (opcional, ver G5)
# ============================================================

def molde(ic_objetivo=IC_OBJETIVO, nivel_S=None):
    """
    (G5) Casquete elipsoidal de referencia, semitransparente. OPCIONAL y de
    valor limitado: es un elipsoide liso derivado del bbox del paciente, no
    corresponde a ninguna anatomia real y tiende a tapar mas que orientar.
    Preserva el producto L*W (el hueso no se estira) y reparte con el IC
    pedido:  L = sqrt(L0*W0*100/IC) ,  W = L*IC/100
    """
    if nivel_S is None:
        nivel_S = _NIVEL_BOVEDA[0]
    b = _bbox_boveda(nivel_S)
    if b is None:
        print("No hay boveda medible. Corre preparar().")
        return

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

    print(f"Molde de referencia (IC {ic_objetivo:.1f}).")
    print(f"  boveda actual : largo {L0:.1f}  ancho {W0:.1f}  -> IC {100.0*W0/L0:.1f}")
    print(f"  objetivo      : largo {L:.1f}  ancho {W:.1f}")
    return nodo


def ocultar_molde():
    nodo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MOLDE)
    if nodo is not None and nodo.GetDisplayNode():
        nodo.GetDisplayNode().SetVisibility(False)
        print("Molde oculto.")


def landmarks():
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
        print(f"Nodo '{NOMBRE_LANDMARKS}' creado con 3 puntos sin colocar.")
    print("Colocalos sobre la linea media y corre metricas().")
    return nodo


def ajustar(nombre, iteraciones=60):
    """ICP rigido contra el molde. Pose inicial, no resultado final: el
    elipsoide es liso y el ICP desliza tangencialmente."""
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
    icp.GetLandmarkTransform().SetModeToRigidBody()
    icp.SetMaximumNumberOfIterations(int(iteraciones))
    icp.SetCheckMeanDistance(1)
    icp.SetMaximumMeanDistance(0.01)
    icp.StartByMatchingCentroidsOff()
    icp.Modified()
    icp.Update()
    m = vtk.vtkMatrix4x4()
    m.DeepCopy(icp.GetMatrix())
    _componer(info['tNode'], m)
    print(f"'{_nombre_de(info)}' ajustada al molde por ICP. Revisa a ojo.")


# ============================================================
# PARTE 7 - EXPORTACION
# ============================================================

def exportar_stl_G(carpeta_destino):
    """Exporta el armado final: cada pieza activa con su transformacion
    aplicada (G1: el hardening ocurre aca, sobre copias)."""
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
    """Borra transforms, molde y landmarks. Los fragmentos NO se tocan."""
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
    _NIVEL_BOVEDA[0] = None
    print("Bloque G limpiado. Los fragmentos del Bloque F siguen intactos.")


# ============================================================

print(f"{CRANIOPLAN_G_VERSION} cargado.")
print("  v1.1: movilidad por contacto con las curvas (no por nombre),")
print("        IC medido solo sobre la boveda, y huecos() en mm.")
print("  IMPORTANTE: no borres las curvas de corte antes de preparar().")
print("Comandos:")
print("  preparar()                     -> inventario + movilidad + ancla")
print("  piezas()                       -> tabla de estado (mira d.CURVA)")
print("  etiquetar('Segment_17','C')    -> vocabulario del quirofano")
print("  fijar('Segment_3')             -> cambia el ancla")
print("  fijar_pieza('Segment_9')       -> marca fija una preexistente")
print("  liberar_pieza('Segment_9')     -> la vuelve movible")
print("  manipular('C')                 -> gizmo interactivo")
print("  mover('C', da=8)               -> traslacion en mm (dr/da/ds)")
print("  rotar('C','LR',15)             -> rotacion en grados")
print("  voltear('C','SI')              -> 180 grados (el flip del cirujano)")
print("  descartar('B') / restaurar('B')")
print("  metricas()                     -> IC de la boveda + transformaciones")
print("  huecos()                       -> separacion entre piezas en mm")
print("  molde(78) / ocultar_molde()    -> referencia opcional")
print("  resetear('C') / resetear()     -> vuelve al post-corte")
print("  exportar_stl_G('C:/salida')")
print("  limpiar_G()")
