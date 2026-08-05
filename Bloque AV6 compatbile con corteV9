# ============================================================
# CRANIOPLAN - BLOQUE A v6  (script de consola para 3D Slicer)
#            COMPATIBLE CON EL CORTE V9
# ============================================================
# Reemplaza al Bloque A V5. Se pega entero en la consola de Python de
# Slicer. Al pegarlo se ejecuta solo y deja los comandos disponibles.
#
# ------------------------------------------------------------
# QUE CAMBIA RESPECTO DE LA VERSION ANTERIOR DE ESTE SCRIPT
# ------------------------------------------------------------
# SOLO los nombres de los nodos que quedan en la escena. El algoritmo es
# identico. El motivo: el Corte V9 busca los nodos por nombre fijo
#
#     NOMBRE_SEGMENTACION   = "Craneo_Automatico"   <- nodo de segmentacion
#     NOMBRE_SEGMENTO_HUESO = "Craneo_Final"        <- segmento adentro
#     NOMBRE_MODELO_CRANEO  = "Craneo_Final"        <- Model node
#
# y la version anterior del Bloque A v6 dejaba "CranioPlan_Segmentacion" y
# "Craneo_Final_Modelo". Por eso cortar() cortaba con
#     ERROR: no encuentro el nodo de segmentacion 'Craneo_Automatico'.
#
# Ahora los tres nombres estan en el bloque CONFIG de abajo. Si alguna vez
# se cambia la config del Corte V9, se cambian aca y listo.
#
# OJO: el segmento se llama "Craneo_Final" y el Model node TAMBIEN se llama
# "Craneo_Final". No chocan: un segmento no es un nodo de la escena, asi que
# slicer.util.getNode("Craneo_Final") devuelve siempre el modelo. Es
# confuso de leer pero es lo que V9 espera; el modulo CranioPlan usa los
# nombres largos y sin ambiguedad.
#
# ------------------------------------------------------------
# PROBLEMAS DE V5 QUE ESTE SCRIPT CORRIGE
# ------------------------------------------------------------
#
# [P1] TC 175 - EL CRANEO DESAPARECIA.
#      Sintoma: "Isla de referencia (mayor volumen, sin tocar borde):
#      vol=2.84cm3" y en el visor solo aparecian vertebras y dos restos.
#      CAUSA: la calota fue una de las islas descartadas por tocar el borde
#      del volumen. En esa serie el FOV recorta el craneo, asi que el veto
#      por borde lo elimino y como referencia quedo una esquirla de 2.84 cm3.
#      FIX: seleccion en DOS PASADAS (ver [P7]).
#
# [P2] TC 559 - MemoryError, nunca generaba el craneo.
#      CAUSA: scipy.ndimage.distance_transform_edt SIEMPRE reserva un feature
#      transform de forma (3, Z, Y, X) en int32 -> 1.84 GB para ese volumen,
#      y V5 lo llamaba UNA VEZ POR RONDA. Ademas guardaba una mascara de
#      volumen completo por isla: 31 islas x 153 MB = 4.7 GB.
#      FIX: cKDTree sobre los voxeles de SUPERFICIE de cada isla, matriz de
#      distancias calculada UNA sola vez, y un unico array de etiquetas.
#
# [P3] TC 403 - "Craneo_Final" se comia la mandibula.
#      CAUSA: la limpieza de malla usaba un umbral RELATIVO al 10% del
#      tamano de la region mayor y conservaba UNA sola region de 454.
#      FIX: umbral ABSOLUTO en cantidad de puntos (MINIMO_PUNTOS_MALLA).
#
# [P4] TC 403 - la almohadilla no se eliminaba entera.
#      CAUSA: (a) la espuma se parte en varios componentes al umbralar;
#      (b) una TC reconstruye un CIRCULO inscripto y rellena el resto con un
#      valor fijo, asi que lo cortado por ese borde nunca llega a una cara
#      plana del array y el test de las 6 caras no lo ve.
#      FIX: filtro de densidad por percentil 95 de HU + deteccion del valor
#      de relleno en las esquinas + detector de losa (tabla de camilla).
#
# [P5] Confusion Craneo_Automatico / Craneo_Final.
#      FIX: nombres explicitos en CONFIG y un unico resultado final.
#
# [P6] Filtro de tamano relativo (5% de la isla mayor): en tu propia corrida
#      del 403 aparecian islas de 7.44, 5.28 y 5.09 cm3, justo en el filo del
#      corte, y una mandibula entera ronda ese orden.
#      FIX: piso ABSOLUTO en cm3.
#
# [P7] TC de JUAN ETCHEVARNE - la almohadilla quedaba como referencia
#      (79.49 cm3, p95=1609 HU, tocando el borde) por delante del craneo
#      (58.74 cm3, p95=1429). En ESE equipo el soporte es tan denso como
#      hueso cortical: la densidad NO es un discriminador universal.
#      Los dos casos son simetricos:
#           Juan   -> la basura toca el borde y el craneo no
#           TC 175 -> el craneo toca el borde y lo unico que no toca es una
#                     esquirla de 2.84 cm3
#      FIX: dos pasadas. Pasada 1 con veto de borde para TODAS (logica de V5,
#      que resuelve Juan). Si lo aceptado no llega a VOLUMEN_CRANEO_MINIMO_CM3,
#      pasada 2 de rescate con la referencia exenta (resuelve la 175).
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (documentados a proposito, no ignorados)
# ------------------------------------------------------------
# - VERTEBRAS: son hueso real. Ningun criterio de material las distingue de
#   la mandibula. Se eliminan solas SOLO si el FOV las recorta. Si entran
#   completas hay que borrarlas a mano con eliminar(n).
# - ALMOHADILLA CON INSERTO METALICO: pasaria densidad y losa.
# - LA DENSIDAD NO ES UNIVERSAL. Confirmado en la TC de Juan: p95=1609 HU en
#   el soporte de inmovilizacion. Quien resuelve esos casos es el veto de
#   borde. Los dos filtros son complementarios.
# - UMBRAL_DENSIDAD_HU = 500 esta elegido a priori. Salvaguarda: si deja el
#   pool vacio se desactiva solo y avisa.
# ============================================================


# ============================================================
# CONFIG
# ============================================================

# --- Nombres de nodos (deben coincidir con la config del Corte V9) ---
NOMBRE_NODO_SEGMENTACION = "Craneo_Automatico"   # V9: NOMBRE_SEGMENTACION
NOMBRE_SEGMENTO_CRANEO   = "Craneo_Final"        # V9: NOMBRE_SEGMENTO_HUESO
NOMBRE_MODELO_CRANEO     = "Craneo_Final"        # V9: NOMBRE_MODELO_CRANEO

# Volumen a usar. None = el ultimo volumen escalar cargado en la escena.
# Se puede poner el nombre exacto del nodo, o un fragmento del nombre.
NOMBRE_VOLUMEN = None

# Umbrales de hueso (Hounsfield).
HU_MIN = 300
HU_MAX = 3000

# Piso ABSOLUTO de volumen para que una isla llegue a la revision (ver [P6]).
VOL_MINIMO_CM3 = 0.5

# Distancia maxima a la masa aceptada para incorporar una isla.
# Cubre huesos separados por suturas abiertas.
MARGEN_PROXIMIDAD_MM = 40.0

# Volumen oseo TOTAL minimo para considerar que lo aceptado es un craneo.
# Cota anatomica, no parametro de ajuste: solo decide si hace falta el
# rescate de la pasada 2 (ver [P7]).
VOLUMEN_CRANEO_MINIMO_CM3 = 20.0

# Filtro de densidad (ver [P4]). Percentil 95 de HU dentro de la isla.
UMBRAL_DENSIDAD_HU = 500.0
UMBRAL_DENSIDAD_DUDOSA_HU = 700.0

# Detector de losa de soporte (tabla de la camilla).
ESPESOR_LOSA_MM = 30.0
FRACCION_LOSA_FOV = 0.85

# Exportacion de malla.
MINIMO_PUNTOS_MALLA = 200      # umbral ABSOLUTO de ruido (ver [P3])
REDUCCION_MALLA = 0.7          # 0.0 desactiva la decimacion
SUAVIZADO_SUPERFICIE = "0.3"   # bajar a "0.1" si se disuelven placas finas

# Rendimiento.
MAX_PUNTOS_KDTREE = 60000      # puntos de superficie por isla para distancias


# ============================================================
# CODIGO
# ============================================================

import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

CRANIOPLAN_BLOQUE_A = "v6 (2026-08-05) consola - compatible Corte V9"

_estado = {
    "volumeNode": None,
    "labels": None,
    "espaciadoZYX": None,
    "filas": [],          # todas las islas, con metricas y decision
    "candidatas": [],     # dicts de las que llegan a la revision
    "segNode": None,
    "colores": {},
    "modelo": None,
}


# ------------------------------------------------------------
# Utilidades puras (numpy) - sin dependencias de Slicer
# ------------------------------------------------------------

def _valorDeRelleno(volArr):
    """
    Valor con el que el tomografo rellena lo que queda FUERA del circulo de
    reconstruccion. Se lee en las cuatro esquinas de un corte central.
    """
    z = volArr.shape[0] // 2
    corte = volArr[z]
    esquinas = {int(corte[0, 0]), int(corte[0, -1]),
                int(corte[-1, 0]), int(corte[-1, -1])}
    if len(esquinas) == 1:
        valor = esquinas.pop()
        if valor < -900:
            return valor
    return None


def _mascaraExterior(volArr, valorRelleno):
    """
    Mascara booleana del "afuera": el relleno fuera del circulo de
    reconstruccion mas las seis caras planas del array. Cubre de una sola vez
    el borde CIRCULAR del FOV y el borde del array.
    """
    exterior = np.zeros(volArr.shape, dtype=bool)
    if valorRelleno is not None:
        exterior |= (volArr == valorRelleno)
    exterior[0, :, :] = True
    exterior[-1, :, :] = True
    exterior[:, 0, :] = True
    exterior[:, -1, :] = True
    exterior[:, :, 0] = True
    exterior[:, :, -1] = True
    return exterior


def _etiquetasTocandoExterior(labels, exterior):
    """
    Etiquetas con al menos un voxel dentro del exterior o pegado a el. Se
    resuelve con seis desplazamientos en vez de dilatar el volumen completo:
    la dilatacion pediria otro array de volumen entero, que es lo que hacia
    fallar a V5 por memoria.
    """
    tocan = set()
    tocan.update(np.unique(labels[exterior]).tolist())
    tocan.update(np.unique(labels[:-1, :, :][exterior[1:, :, :]]).tolist())
    tocan.update(np.unique(labels[1:, :, :][exterior[:-1, :, :]]).tolist())
    tocan.update(np.unique(labels[:, :-1, :][exterior[:, 1:, :]]).tolist())
    tocan.update(np.unique(labels[:, 1:, :][exterior[:, :-1, :]]).tolist())
    tocan.update(np.unique(labels[:, :, :-1][exterior[:, :, 1:]]).tolist())
    tocan.update(np.unique(labels[:, :, 1:][exterior[:, :, :-1]]).tolist())
    tocan.discard(0)
    return tocan


def _esLosaDeSoporte(extMM, fracFOV, espesorMaximoMM, fraccionMinimaFOV):
    """
    Tabla de la camilla: abarca casi todo el campo en algun eje y es delgada
    en su eje mas fino. La caja envolvente del craneo no cumple eso (su eje
    mas fino ronda los 90-100 mm). Una almohadilla CONCAVA tampoco lo cumple
    -a esa la agarra el filtro de densidad-, por eso los dos filtros son
    necesarios y ninguno reemplaza al otro.
    """
    if max(extMM) <= 0.0:
        return False
    return bool(min(extMM) <= espesorMaximoMM and max(fracFOV) >= fraccionMinimaFOV)


def _puntosDeSuperficie(labels, etiqueta, caja, espaciadoZYX, maxPuntos):
    """
    Coordenadas fisicas (mm) de los voxeles de SUPERFICIE de una isla.

    La distancia minima entre dos solidos se alcanza siempre en sus
    superficies: da el mismo resultado que la EDT de volumen completo con dos
    ordenes de magnitud menos de puntos.
    """
    sub = (labels[caja] == etiqueta)
    interior = ndimage.binary_erosion(sub)
    borde = sub & ~interior
    if not borde.any():
        borde = sub
    idx = np.argwhere(borde)
    offset = np.array([caja[0].start, caja[1].start, caja[2].start])
    idx = idx + offset
    if len(idx) > maxPuntos:
        paso = int(np.ceil(len(idx) / float(maxPuntos)))
        idx = idx[::paso]
    return idx.astype(np.float64) * np.asarray(espaciadoZYX, dtype=np.float64)


def analizarVolumen(volArr, espaciadoZYX, cfg):
    """
    Etiqueta las islas y calcula SOLO metricas. No decide nada: la decision
    la toma seleccionarCandidatas(), que puede correrse dos veces con
    criterios distintos sin volver a etiquetar el volumen (que es lo caro).
    """
    mascara = (volArr >= cfg["HU_MIN"]) & (volArr <= cfg["HU_MAX"])
    estructura = ndimage.generate_binary_structure(3, 3)  # 26-conexo
    labels, nIslas = ndimage.label(mascara, structure=estructura)
    del mascara
    if nIslas == 0:
        return None, [], {"nIslas": 0, "valorRelleno": None}
    if nIslas < 32000:
        labels = labels.astype(np.int16)

    conteos = np.bincount(labels.ravel())
    cajas = ndimage.find_objects(labels)

    valorRelleno = _valorDeRelleno(volArr)
    exterior = _mascaraExterior(volArr, valorRelleno)
    etiquetasBorde = _etiquetasTocandoExterior(labels, exterior)
    del exterior

    volVoxCM3 = (espaciadoZYX[0] * espaciadoZYX[1] * espaciadoZYX[2]) / 1000.0

    filas = []
    for k in range(1, nIslas + 1):
        caja = cajas[k - 1]
        if caja is None:
            continue
        vol = int(conteos[k]) * volVoxCM3
        extMM = tuple((caja[e].stop - caja[e].start) * espaciadoZYX[e] for e in range(3))
        fracFOV = tuple((caja[e].stop - caja[e].start) / float(volArr.shape[e])
                        for e in range(3))

        # El percentil solo se calcula si la isla supera el piso de volumen:
        # en la TC de Juan esto evita 7919 percentiles inutiles.
        p95 = None
        if vol >= cfg["VOL_MINIMO_CM3"]:
            sub = (labels[caja] == k)
            valores = volArr[caja][sub]
            if valores.size:
                p95 = float(np.percentile(valores, 95))

        filas.append({
            "etiqueta": k, "vol": vol, "caja": caja,
            "extMM": extMM, "fracFOV": fracFOV, "p95": p95,
            "tocaBorde": k in etiquetasBorde,
            "esLosa": _esLosaDeSoporte(extMM, fracFOV,
                                       cfg["ESPESOR_LOSA_MM"],
                                       cfg["FRACCION_LOSA_FOV"]),
            "motivo": None, "dist": 0.0,
            "referencia": False, "protegida": False,
        })

    return labels, filas, {"nIslas": nIslas, "valorRelleno": valorRelleno}


def _clasificar(filas, cfg, aplicarDensidad=True):
    """
    Asigna el motivo de descarte de cada isla. El veto de borde se aplica
    SIEMPRE aca, a todas por igual; la exencion de la referencia -cuando
    corresponde- se resuelve despues, en seleccionarCandidatas().
    """
    for f in filas:
        f["motivo"] = None
        f["dist"] = 0.0
        f["referencia"] = False
        f["protegida"] = False

    for f in filas:
        if f["vol"] < cfg["VOL_MINIMO_CM3"]:
            f["motivo"] = "volumen chico (%.2f < %.2f cm3)" % (f["vol"], cfg["VOL_MINIMO_CM3"])
        elif f["esLosa"]:
            f["motivo"] = "losa de soporte (tabla de camilla)"
        elif (aplicarDensidad and cfg["UMBRAL_DENSIDAD_HU"] > 0
              and f["p95"] is not None and f["p95"] < cfg["UMBRAL_DENSIDAD_HU"]):
            f["motivo"] = "densidad no osea (p95=%.0f < %.0f HU)" % (
                f["p95"], cfg["UMBRAL_DENSIDAD_HU"])
        elif f["tocaBorde"]:
            f["motivo"] = "toca el borde del FOV o del array"


def _esHuesoPlausible(f, cfg, aplicarDensidad=True):
    """Candidata a REFERENCIA por sus propiedades de material y tamano.
    A proposito no mira el borde: de eso se ocupa quien la llama."""
    if f["vol"] < cfg["VOL_MINIMO_CM3"]:
        return False
    if f["esLosa"]:
        return False
    if (aplicarDensidad and cfg["UMBRAL_DENSIDAD_HU"] > 0
            and f["p95"] is not None and f["p95"] < cfg["UMBRAL_DENSIDAD_HU"]):
        return False
    return True


def seleccionarCandidatas(labels, filas, espaciadoZYX, cfg,
                          exentarReferenciaDelBorde, cache):
    """
    Una pasada completa: clasificar -> elegir referencia -> crecer por
    proximidad.

    exentarReferenciaDelBorde=False
        El veto de borde vale para TODAS las islas. Es la logica de V5.
        Resuelve el caso de Juan.

    exentarReferenciaDelBorde=True
        La referencia puede tocar el borde. Modo de rescate para la TC 175.

    Devuelve (referencia, candidatas, totalCM3, densidadDesactivada).
    """
    aplicarDensidad = True
    _clasificar(filas, cfg, aplicarDensidad=True)

    def elegibles(conDensidad):
        return [f for f in filas
                if _esHuesoPlausible(f, cfg, conDensidad)
                and (exentarReferenciaDelBorde or not f["tocaBorde"])]

    plausibles = elegibles(True)
    densidadDesactivada = False
    if not plausibles:
        # Salvaguarda: hueso muy poco mineralizado. Se apaga la densidad y se
        # avisa, antes que quedarse sin craneo.
        aplicarDensidad = False
        _clasificar(filas, cfg, aplicarDensidad=False)
        plausibles = elegibles(False)
        densidadDesactivada = bool(plausibles)

    if not plausibles:
        return None, [], 0.0, densidadDesactivada

    referencia = max(plausibles, key=lambda f: f["vol"])
    referencia["referencia"] = True
    if referencia["motivo"] == "toca el borde del FOV o del array":
        referencia["motivo"] = None
        referencia["protegida"] = True

    candidatas = crecerPorProximidad(labels, filas, referencia, espaciadoZYX, cfg, cache)
    total = sum(f["vol"] for f in candidatas)
    return referencia, candidatas, total, densidadDesactivada


def crecerPorProximidad(labels, filas, referencia, espaciadoZYX, cfg, cache):
    """
    Incorpora al craneo las islas cercanas a la masa aceptada.

    Reemplaza al distance_transform_edt de V5 (ver [P2]). El cache de puntos
    y arboles se comparte entre las dos pasadas, asi el rescate no cuesta el
    doble.
    """
    pool = [f for f in filas if f["motivo"] is None]
    if len(pool) <= 1:
        return pool

    for f in pool:
        k = f["etiqueta"]
        if k not in cache:
            pts = _puntosDeSuperficie(labels, k, f["caja"], espaciadoZYX,
                                      cfg["MAX_PUNTOS_KDTREE"])
            cache[k] = (pts, cKDTree(pts))

    etiquetas = [f["etiqueta"] for f in pool]
    D = {a: {} for a in etiquetas}
    for i, a in enumerate(etiquetas):
        for b in etiquetas[i + 1:]:
            d, _ = cache[b][1].query(cache[a][0], k=1)
            dmin = float(d.min())
            D[a][b] = dmin
            D[b][a] = dmin

    porEtiqueta = {f["etiqueta"]: f for f in pool}
    aceptadas = {referencia["etiqueta"]}
    pendientes = set(etiquetas) - aceptadas
    ronda = 0
    while pendientes:
        ronda += 1
        nuevas = {}
        for b in pendientes:
            dmin = min(D[b][a] for a in aceptadas)
            if dmin <= cfg["MARGEN_PROXIMIDAD_MM"]:
                nuevas[b] = dmin
        if not nuevas:
            break
        for b, dmin in nuevas.items():
            porEtiqueta[b]["dist"] = dmin
            aceptadas.add(b)
            pendientes.discard(b)
        print("CranioPlan:   ronda %d: +%d isla(s) por proximidad." % (ronda, len(nuevas)))

    for b in pendientes:
        dmin = min(D[b][a] for a in aceptadas)
        porEtiqueta[b]["dist"] = dmin
        porEtiqueta[b]["motivo"] = "lejos de la masa principal (%.0f mm)" % dmin

    return [f for f in filas if f["motivo"] is None]


# ------------------------------------------------------------
# Glue con Slicer
# ------------------------------------------------------------

def _config():
    return {
        "HU_MIN": HU_MIN,
        "HU_MAX": HU_MAX,
        "VOL_MINIMO_CM3": VOL_MINIMO_CM3,
        "MARGEN_PROXIMIDAD_MM": MARGEN_PROXIMIDAD_MM,
        "VOLUMEN_CRANEO_MINIMO_CM3": VOLUMEN_CRANEO_MINIMO_CM3,
        "UMBRAL_DENSIDAD_HU": UMBRAL_DENSIDAD_HU,
        "UMBRAL_DENSIDAD_DUDOSA_HU": UMBRAL_DENSIDAD_DUDOSA_HU,
        "ESPESOR_LOSA_MM": ESPESOR_LOSA_MM,
        "FRACCION_LOSA_FOV": FRACCION_LOSA_FOV,
        "MAX_PUNTOS_KDTREE": MAX_PUNTOS_KDTREE,
    }


def _borrarNodosPorNombre(nombre):
    """
    Borra TODOS los nodos que ya se llamen asi.

    Es imprescindible con nombres fijos: si queda un 'Craneo_Automatico' de
    una corrida anterior, Slicer bautiza al nuevo 'Craneo_Automatico_1' y
    despues slicer.util.getNode('Craneo_Automatico') del Corte V9 agarra el
    VIEJO. El sintoma es un corte que parece no hacer nada.
    """
    nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    while nodo is not None:
        slicer.mrmlScene.RemoveNode(nodo)
        nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)


def _forzarSuperficie(segNode):
    """
    Rehace la superficie cerrada DESDE CERO.

    CreateClosedSurfaceRepresentation() no reconvierte si la representacion ya
    existe: despues de tocar el labelmap puede devolver la malla vieja. Hay
    que borrarla primero. Este es el mismo mecanismo que hacia que, tras
    fusionar, el modelo exportado fuera el de una sola pieza.
    """
    try:
        segNode.GetSegmentation().SetConversionParameter(
            "Smoothing factor", str(SUAVIZADO_SUPERFICIE))
    except Exception:
        pass
    try:
        segNode.RemoveClosedSurfaceRepresentation()
    except Exception:
        pass
    segNode.CreateClosedSurfaceRepresentation()


def listar_volumenes():
    """Lista los volumenes escalares de la escena, con dimensiones y spacing."""
    nodos = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
    if not nodos:
        print("CranioPlan: no hay volumenes cargados.")
        return []
    print("CranioPlan: volumenes en la escena:")
    for i, n in enumerate(nodos):
        dims = n.GetImageData().GetDimensions() if n.GetImageData() else (0, 0, 0)
        esp = n.GetSpacing()
        print("  [%d] %-45s  %dx%dx%d  %.2f/%.2f/%.2f mm"
              % (i, n.GetName(), dims[0], dims[1], dims[2], esp[0], esp[1], esp[2]))
    return nodos


def usar_volumen(referencia=None):
    """Elige el volumen de trabajo. Acepta el indice de listar_volumenes(),
    un fragmento del nombre, o un nodo. Sin argumento usa el ultimo cargado."""
    nodos = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
    if not nodos:
        print("CranioPlan: no hay volumenes cargados.")
        return None
    elegido = None
    if referencia is None:
        elegido = nodos[-1]
    elif isinstance(referencia, int):
        if 0 <= referencia < len(nodos):
            elegido = nodos[referencia]
    elif isinstance(referencia, str):
        exactos = [n for n in nodos if n.GetName() == referencia]
        parciales = [n for n in nodos if referencia.lower() in n.GetName().lower()]
        elegido = exactos[0] if exactos else (parciales[0] if parciales else None)
    else:
        elegido = referencia
    if elegido is None:
        print("CranioPlan: no encontre ese volumen. Usa listar_volumenes().")
        return None
    _estado["volumeNode"] = elegido
    print("CranioPlan: volumen de trabajo -> %s" % elegido.GetName())
    return elegido


def _limpiarEscenaPrevia():
    for clave in ("segNode", "modelo"):
        nodo = _estado.get(clave)
        if nodo is not None:
            try:
                slicer.mrmlScene.RemoveNode(nodo)
            except Exception:
                pass
        _estado[clave] = None
    # Y tambien por nombre, por si quedaron de una corrida anterior o de otro
    # script: con nombres fijos, un duplicado rompe el Corte V9 en silencio.
    _borrarNodosPorNombre(NOMBRE_NODO_SEGMENTACION)
    _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)
    _estado["colores"] = {}
    _estado["candidatas"] = []


_PALETA = [
    (0.85, 0.80, 0.65), (0.55, 0.75, 0.85), (0.75, 0.85, 0.55),
    (0.85, 0.55, 0.75), (0.95, 0.75, 0.45), (0.65, 0.65, 0.90),
    (0.90, 0.65, 0.55), (0.60, 0.85, 0.75),
]


def generar():
    """BLOQUE A v6. Segmenta, filtra y deja las piezas candidatas listas para
    la revision manual. No fusiona nada: eso lo hace confirmar_craneo()."""
    print("=" * 66)
    print("CranioPlan - Bloque A %s" % CRANIOPLAN_BLOQUE_A)

    volumeNode = _estado.get("volumeNode")
    if volumeNode is None or slicer.mrmlScene.GetNodeByID(volumeNode.GetID()) is None:
        volumeNode = usar_volumen(NOMBRE_VOLUMEN)
    if volumeNode is None:
        return None

    dims = volumeNode.GetImageData().GetDimensions()
    esp = volumeNode.GetSpacing()
    print("CranioPlan: volumen   : %s" % volumeNode.GetName())
    print("CranioPlan: dims      : %d x %d x %d voxels" % dims)
    print("CranioPlan: spacing   : %.3f x %.3f x %.3f mm" % esp)
    print("CranioPlan: FOV       : %.1f x %.1f x %.1f mm"
          % (dims[0] * esp[0], dims[1] * esp[1], dims[2] * esp[2]))

    _limpiarEscenaPrevia()

    volArr = slicer.util.arrayFromVolume(volumeNode)
    espaciadoZYX = (esp[2], esp[1], esp[0])
    cfg = _config()

    labels, filas, info = analizarVolumen(volArr, espaciadoZYX, cfg)
    if labels is None or not filas:
        print("CranioPlan: no se encontro ninguna isla en el rango %d-%d HU."
              % (HU_MIN, HU_MAX))
        return None

    print("CranioPlan: islas detectadas: %d" % info["nIslas"])
    if info["valorRelleno"] is not None:
        print("CranioPlan: relleno fuera del FOV detectado en %d HU "
              "(se usa para el borde circular)." % info["valorRelleno"])
    else:
        print("CranioPlan: no pude detectar el relleno fuera del FOV; el borde "
              "circular no se testea, solo las caras del array.")

    # --- Seleccion en dos pasadas (ver [P7]) ---
    cache = {}
    referencia, candidatas, total, densOff = seleccionarCandidatas(
        labels, filas, espaciadoZYX, cfg, exentarReferenciaDelBorde=False, cache=cache)

    print("CranioPlan: pasada 1 (veto de borde estricto, para todas las islas): "
          "%s, %.1f cm3 aceptados."
          % ("sin referencia" if referencia is None
             else "referencia %.2f cm3" % referencia["vol"], total))

    if referencia is None or total < cfg["VOLUMEN_CRANEO_MINIMO_CM3"]:
        print("CranioPlan: eso no alcanza para ser un craneo (minimo %.1f cm3). "
              "Se pasa al modo de RESCATE: la referencia queda exenta del veto "
              "de borde, porque lo mas probable es que el FOV este recortando "
              "el craneo (caso TC 175)."
              % cfg["VOLUMEN_CRANEO_MINIMO_CM3"])
        ref2, cand2, total2, densOff2 = seleccionarCandidatas(
            labels, filas, espaciadoZYX, cfg, exentarReferenciaDelBorde=True, cache=cache)
        if ref2 is not None and total2 > total:
            referencia, candidatas, total, densOff = ref2, cand2, total2, densOff2
            print("CranioPlan: rescate aplicado -> referencia %.2f cm3, "
                  "%.1f cm3 aceptados." % (referencia["vol"], total))
        else:
            print("CranioPlan: el rescate no mejoro nada; se deja la pasada 1.")
            seleccionarCandidatas(labels, filas, espaciadoZYX, cfg,
                                  exentarReferenciaDelBorde=False, cache=cache)

    if referencia is None:
        print("CranioPlan: ninguna isla parece hueso. Revisa los umbrales HU o "
              "el volumen elegido.")
        _imprimirTabla(filas)
        return None
    if densOff:
        print("CranioPlan: ATENCION - el filtro de densidad descartaba todo; se "
              "desactivo solo. Revisa UMBRAL_DENSIDAD_HU (hueso poco "
              "mineralizado?).")
    if referencia["protegida"]:
        print("CranioPlan: la referencia toca el borde del FOV y aun asi se "
              "conserva: es el modo de rescate, no el comportamiento normal.")

    _imprimirTabla(filas)

    if not candidatas:
        print("CranioPlan: no quedo ninguna candidata.")
        return None

    candidatas = sorted(candidatas, key=lambda f: f["vol"], reverse=True)

    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_SEGMENTACION)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentacion = segNode.GetSegmentation()

    lista = []
    for i, fila in enumerate(candidatas, start=1):
        nombre = "Pieza_%d" % i
        segId = segmentacion.AddEmptySegment("", nombre)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            (labels == fila["etiqueta"]).astype(np.uint8), segNode, segId, volumeNode)
        color = _PALETA[(i - 1) % len(_PALETA)]
        segmentacion.GetSegment(segId).SetColor(*color)
        _estado["colores"][segId] = color

        notas = []
        if fila["p95"] is not None and fila["p95"] < UMBRAL_DENSIDAD_DUDOSA_HU:
            notas.append("densidad dudosa (p95=%.0f HU)" % fila["p95"])
        if fila["dist"] > MARGEN_PROXIMIDAD_MM * 0.5:
            notas.append("a %.0f mm de la masa principal" % fila["dist"])
        if fila["protegida"]:
            notas.append("toca el borde, conservada por ser la masa principal")

        lista.append({
            "numero": i, "segId": segId, "vol": fila["vol"],
            "dist": fila["dist"], "p95": fila["p95"],
            "nota": "; ".join(notas), "referencia": fila["referencia"],
        })

    _forzarSuperficie(segNode)

    _estado["labels"] = labels
    _estado["espaciadoZYX"] = espaciadoZYX
    _estado["filas"] = filas
    _estado["candidatas"] = lista
    _estado["segNode"] = segNode

    print("")
    print("--- PIEZAS CANDIDATAS (revisar antes de confirmar) ---")
    for d in lista:
        p95 = "n/d" if d["p95"] is None else "%.0f" % d["p95"]
        marca = " *referencia*" if d["referencia"] else ""
        extra = ("   <- " + d["nota"]) if d["nota"] else ""
        print("  [%d] vol=%7.2f cm3   p95=%5s HU   dist=%4.0f mm%s%s"
              % (d["numero"], d["vol"], p95, d["dist"], marca, extra))
    print("-" * 54)
    ayuda()
    return lista


def _imprimirTabla(filas):
    print("")
    print("--- ISLAS DETECTADAS (todas, con motivo de descarte) ---")
    print("   vol(cm3)  p95(HU)  espMin(mm)  %FOV  borde  estado")
    for f in sorted(filas, key=lambda f: f["vol"], reverse=True):
        if f["vol"] < VOL_MINIMO_CM3 and f["motivo"]:
            continue  # no llenamos la consola con motitas
        p95 = "n/d" if f["p95"] is None else "%.0f" % f["p95"]
        estado = "CANDIDATA" if f["motivo"] is None else f["motivo"]
        if f["referencia"]:
            estado += "  *referencia*"
        print("  %9.2f  %7s  %10.1f  %4.0f%%  %-5s  %s"
              % (f["vol"], p95, min(f["extMM"]), max(f["fracFOV"]) * 100,
                 "si" if f["tocaBorde"] else "no", estado))
    ocultas = sum(1 for f in filas if f["vol"] < VOL_MINIMO_CM3 and f["motivo"])
    if ocultas:
        print("  (%d isla(s) por debajo de %.2f cm3 no se listan)"
              % (ocultas, VOL_MINIMO_CM3))
    print("-" * 54)


def resaltar(numero):
    """Pinta esa pieza de rojo; el resto vuelve a su color."""
    segNode = _estado.get("segNode")
    if segNode is None:
        print("CranioPlan: primero corre generar().")
        return
    segmentacion = segNode.GetSegmentation()
    for d in _estado["candidatas"]:
        seg = segmentacion.GetSegment(d["segId"])
        if seg is None:
            continue
        if d["numero"] == numero:
            seg.SetColor(1.0, 0.15, 0.15)
        else:
            seg.SetColor(*_estado["colores"].get(d["segId"], (0.6, 0.6, 0.6)))


def mostrar_todas():
    """Restaura los colores originales de todas las piezas."""
    segNode = _estado.get("segNode")
    if segNode is None:
        return
    segmentacion = segNode.GetSegmentation()
    for d in _estado["candidatas"]:
        seg = segmentacion.GetSegment(d["segId"])
        if seg is not None:
            seg.SetColor(*_estado["colores"].get(d["segId"], (0.6, 0.6, 0.6)))


def eliminar(numero):
    """Borra esa pieza de la revision."""
    segNode = _estado.get("segNode")
    if segNode is None:
        print("CranioPlan: primero corre generar().")
        return
    d = next((x for x in _estado["candidatas"] if x["numero"] == numero), None)
    if d is None:
        print("CranioPlan: no existe la pieza %s." % numero)
        return
    segNode.GetSegmentation().RemoveSegment(d["segId"])
    _estado["candidatas"] = [x for x in _estado["candidatas"] if x["numero"] != numero]
    print("CranioPlan: pieza %d eliminada. Quedan %d."
          % (numero, len(_estado["candidatas"])))


def confirmar_craneo():
    """
    Fusiona las piezas que quedaron en un unico segmento, con el nombre que
    espera el Corte V9 (NOMBRE_SEGMENTO_CRANEO).

    Se hace con un OR de mascaras en numpy y UNA sola escritura, no con el
    efecto "Logical operators" del Segment Editor: ese efecto arrastra las
    reglas de masking del editor (estado oculto) y ademas dejaba viva la
    representacion de superficie vieja de cada segmento.
    """
    segNode = _estado.get("segNode")
    volumeNode = _estado.get("volumeNode")
    if segNode is None or volumeNode is None:
        print("CranioPlan: primero corre generar().")
        return False
    if not _estado["candidatas"]:
        print("CranioPlan: no queda ninguna pieza.")
        return False

    segmentacion = segNode.GetSegmentation()
    esp = volumeNode.GetSpacing()
    volVoxCM3 = (esp[0] * esp[1] * esp[2]) / 1000.0

    union = None
    sumaPiezas = 0.0
    leidas = 0
    ids = [d["segId"] for d in _estado["candidatas"]
           if segmentacion.GetSegment(d["segId"]) is not None]
    for segId in ids:
        arr = slicer.util.arrayFromSegmentBinaryLabelmap(segNode, segId, volumeNode)
        if arr is None:
            print("CranioPlan: ADVERTENCIA - no pude leer %s, queda afuera." % segId)
            continue
        arr = arr.astype(bool)
        sumaPiezas += int(np.count_nonzero(arr)) * volVoxCM3
        union = arr.copy() if union is None else np.logical_or(union, arr)
        leidas += 1

    if union is None or not union.any():
        print("CranioPlan: la union quedo vacia.")
        return False

    volUnion = int(np.count_nonzero(union)) * volVoxCM3

    base = ids[0]
    slicer.util.updateSegmentBinaryLabelmapFromArray(
        union.astype(np.uint8), segNode, base, volumeNode)
    for i in range(segmentacion.GetNumberOfSegments() - 1, -1, -1):
        segId = segmentacion.GetNthSegmentID(i)
        if segId != base:
            segmentacion.RemoveSegment(segId)
    segmentacion.GetSegment(base).SetName(NOMBRE_SEGMENTO_CRANEO)
    segmentacion.GetSegment(base).SetColor(0.9, 0.8, 0.6)

    _forzarSuperficie(segNode)

    verif = slicer.util.arrayFromSegmentBinaryLabelmap(segNode, base, volumeNode)
    volNodo = (int(np.count_nonzero(verif)) * volVoxCM3) if verif is not None else float('nan')
    estructura = ndimage.generate_binary_structure(3, 3)
    _, piezasVoxel = ndimage.label(union, structure=estructura)

    poly = vtk.vtkPolyData()
    segNode.GetClosedSurfaceRepresentation(base, poly)
    regiones = _contarRegiones(poly)

    print("")
    print("--- CONFIRMACION ---")
    print("  piezas fusionadas       : %d de %d" % (leidas, len(ids)))
    print("  suma de las piezas      : %.2f cm3" % sumaPiezas)
    print("  union (sin solapamiento): %.2f cm3" % volUnion)
    print("  escrito en el segmento  : %.2f cm3" % volNodo)
    print("  piezas conexas (voxels) : %d" % piezasVoxel)
    print("  regiones en la malla    : %d" % regiones)
    print("--------------------")
    if regiones < piezasVoxel:
        print("CranioPlan: ADVERTENCIA - la malla tiene menos regiones que los "
              "voxels. Alguna placa fina se disolvio en el suavizado: bajar "
              "SUAVIZADO_SUPERFICIE a '0.1' y repetir.")
    print("CranioPlan: listo. Segmento '%s' dentro del nodo '%s'. "
          "Ahora: enviar_a_planner()"
          % (NOMBRE_SEGMENTO_CRANEO, NOMBRE_NODO_SEGMENTACION))
    return True


def _contarRegiones(polyData):
    if polyData is None or polyData.GetNumberOfPoints() == 0:
        return 0
    c = vtk.vtkPolyDataConnectivityFilter()
    c.SetInputData(polyData)
    c.SetExtractionModeToAllRegions()
    c.Update()
    return int(c.GetNumberOfExtractedRegions())


def _limpiarRuidoMalla(polyData, minimoPuntos):
    """
    Quita SOLO el ruido de marching cubes, con umbral ABSOLUTO en cantidad de
    puntos.

    V5 usaba un umbral RELATIVO al 10% de la region mayor y conservaba UNA
    sola region de 454: por eso 'Craneo_Final' aparecia sin mandibula. Con un
    umbral relativo, cualquier placa craneal desconectada por una sutura
    abierta que mida menos del 10% de la boveda se borra sola.
    """
    if polyData is None or polyData.GetNumberOfPoints() == 0:
        return polyData

    con = vtk.vtkPolyDataConnectivityFilter()
    con.SetInputData(polyData)
    con.SetExtractionModeToAllRegions()
    con.ColorRegionsOn()
    con.Update()
    nRegiones = int(con.GetNumberOfExtractedRegions())
    if nRegiones <= 1:
        limpiar = vtk.vtkCleanPolyData()
        limpiar.SetInputData(polyData)
        limpiar.Update()
        salida = vtk.vtkPolyData()
        salida.DeepCopy(limpiar.GetOutput())
        return salida

    arr = con.GetOutput().GetPointData().GetArray("RegionId")
    if arr is None:
        return polyData
    conteo = {}
    for i in range(arr.GetNumberOfTuples()):
        rid = int(arr.GetTuple1(i))
        conteo[rid] = conteo.get(rid, 0) + 1

    conservar = [rid for rid, c in conteo.items() if c >= minimoPuntos]
    descartadas = sorted([c for rid, c in conteo.items() if c < minimoPuntos],
                         reverse=True)
    if not conservar:
        conservar = [max(conteo, key=conteo.get)]
        descartadas = []

    con.SetExtractionModeToSpecifiedRegions()
    con.InitializeSpecifiedRegionList()
    for rid in conservar:
        con.AddSpecifiedRegion(rid)
    con.Update()
    limpiar = vtk.vtkCleanPolyData()
    limpiar.SetInputConnection(con.GetOutputPort())
    limpiar.Update()
    salida = vtk.vtkPolyData()
    salida.DeepCopy(limpiar.GetOutput())

    print("CranioPlan: limpieza de malla: %d region(es) -> %d conservada(s) "
          "(umbral ABSOLUTO %d puntos)." % (nRegiones, len(conservar), minimoPuntos))
    if descartadas:
        print("CranioPlan:   descartadas (puntos): %s%s"
              % (descartadas[:12], " ..." if len(descartadas) > 12 else ""))
        if descartadas[0] >= minimoPuntos * 0.5:
            print("CranioPlan:   ATENCION - la mayor descartada esta cerca del "
                  "umbral. Si falta hueso, bajar MINIMO_PUNTOS_MALLA.")
    return salida


def _decimarMalla(polyData, reduccion):
    """Decima verificando que no desaparezcan regiones conexas. Si la cuenta
    baja, reintenta con la mitad y si sigue bajando devuelve la malla sin
    decimar: preferimos render pesado antes que perder hueso."""
    if not reduccion or reduccion <= 0.0 or polyData is None:
        return polyData
    if polyData.GetNumberOfCells() == 0:
        return polyData
    antes = _contarRegiones(polyData)
    intento = float(reduccion)
    for _ in range(2):
        d = vtk.vtkDecimatePro()
        d.SetInputData(polyData)
        d.SetTargetReduction(intento)
        d.PreserveTopologyOn()
        d.Update()
        res = vtk.vtkPolyData()
        res.DeepCopy(d.GetOutput())
        if res.GetNumberOfPoints() == 0:
            break
        if _contarRegiones(res) >= antes:
            return res
        print("CranioPlan: la decimacion al %.0f%% perdio regiones; reintento."
              % (intento * 100))
        intento /= 2.0
    print("CranioPlan: se devuelve la malla SIN decimar para no perder piezas.")
    return polyData


def enviar_a_planner():
    """
    Exporta el segmento del craneo a un Model node con el nombre que espera
    el Corte V9 (NOMBRE_MODELO_CRANEO).

    Se fuerza la reconversion de la superficie antes de leerla: si el labelmap
    cambio despues de que la representacion existia,
    CreateClosedSurfaceRepresentation() NO reconvierte y se lee una malla
    vieja.
    """
    segNode = _estado.get("segNode")
    if segNode is None:
        print("CranioPlan: primero corre generar() y confirmar_craneo().")
        return None
    segmentacion = segNode.GetSegmentation()
    segId = segmentacion.GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_CRANEO)
    if not segId:
        print("CranioPlan: todavia no existe el segmento '%s'. Corre "
              "confirmar_craneo()." % NOMBRE_SEGMENTO_CRANEO)
        return None

    _forzarSuperficie(segNode)

    bruto = vtk.vtkPolyData()
    segNode.GetClosedSurfaceRepresentation(segId, bruto)
    if bruto.GetNumberOfPoints() == 0:
        print("CranioPlan: la conversion a superficie no devolvio geometria.")
        return None

    regBruto, celBruto = _contarRegiones(bruto), bruto.GetNumberOfCells()
    malla = _limpiarRuidoMalla(bruto, MINIMO_PUNTOS_MALLA)
    regLimpio, celLimpio = _contarRegiones(malla), malla.GetNumberOfCells()
    malla = _decimarMalla(malla, REDUCCION_MALLA)

    modelo = _estado.get("modelo")
    if modelo is not None:
        try:
            slicer.mrmlScene.RemoveNode(modelo)
        except Exception:
            pass
    # Nombre EXACTO, sin GenerateUniqueName: si Slicer le agrega un sufijo, el
    # Corte V9 no lo encuentra (o encuentra el modelo viejo).
    _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)
    modelo = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLModelNode', NOMBRE_MODELO_CRANEO)
    modelo.SetAndObservePolyData(malla)
    modelo.CreateDefaultDisplayNodes()
    modelo.GetDisplayNode().SetColor(0.9, 0.8, 0.6)
    modelo.GetDisplayNode().SetScalarVisibility(False)
    _estado["modelo"] = modelo

    display = segNode.GetDisplayNode()
    if display is not None:
        display.SetVisibility(False)

    print("")
    print("--- EXPORTACION A MODELO ---")
    print("  bruto  : %7d triangulos, %3d region(es)" % (celBruto, regBruto))
    print("  limpio : %7d triangulos, %3d region(es)" % (celLimpio, regLimpio))
    print("  final  : %7d triangulos, %3d region(es)"
          % (malla.GetNumberOfCells(), _contarRegiones(malla)))
    print("----------------------------")
    if modelo.GetName() != NOMBRE_MODELO_CRANEO:
        print("CranioPlan: ATENCION - el modelo quedo como '%s' y no como '%s'. "
              "Habia otro nodo con ese nombre: borralo y repeti, o el Corte V9 "
              "va a usar el equivocado." % (modelo.GetName(), NOMBRE_MODELO_CRANEO))
    else:
        print("CranioPlan: modelo '%s' creado." % modelo.GetName())
    print("CranioPlan: se oculto la segmentacion para no renderizar dos veces "
          "el mismo craneo (causa tipica del lag al rotar).")
    print("CranioPlan: ya podes dibujar las curvas y correr cortar() del "
          "Corte V9.")
    return modelo


def verificar_compatibilidad_v9():
    """
    Chequea que la escena tenga EXACTAMENTE lo que el Corte V9 va a buscar.
    Corrarlo antes de cortar() ahorra el 'no encuentro el nodo...'.
    """
    print("")
    print("--- COMPATIBILIDAD CON EL CORTE V9 ---")
    ok = True

    segs = [n for n in slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
            if n.GetName() == NOMBRE_NODO_SEGMENTACION]
    if len(segs) == 1:
        print("  [OK]    nodo de segmentacion '%s'" % NOMBRE_NODO_SEGMENTACION)
        sid = segs[0].GetSegmentation().GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_CRANEO)
        if sid:
            print("  [OK]    segmento '%s' adentro" % NOMBRE_SEGMENTO_CRANEO)
        else:
            print("  [FALTA] segmento '%s': corre confirmar_craneo()"
                  % NOMBRE_SEGMENTO_CRANEO)
            ok = False
    elif not segs:
        print("  [FALTA] nodo de segmentacion '%s': corre generar()"
              % NOMBRE_NODO_SEGMENTACION)
        ok = False
    else:
        print("  [ERROR] hay %d nodos llamados '%s'. Borra los sobrantes."
              % (len(segs), NOMBRE_NODO_SEGMENTACION))
        ok = False

    mods = [n for n in slicer.util.getNodesByClass("vtkMRMLModelNode")
            if n.GetName() == NOMBRE_MODELO_CRANEO]
    if len(mods) == 1:
        print("  [OK]    Model node '%s'" % NOMBRE_MODELO_CRANEO)
    elif not mods:
        print("  [FALTA] Model node '%s': corre enviar_a_planner()"
              % NOMBRE_MODELO_CRANEO)
        ok = False
    else:
        print("  [ERROR] hay %d modelos llamados '%s'. Borra los sobrantes."
              % (len(mods), NOMBRE_MODELO_CRANEO))
        ok = False

    nCer = len(slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode'))
    nAbi = len([n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')])
    print("  curvas: %d cerrada(s), %d abierta(s)" % (nCer, nAbi))
    if nCer + nAbi == 0:
        print("  [FALTA] dibuja al menos una curva de corte.")
        ok = False

    print("  -> %s" % ("todo listo para cortar()" if ok else "falta algo, ver arriba"))
    print("--------------------------------------")
    return ok


def ayuda():
    print("")
    print("Comandos disponibles:")
    print("  listar_volumenes()          -> lista las series cargadas")
    print("  usar_volumen(i)             -> elige serie por indice o nombre")
    print("  generar()                   -> corre el Bloque A sobre esa serie")
    print("  resaltar(n)                 -> pinta esa pieza de rojo")
    print("  eliminar(n)                 -> borra esa pieza de la revision")
    print("  mostrar_todas()             -> restaura colores")
    print("  confirmar_craneo()          -> fusiona en el segmento '%s'"
          % NOMBRE_SEGMENTO_CRANEO)
    print("  enviar_a_planner()          -> exporta el modelo '%s'"
          % NOMBRE_MODELO_CRANEO)
    print("  verificar_compatibilidad_v9() -> chequea la escena antes de cortar()")
    print("")
    print("Nodos que quedan (los que busca el Corte V9):")
    print("  segmentacion = '%s'  (segmento '%s')"
          % (NOMBRE_NODO_SEGMENTACION, NOMBRE_SEGMENTO_CRANEO))
    print("  modelo       = '%s'" % NOMBRE_MODELO_CRANEO)


generar()
