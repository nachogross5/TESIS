# ============================================================
# CRANIOPLAN - BLOQUE A v7.2  (script de consola para 3D Slicer)
#            COMPATIBLE CON EL CORTE V9
# ============================================================
# Se pega entero en la consola de Python de Slicer. Al pegarlo se ejecuta.
#
# ============================================================
# DE DONDE SALEN LOS UMBRALES
# ============================================================
# De 206 islas etiquetadas a mano sobre 9 series de 7 pacientes del Hospital
# Garrahan: JUAN 702, RUIZ 175, AVILA 403, LUNA 167, JUSTINA 559/381 (dos
# reconstrucciones), MESKY 912, MORENO 694 (preoperatorio y postoperatorio).
# El relevamiento se hizo con BloqueA_v72_DIAGNOSTICO.py, que mide lo mismo
# pero no elimina nada.
#
# ------------------------------------------------------------
# QUE CAMBIA RESPECTO DE v7.1
# ------------------------------------------------------------
# [1] LA REFERENCIA SE ELIGE POR CORTICAL + COMPACIDAD + VOLUMEN, y se saco
#     la cota fija de 180 mm del camino critico.
#
#     El problema de la cota fija: el craneo mas grande que medimos (MESKY)
#     llega a 166 mm. Con la cota en 180 quedaban 14 mm de holgura, un 8 %.
#     Si llega un chico de 3-4 anios con escafocefalia diagnosticada tarde y
#     su craneo pasa los 180 mm, el craneo queda EXCLUIDO de ser referencia y
#     todo se derrumba: gana otra cosa y encima el criterio de extension
#     empieza a borrar el craneo por "muy grande".
#
#     La compacidad (volumen / eje_mayor^3, ver _compacidad) resuelve lo
#     mismo sin depender de la edad, porque es adimensional. Craneos medidos:
#     24 a 78. Impostores: almohadilla de JUAN 11, riel de camilla de
#     JUSTINA601 1. El hueco entre 11 y 24 es el margen.
#
# [2] EL VETO DE BORDE DEL FOV DESAPARECIO POR COMPLETO. En v7.1 todavia se
#     usaba para elegir la referencia en una primera pasada, con un modo de
#     rescate para cuando eso dejaba al craneo afuera (caso RUIZ). Lo unico
#     que ese veto resolvia era JUAN, y ahora lo resuelve la compacidad. Se
#     fueron las dos pasadas y el caso especial. Tocar el borde ahora solo
#     genera una NOTA DE CALIDAD para el cirujano: el craneo puede venir
#     recortado.
#
# [3] EL FILTRO DE MATERIAL YA NO ES BORRADO DURO. En MORENO aparecio hueso
#     craneal real con 6.0 % de cortical contra un umbral de 5 %, y el
#     chupete de JUAN tiene 3.88 %. Los rangos se tocan: no existe ningun
#     umbral que deje el hueso adentro y el chupete afuera. Entonces una isla
#     de cortical baja que este a <= DIST_ACEPTACION_MM de la referencia pasa
#     a DUDOSA en vez de eliminarse. Cuesta un click y no puede perder una
#     lamina de hueso pegada al craneo.
#
# [4] ELIMINACION POR EXTENSION RELATIVA. Una isla que no es la referencia y
#     mide mas de FACTOR_EXTENSION_MAXIMA veces el eje mayor de la referencia
#     no puede ser un pedazo de craneo: un fragmento no es mas grande que el
#     craneo del que salio. Limite relativo, asi que escala con el paciente.
#     Peor hueso no-referencia medido: 0.50x la referencia, contra un limite
#     de 1.25x.
#
# [5] REFERENCIA VISIBLE Y FORZABLE. Se imprime la tabla de candidatas con
#     sus metricas, hay usar_referencia(n) para corregir a mano, y un aviso
#     fuerte cuando se elimina por extension una isla grande y densa (que es
#     justo el caso donde, si el algoritmo se equivoco, se equivoco feo).
#
# [6] MODO POSTOPERATORIO. En un craneo ya operado los colgajos estan
#     genuinamente separados: en MORENO 501 quedaron a 2.2, 2.6, 5.1, 11.7 y
#     16.7 mm. Ningun umbral preoperatorio razonable los recupera, y no
#     deberia: son dos problemas distintos. postop() sube la distancia de
#     aceptacion y avisa que la revision manual pasa a ser obligatoria.
#
# [7] COMPACIDAD MINIMA DE UNA PIEZA. Verificando el v7.2 sobre los 6 casos
#     aparecio ruido de reconstruccion que quedaba DUDOSO: nubes de motas
#     alrededor del craneo, con cortical casi nula pero pegadas a el, asi que
#     la regla blanda de material [3] las conservaba. En JUAN son dos islas de
#     3.56 y 0.80 cm3 que se ven como salpicado naranja sobre toda la calota.
#     Se separan por FORMA, no por densidad: compacidad 1.5-2.1 el ruido y los
#     soportes, contra 14.3-157 el hueso real no-referencia en preoperatorio.
#     Regla: isla que no es la referencia con compacidad < 6 se elimina.
#     Desactivada en postoperatorio (alli el colgajo mas flaco da 8.5).
#     Efecto medido: JUAN 14 -> 10 dudosas, AVILA 3 -> 1, JUSTINA601 6 -> 4,
#     RUIZ 6 -> 5, y ademas se va el riel de camilla de JUSTINA.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (documentados a proposito, no ignorados)
# ------------------------------------------------------------
# [C1] LA COMPACIDAD ESTA CALIBRADA CON DOS IMPOSTORES. Siete craneos de un
#      lado (24 a 78) y dos objetos del otro (11 y 1). El hueco es ancho pero
#      n=2. Mitigacion: la eleccion de referencia es visible y forzable.
#
# [C2] EL CHUPETE NO ES SEPARABLE POR MATERIAL. JUAN 3.88 % de cortical,
#      MESKY 0.0 %. El hueso mas marginal medido tiene 6.0 %. Con el umbral
#      en 5 %, el chupete de JUAN cae del lado del hueso y ademas esta a
#      1.69 mm del craneo, asi que en v7.2 queda DUDOSO. Es un click, a
#      cambio de no poder perder una lamina fina.
#
# [C3] LAS VERTEBRAS SON HUESO REAL. Ningun criterio de material ni de
#      distancia las distingue de la mandibula: estan pegadas a la base por
#      los condilos. Entran como aceptadas o dudosas segun el caso. Se borran
#      con eliminar(n) o de una vez con eliminar_dudosas(). NO se eliminan
#      solas a proposito.
#
# [C4] SI LA BOVEDA VIENE PARTIDA EN DOS CASCARAS FINAS y ninguna sola llega
#      a compacidad 16, se relaja ese filtro y se avisa. Nunca se observo,
#      pero la salvaguarda existe.
#
# [C5] NUEVE SERIES DE SIETE PACIENTES, UN SOLO HOSPITAL. Los umbrales son
#      defendibles sobre estos datos y nada mas. Cada caso nuevo conviene
#      pasarlo primero por el diagnostico y sumarlo al CSV.
#
# [C6] LA ENVOLVENCIA NO DECIDE NADA. Quedo descartada como criterio con
#      datos: las almohadillas de MORENO dieron 0.3 a 1.7 %, por debajo de
#      las vertebras, porque venian fragmentadas. Se sigue midiendo por si
#      sirve para otra cosa.
# ============================================================


# ============================================================
# CONFIG
# ============================================================

# --- Nombres de nodos (deben coincidir con la config del Corte V9) ---
NOMBRE_NODO_SEGMENTACION = "Craneo_Automatico"   # V9: NOMBRE_SEGMENTACION
NOMBRE_SEGMENTO_CRANEO   = "Craneo_Final"        # V9: NOMBRE_SEGMENTO_HUESO
NOMBRE_MODELO_CRANEO     = "Craneo_Final"        # V9: NOMBRE_MODELO_CRANEO

NOMBRE_VOLUMEN = None          # None = ultimo volumen escalar cargado

HU_MIN = 300
HU_MAX = 3000

# --- Material (ver [3] y [C2]) ---
UMBRAL_CORTICAL_HU = 700.0
FRACCION_CORTICAL_MINIMA = 0.05     # hueso craneal mas marginal medido: 6.0 %

# --- Forma: compacidad para elegir referencia (ver [1] y [C1]) ---
COMPACIDAD_MINIMA_REFERENCIA = 16.0   # craneos 24-78 | almohadilla 11 | riel 1

# --- Extension relativa para eliminar (ver [4]) ---
FACTOR_EXTENSION_MAXIMA = 1.25        # peor hueso no-referencia: 0.50x

# --- Compacidad minima de una PIEZA (ver [7]) ---
# Una isla que no es la referencia y es difusa -no compacta- no es un pedazo
# de hueso: es ruido de reconstruccion, un riel o una almohadilla.
#   ruido difuso medido ............ 1.5 - 2.1
#   rieles y almohadillas .......... 0.5 - 4.8
#   hueso real no-referencia PREOP . 14.3 - 157
# Se desactiva en modo postoperatorio: alli el colgajo mas flaco de MORENO
# 501 da 8.5 y no queremos un margen de 1.4x sobre hueso real.
COMPACIDAD_MINIMA_PIEZA = 6.0

# --- Distancias ---
DIST_ACEPTACION_MM = 2.0     # a la REFERENCIA. Hueso craneal medido: max 1.70
DIST_DUDOSA_MM     = 15.0    # a la masa aceptada. Mas alla, se elimina.

# --- Modo postoperatorio (ver [6]) ---
MODO_POSTOP = False
DIST_POSTOP_ACEPTACION_MM = 60.0

# --- Pisos y cotas ---
VOL_MINIMO_CM3 = 0.05                 # debajo de esto no se analiza
VOLUMEN_CRANEO_MINIMO_CM3 = 20.0      # solo para avisar si algo salio mal
ESPESOR_LOSA_MM = 30.0
FRACCION_LOSA_FOV = 0.85

# --- Envolvencia: se mide, no decide (ver [C6]) ---
ENV_BINS_COSTHETA = 12
ENV_BINS_PHI      = 24

# --- Exportacion de malla ---
MINIMO_PUNTOS_MALLA = 200      # umbral ABSOLUTO de ruido de marching cubes
REDUCCION_MALLA = 0.7          # 0.0 desactiva la decimacion
SUAVIZADO_SUPERFICIE = "0.3"   # bajar a "0.1" si se disuelven placas finas

MAX_PUNTOS_KDTREE = 40000


# ============================================================
# CODIGO
# ============================================================

import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

CRANIOPLAN_BLOQUE_A = "v7.2a (2026-08-18) consola - calibrado sobre 9 series"

_estado = {
    "volumeNode": None, "labels": None, "espaciadoZYX": None,
    "filas": [], "candidatas": [], "segNode": None, "colores": {},
    "modelo": None, "referencia": None,
}

COLOR_ACEPTADA = (0.90, 0.80, 0.60)
COLOR_DUDOSA   = (1.00, 0.60, 0.10)


# ============================================================
# NUCLEO DE DECISION - identico en el modulo y en el diagnostico
# ============================================================
# Este bloque es byte por byte el mismo en BloqueA_v72_consola.py y en
# BloqueA_v72_DIAGNOSTICO.py. Es a proposito: si el diagnostico usara una
# copia parecida pero no igual, con el tiempo las dos derivarian y el
# relevamiento dejaria de medir lo que el modulo realmente hace.
# ============================================================


def _valorDeRelleno(volArr):
    """Valor con el que el tomografo rellena lo que queda FUERA del circulo
    de reconstruccion. Se lee en las cuatro esquinas de un corte central."""
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
    """El 'afuera': el relleno fuera del circulo de reconstruccion mas las
    seis caras planas del array."""
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
    """Seis desplazamientos en vez de dilatar el volumen completo: la
    dilatacion pediria otro array de volumen entero, que es lo que hacia
    fallar al v5 por memoria.

    OJO: en v7.2 tocar el borde NO elimina ni condiciona nada. Se mide y se
    informa porque es un dato de CALIDAD DEL ESTUDIO util para el cirujano
    (el craneo puede venir recortado), no un criterio de segmentacion."""
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


def _esLosaDeSoporte(extMM, fracFOV):
    """Tabla de la camilla: abarca casi todo el campo en algun eje y es
    delgada en su eje mas fino. La caja del craneo no cumple eso (su eje mas
    fino ronda los 90-150 mm)."""
    if max(extMM) <= 0.0:
        return False
    return bool(min(extMM) <= ESPESOR_LOSA_MM and max(fracFOV) >= FRACCION_LOSA_FOV)


def _compacidad(volCM3, extMM):
    """
    Cuanto volumen ocupa la isla en relacion al cubo de su eje mas largo,
    x1e6 para que de numeros legibles.

    ES EL CRITERIO DE FORMA QUE REEMPLAZA A LA COTA FIJA DE 180 mm.

    Por que sirve: un craneo es una masa COMPACTA, mucho volumen metido en
    su propio tamano. Una almohadilla o un riel de camilla son ALARGADOS: se
    estiran mucho y ocupan poco. Medido sobre 7 craneos y 2 impostores:

        riel de camilla (JUSTINA601) ....  1
        almohadilla envolvente (JUAN) ... 11
        craneo de JUAN .................. 24   <- el mas flaco de los craneos
        craneo de MESKY ................. 46
        craneo de JUSTINA ............... 78

    Y lo importante: la magnitud es ADIMENSIONAL, o sea invariante de escala.
    El craneo de un chico de 4 anios es mas grande pero igual de compacto; la
    almohadilla sigue siendo alargada. Por eso este criterio no se rompe con
    la edad, que era el problema de la cota fija de 180 mm (el craneo mas
    grande que medimos llega a 166 mm: solo 14 mm de holgura).

    LIMITE CONOCIDO: el corte esta calibrado con DOS impostores contra siete
    craneos. El hueco entre 11 y 24 es ancho, pero no es lo mismo que tenerlo
    medido en veinte almohadillas. Por eso la eleccion de referencia queda
    ademas visible y forzable con usar_referencia(n).
    """
    e = max(extMM)
    if e <= 0.0:
        return 0.0
    return float(volCM3) / (e ** 3) * 1.0e6


def _envolvencia(ptsMM, centroMM):
    """
    Fraccion de direcciones del espacio, vistas desde el centro del craneo,
    en las que aparece esta isla. Bins de AREA IGUAL (uniformes en cos(theta)
    y en phi), no grilla lat-lon, que concentraria bins diminutos en los polos.

    NO DECIDE NADA, en ninguna version. Se dejo midiendo porque cuesta poco.
    Como detector de almohadillas quedo DESCARTADA con datos: las tres
    almohadillas de MORENO dieron 0.3 a 1.7 %, por debajo de las vertebras,
    porque venian partidas en fragmentos y un fragmento no envuelve nada. La
    metrica describe la forma de la ISLA, no la del objeto fisico.
    """
    if ptsMM is None or len(ptsMM) == 0:
        return 0.0
    v = np.asarray(ptsMM, dtype=np.float64) - np.asarray(centroMM, dtype=np.float64)
    r = np.linalg.norm(v, axis=1)
    ok = r > 1e-6
    if not ok.any():
        return 0.0
    v, r = v[ok], r[ok]
    cosTheta = v[:, 0] / r
    phi = np.arctan2(v[:, 2], v[:, 1])
    i = np.clip(((cosTheta + 1.0) * 0.5 * ENV_BINS_COSTHETA).astype(np.int32),
                0, ENV_BINS_COSTHETA - 1)
    j = np.clip(((phi + np.pi) / (2.0 * np.pi) * ENV_BINS_PHI).astype(np.int32),
                0, ENV_BINS_PHI - 1)
    ocupados = np.unique(i.astype(np.int64) * ENV_BINS_PHI + j.astype(np.int64))
    return float(len(ocupados)) / float(ENV_BINS_COSTHETA * ENV_BINS_PHI)


def _puntosDeSuperficie(labels, etiqueta, caja, espaciadoZYX, maxPuntos):
    """Coordenadas fisicas (mm) de los voxeles de SUPERFICIE de una isla. La
    distancia minima entre dos solidos se alcanza siempre en sus superficies:
    mismo resultado que la EDT de volumen completo con dos ordenes de
    magnitud menos de puntos."""
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


def analizarVolumen(volArr, espaciadoZYX, huMin, huMax):
    """Etiqueta las islas y calcula metricas. No decide nada."""
    mascara = (volArr >= huMin) & (volArr <= huMax)
    estructura = ndimage.generate_binary_structure(3, 3)   # 26-conexo
    labels, nIslas = ndimage.label(mascara, structure=estructura)
    del mascara
    if nIslas == 0:
        return None, [], {"nIslas": 0, "valorRelleno": None, "totalCM3": 0.0}
    if nIslas < 32000:
        labels = labels.astype(np.int16)

    conteos = np.bincount(labels.ravel())
    cajas = ndimage.find_objects(labels)
    valorRelleno = _valorDeRelleno(volArr)
    exterior = _mascaraExterior(volArr, valorRelleno)
    etiquetasBorde = _etiquetasTocandoExterior(labels, exterior)
    del exterior

    volVoxCM3 = (espaciadoZYX[0] * espaciadoZYX[1] * espaciadoZYX[2]) / 1000.0
    totalCM3 = float(conteos[1:].sum()) * volVoxCM3

    filas = []
    for k in range(1, nIslas + 1):
        caja = cajas[k - 1]
        if caja is None:
            continue
        vol = int(conteos[k]) * volVoxCM3
        if vol < VOL_MINIMO_CM3:
            continue
        extMM = tuple((caja[e].stop - caja[e].start) * espaciadoZYX[e] for e in range(3))
        fracFOV = tuple((caja[e].stop - caja[e].start) / float(volArr.shape[e])
                        for e in range(3))
        sub = (labels[caja] == k)
        valores = volArr[caja][sub]
        idx = np.argwhere(sub)
        filas.append({
            "etiqueta": k, "vol": vol, "caja": caja,
            "extMM": extMM, "fracFOV": fracFOV,
            "p50": float(np.percentile(valores, 50)),
            "p95": float(np.percentile(valores, 95)),
            "p99": float(np.percentile(valores, 99)),
            "max": float(valores.max()),
            "fracCort": float(np.count_nonzero(valores >= UMBRAL_CORTICAL_HU))
                        / float(valores.size),
            "compac": _compacidad(vol, extMM),
            "tocaBorde": k in etiquetasBorde,
            "esLosa": _esLosaDeSoporte(extMM, fracFOV),
            "ijk": (idx[:, 2].mean() + caja[2].start,
                    idx[:, 1].mean() + caja[1].start,
                    idx[:, 0].mean() + caja[0].start),
            "distRef": None, "distMasa": None, "env": None,
            "veredicto": None, "categoria": "", "motivo": "",
            "referencia": False, "admisible": False,
            "etiqueta_real": "", "nota_real": "",
        })
    filas.sort(key=lambda f: f["vol"], reverse=True)
    return labels, filas, {"nIslas": nIslas, "valorRelleno": valorRelleno,
                           "totalCM3": totalCM3}


def candidatasAReferencia(filas):
    """
    Admisibilidad para competir por REFERENCIA. Dos filtros, cada uno tapa
    el punto ciego del otro:

      cortical >= 5 %   -> saca espuma, gel, plastico, ruido y el riel de la
                           camilla de JUSTINA601 (cortical 23.7 % pasaria por
                           poco, pero cae por compacidad).
      compacidad >= 16  -> saca la almohadilla envolvente de JUAN, que es el
                           unico impostor del dataset que le GANA en volumen
                           al craneo (79.5 vs 58.7 cm3). Con el criterio del
                           v6, "mayor volumen", la referencia de JUAN habria
                           sido la almohadilla.

    Entre las admisibles gana la de mayor volumen.

    NOTA: aca ya NO se usa el veto de borde del FOV que tenian v6 y v7.1.
    Era innecesario: lo unico que resolvia era JUAN, y eso ahora lo resuelve
    la compacidad. Y hacia falta un modo de rescate en dos pasadas porque en
    la TC 175 (RUIZ) el craneo SI toca el borde. Sacarlo elimina las dos
    pasadas y arregla RUIZ sin caso especial.
    """
    for f in filas:
        f["admisible"] = (f["fracCort"] >= FRACCION_CORTICAL_MINIMA
                          and f["compac"] >= COMPACIDAD_MINIMA_REFERENCIA)
    return [f for f in filas if f["admisible"]]


def _elegirReferencia(filas, avisar=print):
    """Devuelve (referencia, relajado). 'relajado' dice si hubo que soltar
    algun filtro, y en ese caso hay que mirar la eleccion con desconfianza."""
    cand = candidatasAReferencia(filas)
    if cand:
        return max(cand, key=lambda f: f["vol"]), None

    # Salvaguarda 1: ningun candidato compacto. Puede pasar si la boveda esta
    # partida en dos cascaras finas y ninguna sola es compacta.
    porCortical = [f for f in filas if f["fracCort"] >= FRACCION_CORTICAL_MINIMA]
    if porCortical:
        avisar("CranioPlan: ATENCION - ninguna isla llega a compacidad %.0f. Se "
               "relaja ese filtro y se elige por volumen entre las de cortical "
               "suficiente. REVISA la referencia elegida."
               % COMPACIDAD_MINIMA_REFERENCIA)
        for f in porCortical:
            f["admisible"] = True
        return max(porCortical, key=lambda f: f["vol"]), "compacidad"

    # Salvaguarda 2: hueso muy poco mineralizado o serie de kernel muy blando.
    if filas:
        avisar("CranioPlan: ATENCION - ninguna isla llega a cortical %.0f%%. Se "
               "relajan TODOS los filtros de referencia. Es muy probable que la "
               "serie elegida no sirva (kernel blando?). REVISA todo."
               % (FRACCION_CORTICAL_MINIMA * 100.0))
        for f in filas:
            f["admisible"] = True
        return max(filas, key=lambda f: f["vol"]), "cortical"
    return None, None


def decidir(labels, filas, espaciadoZYX, referenciaForzada=None, avisar=print):
    """
    Asigna veredicto a cada isla. ES LA MISMA FUNCION que usa el diagnostico.

    referenciaForzada: numero de isla (1-based sobre filas ordenadas por
    volumen) para saltear la eleccion automatica.

    Devuelve (referencia, aceptadas, dudosas, relajado).
    """
    for f in filas:
        f.update({"distRef": None, "distMasa": None, "env": None,
                  "veredicto": None, "categoria": "", "motivo": "",
                  "referencia": False, "admisible": False})

    # --- 1. referencia ---
    relajado = None
    if referenciaForzada is not None:
        if not (1 <= referenciaForzada <= len(filas)):
            avisar("CranioPlan: no existe la isla %s." % referenciaForzada)
            return None, [], [], None
        candidatasAReferencia(filas)          # deja marcada la admisibilidad
        referencia = filas[referenciaForzada - 1]
        avisar("CranioPlan: referencia FORZADA a mano -> isla %d (%.2f cm3)."
               % (referenciaForzada, referencia["vol"]))
    else:
        referencia, relajado = _elegirReferencia(filas, avisar=avisar)
        if referencia is None:
            return None, [], [], None
    referencia["referencia"] = True

    distAceptacion = (DIST_POSTOP_ACEPTACION_MM if MODO_POSTOP
                      else DIST_ACEPTACION_MM)

    # --- 2. geometria: puntos y arbol de cada isla ---
    # Se calcula para TODAS, incluidas las de cortical baja, porque el filtro
    # de material dejo de ser un borrado ciego y ahora necesita la distancia.
    cache = {}
    for f in filas:
        k = f["etiqueta"]
        pts = _puntosDeSuperficie(labels, k, f["caja"], espaciadoZYX,
                                  MAX_PUNTOS_KDTREE)
        cache[k] = (pts, cKDTree(pts))

    kRef = referencia["etiqueta"]
    arbolRef = cache[kRef][1]
    referencia["distRef"] = 0.0
    for f in filas:
        if f["etiqueta"] == kRef:
            continue
        d, _ = arbolRef.query(cache[f["etiqueta"]][0], k=1)
        f["distRef"] = float(d.min())

    centroRef = np.mean(cache[kRef][0], axis=0)
    for f in filas:
        f["env"] = _envolvencia(cache[f["etiqueta"]][0], centroRef)

    # --- 3. primera vuelta de veredictos ---
    limiteExt = FACTOR_EXTENSION_MAXIMA * max(referencia["extMM"])
    aceptadas, pendientes = [], []

    for f in filas:
        if f["referencia"]:
            f["veredicto"] = "ACEPTADA"
            f["categoria"] = "referencia"
            aceptadas.append(f)
            continue

        if max(f["extMM"]) > limiteExt:
            # Un fragmento de craneo NO puede ser mas grande que el craneo del
            # que salio. Limite RELATIVO a la referencia: escala con el
            # paciente, asi que no se rompe con la edad.
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "EXTENSION"
            f["motivo"] = ("mide %.0f mm, mas de %.2fx la referencia (%.0f mm)"
                           % (max(f["extMM"]), FACTOR_EXTENSION_MAXIMA, limiteExt))
            if f["vol"] >= 5.0 and f["fracCort"] >= 0.30:
                avisar("CranioPlan: *** AVISO FUERTE - se elimino por extension "
                       "una isla de %.1f cm3 con %.0f%% de cortical. Si eso era "
                       "hueso, la REFERENCIA esta mal elegida: mira la tabla de "
                       "candidatas y corregi con usar_referencia(n)."
                       % (f["vol"], f["fracCort"] * 100.0))
            continue

        if (not MODO_POSTOP) and f["compac"] < COMPACIDAD_MINIMA_PIEZA:
            # Nube difusa: mucha extension y casi nada de volumen. El ruido de
            # JUAN mide 127 mm de largo con 3.56 cm3, o sea es MAS LARGO que
            # el craneo y pesa dieciseis veces menos. Un fragmento de hueso
            # real, por fino que sea, es un objeto compacto.
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "DIFUSA"
            f["motivo"] = ("difusa: compacidad %.1f < %.1f (mide %.0f mm y solo "
                           "%.2f cm3)" % (f["compac"], COMPACIDAD_MINIMA_PIEZA,
                                          max(f["extMM"]), f["vol"]))
            continue

        if f["esLosa"]:
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "LOSA"
            f["motivo"] = ("losa de soporte: lado corto %.1f mm y %.0f%% del campo"
                           % (min(f["extMM"]), max(f["fracFOV"]) * 100.0))
            continue

        if f["fracCort"] < FRACCION_CORTICAL_MINIMA:
            # CAMBIO CLAVE DE v7.2: el filtro de material ya no es borrado
            # duro. Motivo: en MORENO apareció hueso craneal real con 6.0 %
            # de cortical contra un umbral de 5 %, y el chupete de JUAN tiene
            # 3.88 %. Los rangos se tocan: NO existe umbral que deje el hueso
            # adentro y el chupete afuera. Entonces lo que esta pegado al
            # craneo no se borra, se manda a revision.
            if f["distRef"] is not None and f["distRef"] <= distAceptacion:
                f["veredicto"] = "DUDOSA"
                f["categoria"] = "MATERIAL_PEGADO"
                f["motivo"] = ("cortical %.1f%% (baja) pero a solo %.1f mm de la "
                               "referencia" % (f["fracCort"] * 100.0, f["distRef"]))
            else:
                f["veredicto"] = "ELIMINADA"
                f["categoria"] = "MATERIAL"
                f["motivo"] = ("cortical %.1f%% < %.0f%% y a %.0f mm"
                               % (f["fracCort"] * 100.0,
                                  FRACCION_CORTICAL_MINIMA * 100.0, f["distRef"]))
            continue

        if f["distRef"] is not None and f["distRef"] <= distAceptacion:
            f["veredicto"] = "ACEPTADA"
            f["categoria"] = "pegada"
            aceptadas.append(f)
        else:
            pendientes.append(f)

    # --- 4. distancia a la masa aceptada, para las pendientes ---
    for f in filas:
        if f in aceptadas:
            f["distMasa"] = 0.0 if f["referencia"] else f["distRef"]
            continue
        dmin = None
        for a in aceptadas:
            d, _ = cache[a["etiqueta"]][1].query(cache[f["etiqueta"]][0], k=1)
            dm = float(d.min())
            dmin = dm if dmin is None else min(dmin, dm)
        f["distMasa"] = dmin

    dudosas = [f for f in filas if f["veredicto"] == "DUDOSA"]
    for f in pendientes:
        if f["distMasa"] is not None and f["distMasa"] <= DIST_DUDOSA_MM:
            f["veredicto"] = "DUDOSA"
            f["categoria"] = "DISTANCIA"
            f["motivo"] = ("a %.1f mm de la referencia (limite %.1f mm)"
                           % (f["distRef"], distAceptacion))
            dudosas.append(f)
        else:
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "LEJOS"
            f["motivo"] = "a %.0f mm de la masa aceptada" % f["distMasa"]

    aceptadas.sort(key=lambda f: f["vol"], reverse=True)
    dudosas.sort(key=lambda f: f["vol"], reverse=True)

    totalAceptado = sum(f["vol"] for f in aceptadas)
    if totalAceptado < VOLUMEN_CRANEO_MINIMO_CM3:
        avisar("CranioPlan: *** AVISO FUERTE - lo aceptado suma %.1f cm3, menos "
               "que el minimo anatomico de %.1f cm3. La referencia puede estar "
               "mal elegida o la serie no servir."
               % (totalAceptado, VOLUMEN_CRANEO_MINIMO_CM3))
    if referencia["tocaBorde"]:
        avisar("CranioPlan: NOTA DE CALIDAD - la masa principal toca el borde "
               "del campo de vision. El craneo puede venir RECORTADO y el "
               "modelo incompleto. Esto no cambia la segmentacion; es un dato "
               "para el cirujano.")
    return referencia, aceptadas, dudosas, relajado


def tablaDeCandidatas(filas, referencia, cuantas=3):
    """Las mejores candidatas a referencia, para que la decision sea visible
    y no solo su resultado."""
    orden = sorted(filas, key=lambda f: f["vol"], reverse=True)[:max(cuantas, 3)]
    lineas = ["  Candidatas a referencia (por volumen):",
              "    #   vol(cm3)   ext(mm)  cort%   compac   estado"]
    for i, f in enumerate(orden, start=1):
        if referencia is not None and f is referencia:
            estado = "<= ELEGIDA"
        elif not f["admisible"]:
            motivos = []
            if f["fracCort"] < FRACCION_CORTICAL_MINIMA:
                motivos.append("cortical baja")
            if f["compac"] < COMPACIDAD_MINIMA_REFERENCIA:
                motivos.append("alargada")
            estado = "excluida: " + " y ".join(motivos)
        else:
            estado = "admisible"
        lineas.append("   %2d  %9.2f  %8.1f  %5.1f  %7.1f   %s"
                      % (i, f["vol"], max(f["extMM"]), f["fracCort"] * 100.0,
                         f["compac"], estado))
    return "\n".join(lineas)
# ------------------------------------------------------------
# Glue con Slicer
# ------------------------------------------------------------

def _borrarNodosPorNombre(nombre):
    """Con nombres fijos, un duplicado rompe el Corte V9 en silencio: si queda
    un 'Craneo_Automatico' viejo, Slicer bautiza al nuevo
    'Craneo_Automatico_1' y getNode() agarra el VIEJO."""
    nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    while nodo is not None:
        slicer.mrmlScene.RemoveNode(nodo)
        nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)


def _forzarSuperficie(segNode):
    """CreateClosedSurfaceRepresentation() no reconvierte si la representacion
    ya existe: despues de tocar el labelmap puede devolver la malla vieja."""
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
    nodos = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
    if not nodos:
        print("CranioPlan: no hay volumenes cargados.")
        return []
    print("CranioPlan: volumenes en la escena:")
    for i, n in enumerate(nodos):
        dims = n.GetImageData().GetDimensions() if n.GetImageData() else (0, 0, 0)
        esp = n.GetSpacing()
        print("  [%d] %-45s  %dx%dx%d  %.3f/%.3f/%.3f mm"
              % (i, n.GetName(), dims[0], dims[1], dims[2], esp[0], esp[1], esp[2]))
    return nodos


def usar_volumen(referencia=None):
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


def _encabezado(volumeNode):
    dims = volumeNode.GetImageData().GetDimensions()
    esp = volumeNode.GetSpacing()
    print("CranioPlan: volumen : %s" % volumeNode.GetName())
    print("CranioPlan: dims    : %d x %d x %d voxels" % dims)
    print("CranioPlan: spacing : %.3f x %.3f x %.3f mm" % esp)
    print("CranioPlan: FOV     : %.1f x %.1f x %.1f mm"
          % (dims[0] * esp[0], dims[1] * esp[1], dims[2] * esp[2]))
    n = volumeNode.GetName().lower()
    if ("cerebro" in n or "brain" in n) and "hueso" not in n:
        print("CranioPlan: AVISO - serie de kernel BLANDO. El hueso delgado "
              "pierde HU. Si el estudio tiene serie de HUESO, preferila.")


def _limpiarEscenaPrevia():
    for clave in ("segNode", "modelo"):
        nodo = _estado.get(clave)
        if nodo is not None:
            try:
                slicer.mrmlScene.RemoveNode(nodo)
            except Exception:
                pass
        _estado[clave] = None
    _borrarNodosPorNombre(NOMBRE_NODO_SEGMENTACION)
    _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)
    _estado["colores"] = {}
    _estado["candidatas"] = []


def postop(activar=True):
    """Activa el modo postoperatorio y vuelve a decidir sin re-segmentar.

    En un craneo ya operado los colgajos estan genuinamente separados
    (MORENO 501: 2.2 a 16.7 mm de la masa principal). Subir la distancia de
    aceptacion los recupera, a costa de que entren tambien las vertebras: en
    este modo la revision manual pasa a ser OBLIGATORIA, no opcional."""
    global MODO_POSTOP
    MODO_POSTOP = bool(activar)
    print("CranioPlan: MODO POSTOPERATORIO %s (distancia de aceptacion %.0f mm)."
          % ("ACTIVADO" if MODO_POSTOP else "desactivado",
             DIST_POSTOP_ACEPTACION_MM if MODO_POSTOP else DIST_ACEPTACION_MM))
    if MODO_POSTOP:
        print("CranioPlan: en este modo van a entrar vertebras y fragmentos que "
              "en preoperatorio quedarian afuera. Revisa TODO antes de "
              "confirmar_craneo().")
    if _estado.get("filas"):
        return _redecidir()
    return None


def usar_referencia(numero):
    """Fuerza la referencia a la isla `numero` de la tabla y vuelve a decidir,
    sin re-segmentar el volumen (que es la parte lenta)."""
    if not _estado.get("filas"):
        print("CranioPlan: primero corre generar().")
        return None
    return _redecidir(referenciaForzada=int(numero))


def _redecidir(referenciaForzada=None):
    """Rehace la decision y reconstruye la segmentacion con las filas ya
    analizadas."""
    labels = _estado.get("labels")
    filas = _estado.get("filas")
    volumeNode = _estado.get("volumeNode")
    if labels is None or not filas or volumeNode is None:
        print("CranioPlan: no hay un analisis previo en memoria.")
        return None
    referencia, aceptadas, dudosas, relajado = decidir(
        labels, filas, _estado["espaciadoZYX"], referenciaForzada=referenciaForzada)
    if referencia is None:
        return None
    print("")
    print(tablaDeCandidatas(filas, referencia))
    _imprimirTabla(filas)
    return _construirEscena(labels, volumeNode, referencia, aceptadas, dudosas)


def generar():
    """BLOQUE A v7.2. Segmenta, decide y deja las piezas listas para la
    revision. No fusiona nada: eso lo hace confirmar_craneo()."""
    print("=" * 66)
    print("CranioPlan - Bloque A %s" % CRANIOPLAN_BLOQUE_A)
    if MODO_POSTOP:
        print("CranioPlan: *** MODO POSTOPERATORIO ACTIVO ***")

    volumeNode = _estado.get("volumeNode")
    if volumeNode is None or slicer.mrmlScene.GetNodeByID(volumeNode.GetID()) is None:
        volumeNode = usar_volumen(NOMBRE_VOLUMEN)
    if volumeNode is None:
        return None
    _encabezado(volumeNode)
    _limpiarEscenaPrevia()

    volArr = slicer.util.arrayFromVolume(volumeNode)
    esp = volumeNode.GetSpacing()
    espaciadoZYX = (esp[2], esp[1], esp[0])

    labels, filas, info = analizarVolumen(volArr, espaciadoZYX, HU_MIN, HU_MAX)
    if labels is None or not filas:
        print("CranioPlan: no se encontro ninguna isla util en %d-%d HU."
              % (HU_MIN, HU_MAX))
        return None

    print("CranioPlan: islas detectadas: %d  (%d por encima de %.2f cm3)"
          % (info["nIslas"], len(filas), VOL_MINIMO_CM3))
    print("CranioPlan: total en mascara: %.2f cm3" % info["totalCM3"])
    if info["valorRelleno"] is not None:
        print("CranioPlan: relleno fuera del FOV en %d HU (borde circular "
              "detectado)." % info["valorRelleno"])

    _estado["labels"] = labels
    _estado["espaciadoZYX"] = espaciadoZYX
    _estado["filas"] = filas

    referencia, aceptadas, dudosas, relajado = decidir(labels, filas, espaciadoZYX)
    if referencia is None:
        print("CranioPlan: ninguna isla puede ser referencia. Revisa los "
              "umbrales HU o el volumen elegido.")
        return None

    print("")
    print(tablaDeCandidatas(filas, referencia))
    print("CranioPlan: referencia -> %.2f cm3, cortical %.1f%%, compacidad %.1f, "
          "ext %.0f mm" % (referencia["vol"], referencia["fracCort"] * 100.0,
                           referencia["compac"], max(referencia["extMM"])))
    if relajado:
        print("CranioPlan: OJO - se relajo el filtro de %s para poder elegir "
              "referencia. Mirala con desconfianza." % relajado)

    _imprimirTabla(filas)
    return _construirEscena(labels, volumeNode, referencia, aceptadas, dudosas)


def _construirEscena(labels, volumeNode, referencia, aceptadas, dudosas):
    _limpiarEscenaPrevia()
    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_SEGMENTACION)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentacion = segNode.GetSegmentation()

    buf = np.zeros(labels.shape, dtype=np.uint8)
    lista = []
    for i, f in enumerate(aceptadas + dudosas, start=1):
        etiquetaZona = "OK" if f["veredicto"] == "ACEPTADA" else "DUDOSA"
        segId = segmentacion.AddEmptySegment(
            "", "%02d_%s_%.2fcm3" % (i, etiquetaZona, f["vol"]))
        caja = f["caja"]
        buf[caja] = (labels[caja] == f["etiqueta"]).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(buf, segNode, segId, volumeNode)
        buf[caja] = 0
        color = COLOR_ACEPTADA if f["veredicto"] == "ACEPTADA" else COLOR_DUDOSA
        segmentacion.GetSegment(segId).SetColor(*color)
        _estado["colores"][segId] = color
        lista.append({
            "numero": i, "segId": segId, "vol": f["vol"],
            "distRef": f["distRef"], "fracCort": f["fracCort"],
            "compac": f["compac"], "env": f["env"], "zona": f["veredicto"],
            "categoria": f["categoria"], "nota": f["motivo"],
            "referencia": f["referencia"],
        })
    del buf
    _forzarSuperficie(segNode)

    _estado["segNode"] = segNode
    _estado["candidatas"] = lista
    _estado["referencia"] = referencia

    print("")
    print("--- PIEZAS EN LA ESCENA ---")
    for d in lista:
        marca = " *referencia*" if d["referencia"] else ""
        extra = ("   <- " + d["nota"]) if d["nota"] else ""
        print("  [%2d] %-8s vol=%8.2f cm3  cort=%5.1f%%  dRef=%5.1f mm%s%s"
              % (d["numero"], d["zona"], d["vol"], d["fracCort"] * 100.0,
                 d["distRef"] if d["distRef"] is not None else -1, marca, extra))
    nd = sum(1 for d in lista if d["zona"] == "DUDOSA")
    print("-" * 60)
    print("CranioPlan: %d aceptada(s) + %d dudosa(s) (naranja en el visor)."
          % (len(lista) - nd, nd))
    if nd:
        print("CranioPlan: las dudosas se CONSERVAN a proposito. Miralas con "
              "resaltar(n) y borra las que no sirvan con eliminar(n), o todas "
              "juntas con eliminar_dudosas().")
    ayuda()
    return lista


def _imprimirTabla(filas):
    print("")
    print("--- TODAS LAS ISLAS ANALIZADAS ---")
    print("   vol(cm3)  p95(HU)  cort%  compac   env%   dRef  dMasa  espMin  ext"
          "   %FOV  borde  veredicto  motivo")
    for f in filas:
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "    -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "    -" if f["distMasa"] is None else "%5.1f" % f["distMasa"]
        marca = "  *ref*" if f["referencia"] else ""
        print("  %9.2f  %7.0f  %5.1f  %6.1f  %s  %s  %s  %6.1f %5.0f %4.0f%%"
              "  %-5s  %-9s %s%s"
              % (f["vol"], f["p95"], f["fracCort"] * 100.0, f["compac"], env,
                 dr, dm, min(f["extMM"]), max(f["extMM"]),
                 max(f["fracFOV"]) * 100, "si" if f["tocaBorde"] else "no",
                 f["veredicto"] or "?", f["motivo"], marca))
    print("-" * 60)


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
    segNode = _estado.get("segNode")
    if segNode is None:
        return
    disp = segNode.GetDisplayNode()
    segmentacion = segNode.GetSegmentation()
    for d in _estado["candidatas"]:
        disp.SetSegmentVisibility(d["segId"], True)
        seg = segmentacion.GetSegment(d["segId"])
        if seg is not None:
            seg.SetColor(*_estado["colores"].get(d["segId"], (0.6, 0.6, 0.6)))


def solo_dudosas():
    """Deja a la vista unicamente lo que el algoritmo no supo decidir."""
    segNode = _estado.get("segNode")
    if segNode is None:
        print("CranioPlan: primero corre generar().")
        return
    disp = segNode.GetDisplayNode()
    n = 0
    for d in _estado["candidatas"]:
        v = d["zona"] == "DUDOSA"
        disp.SetSegmentVisibility(d["segId"], v)
        n += 1 if v else 0
    print("CranioPlan: %d dudosa(s) a la vista." % n)


def eliminar(numero):
    segNode = _estado.get("segNode")
    if segNode is None:
        print("CranioPlan: primero corre generar().")
        return
    d = next((x for x in _estado["candidatas"] if x["numero"] == numero), None)
    if d is None:
        print("CranioPlan: no existe la pieza %s." % numero)
        return
    if d["referencia"]:
        print("CranioPlan: esa es la masa principal del craneo. Si de verdad "
              "queres cambiarla, usa usar_referencia(n) en vez de borrarla.")
        return
    segNode.GetSegmentation().RemoveSegment(d["segId"])
    _estado["candidatas"] = [x for x in _estado["candidatas"] if x["numero"] != numero]
    print("CranioPlan: pieza %d eliminada. Quedan %d."
          % (numero, len(_estado["candidatas"])))


def eliminar_dudosas():
    """Borra de una sola vez todas las marcadas como dudosas. En los 9 casos
    relevados casi todas eran vertebras cervicales."""
    dudosas = [d["numero"] for d in _estado.get("candidatas", [])
               if d["zona"] == "DUDOSA"]
    if not dudosas:
        print("CranioPlan: no hay dudosas.")
        return
    for n in sorted(dudosas, reverse=True):
        eliminar(n)
    print("CranioPlan: %d dudosa(s) eliminadas." % len(dudosas))
def confirmar_craneo():
    """Fusiona las piezas que quedaron en un unico segmento con el nombre que
    espera el Corte V9.

    Se hace con un OR de mascaras en numpy y UNA sola escritura, no con el
    efecto 'Logical operators' del Segment Editor: ese efecto arrastra las
    reglas de masking del editor (estado oculto) y dejaba viva la
    representacion de superficie vieja de cada segmento."""
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
    segmentacion.GetSegment(base).SetColor(*COLOR_ACEPTADA)

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
    print("CranioPlan: listo. Ahora: enviar_a_planner()")
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
    """Umbral ABSOLUTO en cantidad de puntos. Con un umbral RELATIVO,
    cualquier placa craneal desconectada por una sutura abierta que mida
    menos del 10 % de la boveda se borra sola: por eso 'Craneo_Final'
    aparecia sin mandibula en el 403."""
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
    if descartadas and descartadas[0] >= minimoPuntos * 0.5:
        print("CranioPlan:   ATENCION - la mayor descartada esta cerca del "
              "umbral. Si falta hueso, bajar MINIMO_PUNTOS_MALLA.")
    return salida


def _decimarMalla(polyData, reduccion):
    """Decima verificando que no desaparezcan regiones conexas: preferimos
    render pesado antes que perder hueso."""
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
    """Exporta el segmento a un Model node con el nombre que espera el
    Corte V9."""
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
    _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)
    modelo = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLModelNode', NOMBRE_MODELO_CRANEO)
    modelo.SetAndObservePolyData(malla)
    modelo.CreateDefaultDisplayNodes()
    modelo.GetDisplayNode().SetColor(*COLOR_ACEPTADA)
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
        print("CranioPlan: ATENCION - el modelo quedo como '%s'. Habia otro "
              "nodo con ese nombre: el Corte V9 va a usar el equivocado."
              % modelo.GetName())
    else:
        print("CranioPlan: modelo '%s' creado." % modelo.GetName())
    print("CranioPlan: ya podes dibujar las curvas y correr cortar().")
    return modelo


def verificar_compatibilidad_v9():
    """Chequea que la escena tenga EXACTAMENTE lo que el Corte V9 va a
    buscar. Correrlo antes de cortar() ahorra el 'no encuentro el nodo...'."""
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
        print("  [FALTA] nodo de segmentacion: corre generar()")
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
        print("  [FALTA] Model node: corre enviar_a_planner()")
        ok = False
    else:
        print("  [ERROR] hay %d modelos llamados '%s'."
              % (len(mods), NOMBRE_MODELO_CRANEO))
        ok = False

    nCer = len(slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode'))
    nAbi = len([n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')])
    print("  curvas: %d cerrada(s), %d abierta(s)" % (nCer, nAbi))
    if nCer + nAbi == 0:
        print("  [FALTA] dibuja al menos una curva de corte.")
        ok = False
    print("  -> %s" % ("todo listo para cortar()" if ok else "falta algo"))
    print("--------------------------------------")
    return ok




def ayuda():
    print("")
    print("Comandos disponibles:")
    print("  listar_volumenes() / usar_volumen(i)")
    print("  generar()                   -> corre el Bloque A v7.2")
    print("  usar_referencia(n)          -> forzar la referencia a mano")
    print("  postop()                    -> modo postoperatorio (craneo ya cortado)")
    print("  solo_dudosas()              -> ver solo lo no decidido")
    print("  resaltar(n) / mostrar_todas()")
    print("  eliminar(n) / eliminar_dudosas()")
    print("  confirmar_craneo()          -> fusiona en '%s'" % NOMBRE_SEGMENTO_CRANEO)
    print("  enviar_a_planner()          -> exporta el modelo '%s'" % NOMBRE_MODELO_CRANEO)
    print("  verificar_compatibilidad_v9() -> chequea antes de cortar()")


generar()
