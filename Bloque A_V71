# ============================================================
# CRANIOPLAN - BLOQUE A v7.1  (script de consola para 3D Slicer)
#            COMPATIBLE CON EL CORTE V9
# ============================================================
# Reemplaza al Bloque A v6. Se pega entero en la consola de Python de
# Slicer. Al pegarlo se ejecuta solo.
#
# ============================================================
# DE DONDE SALEN LOS UMBRALES DE ESTE ARCHIVO
# ============================================================
# De 143 islas etiquetadas a mano sobre 6 series de 5 pacientes del
# Hospital Garrahan (JUAN 702, RUIZ 175, AVILA 403, LUNA 167,
# JUSTINA 559/381 en sus dos reconstrucciones). El relevamiento se hizo
# con el script BloqueA_v6_DIAGNOSTICO.py, que lista sin eliminar.
#
# Sobre las 80 islas con referencia valida (se excluye RUIZ, ver [C1]):
#
#            distancia a la masa      fraccion cortical
#   CRANEO    0.00 - 1.70 mm          8.6 % - 59.0 %     (n=17)
#   BASURA    0.96 - 53.6 mm          0.0 % - 57.2 %     (n=63)
#
# CONCLUSION 1: la fraccion cortical NO separa hueso de basura. La
# camilla de RUIZ tiene 69.8 % de cortical y la de JUAN 57 %; la base de
# craneo de AVILA tiene 8.6 %. Los rangos se superponen enteros. En ese
# equipo el soporte de inmovilizacion es tan denso como hueso cortical.
# Sirve solo para descartar espuma, plastico y ruido (cortical ~ 0 %).
#
# CONCLUSION 2: la DISTANCIA separa casi perfecto. TODO el hueso craneal
# real -calota, base, mandibula, en los 6 estudios- esta a menos de
# 2 mm de la masa principal. Sin excepciones.
#
# REGLA RESULTANTE:  cortical >= 5 %  Y  distancia <= 2 mm
#   -> 0 falsos negativos y 0 falsos positivos sobre las 80 islas.
#   Margenes: el hueso mas marginal tiene 8.6 % de cortical contra un
#   umbral de 5; la basura mas cercana con cortical alta (la almohadilla
#   de JUAN) esta a 2.51 mm contra un umbral de 2.
#
# ------------------------------------------------------------
# QUE CAMBIA RESPECTO DEL v6
# ------------------------------------------------------------
# [1] EL VETO DE BORDE YA NO ELIMINA NADA. Era la causa dominante de
#     falsos negativos: se comio la calota entera de RUIZ (103 cm3), las
#     dos mandibulas de LUNA, las dos de JUSTINA 584 y la de JUSTINA 601
#     (16 cm3). Ahora solo se usa para ELEGIR LA REFERENCIA en la pasada
#     1, que es donde si funciona (resuelve el caso de JUAN).
#
# [2] EL PISO DE VOLUMEN BAJA DE 0.5 A 0.05 cm3. El piso de 0.5 se comia
#     la placa craneal de 0.09 cm3 de JUAN (cortical 19 %, a 1.3 mm).
#     La decision la toma ahora material + distancia, no el tamano.
#
# [3] EL MARGEN DE PROXIMIDAD BAJA DE 40 A 2 mm. 40 mm era el numero que
#     dejaba entrar almohadillas densas (JUSTINA 584: dos islas de
#     almohadilla aceptadas a 13 y 16 mm). 2 mm es lo que los datos
#     muestran como rango real de las suturas abiertas.
#
# [4] APARECE LA ZONA DUDOSA. Nada que sea plausiblemente hueso se
#     elimina solo. Entre 2 y 15 mm la isla se conserva, se pinta
#     distinto y se lista para que la borre el cirujano. Son 3 a 11
#     islas por caso, casi todas vertebras cervicales.
#
# [5] EL FILTRO DE DENSIDAD DEJA DE USAR p95 Y USA FRACCION CORTICAL.
#     El p95 de una lamina de 1-2 mm esta arruinado por volumen parcial;
#     la cola superior de HU no. Ademas el umbral p95 < 500 decidia por
#     margenes de 20 HU (JUSTINA: 464, 481, 475 contra 500), que es del
#     orden del ruido.
#
# ------------------------------------------------------------
# QUE CAMBIA EN v7.1 RESPECTO DE v7
# ------------------------------------------------------------
# [6] SE ELIMINA EL CRECIMIENTO TRANSITIVO. El v7 aceptaba una isla si
#     estaba a <= 2 mm de CUALQUIER isla ya aceptada, asi que una isla
#     aceptada servia de trampolin para la siguiente. Consecuencia
#     observada en la TC de JUAN: hay una vertebra a 1.57 mm del craneo;
#     entra, y desde ella entra la siguiente, y la columna cervical sube
#     entera como ACEPTADA. En LUNA, AVILA y JUSTINA la primera vertebra
#     queda a mas de 2 mm, la cadena no arranca, y todas quedan DUDOSAS.
#     El mismo algoritmo daba resultados distintos segun la flexion del
#     cuello del paciente, que es justo el tipo de dependencia que este
#     rediseno vino a sacar.
#
#     Ademas, el crecimiento transitivo NO es la regla que se valido: la
#     columna 'dist' del CSV mide distancia a la isla de REFERENCIA, y es
#     con esa definicion que la separacion dio perfecta (todo el hueso
#     craneal a <= 1.70 mm de la referencia, directo, sin intermediarios).
#     v7.1 acepta unicamente lo que esta a <= DIST_ACEPTACION_MM de la
#     REFERENCIA. Lo demas que parezca hueso y este dentro de
#     DIST_DUDOSA_MM de la masa aceptada queda DUDOSO.
#
#     Riesgo asumido: una calota partida en tres, donde la tercera placa
#     toca a la segunda pero esta lejos de la referencia, quedaria DUDOSA
#     en vez de aceptada. En los 6 estudios no ocurrio nunca. Y si
#     ocurriera, el costo es un click, no hueso perdido.
#
# [7] SE MIDE LA ENVOLVENCIA, PERO NO DECIDE NADA TODAVIA. Es la fraccion
#     de direcciones del espacio, vistas desde el centro del craneo, en
#     las que aparece la isla. Mide "esta cosa me rodea?", que es lo que
#     define fisicamente a una almohadilla de inmovilizacion y no depende
#     ni del encuadre del FOV, ni de la orientacion del paciente, ni del
#     kernel de reconstruccion.
#
#     Se imprime en la columna env% de la tabla. NO se usa como criterio
#     de descarte porque hay 3 almohadillas en todo el dataset y ningun
#     valor medido de esta metrica: poner un umbral ahora seria repetir
#     el error del p95 < 500 HU y del 5% relativo. Cuando haya 6 u 8
#     almohadillas medidas, el umbral sale de los datos.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (documentados a proposito, no ignorados)
# ------------------------------------------------------------
# [C1] LAS DISTANCIAS DE RUIZ NO SE USARON PARA CALIBRAR. En ese caso el
#      veto de borde del v6 elimino la calota y la referencia termino
#      siendo una vertebra de 2.38 cm3, asi que las distancias de esa
#      tabla estan medidas contra la vertebra y no contra el craneo (las
#      mandibulas figuran a 19-20 mm). Es un artefacto del bug, no un
#      dato. Con la referencia corregida RUIZ entra en el mismo patron,
#      pero conviene volver a correr el diagnostico sobre ese caso para
#      confirmarlo con numeros en vez de con un argumento.
#
# [C2] EL CHUPETE NO ES SEPARABLE POR MATERIAL. El chupete de JUAN tiene
#      7.59 cm3, cortical 3.88 %, a 1.69 mm del craneo. El hueso mas
#      marginal tiene 8.6 %. El umbral de 5 % lo separa, pero el margen
#      real es de 1.1 puntos porcentuales de un lado y 4.7 del otro. Con
#      UN solo ejemplo de chupete no se puede afirmar que sea estable.
#      Si aparece un chupete que pase el filtro, va a quedar aceptado y
#      hay que borrarlo a mano.
#
# [C3] LAS VERTEBRAS SON HUESO REAL. Ningun criterio de material o de
#      distancia las distingue de la mandibula: estan pegadas a la base
#      del craneo por los condilos occipitales. Entran como aceptadas o
#      como dudosas segun el caso. Se borran con eliminar(n) o de una
#      sola vez con eliminar_dudosas(). NO se eliminan solas a proposito:
#      seria eliminar hueso por un criterio anatomico que no tenemos.
#
# [C4] LA ALMOHADILLA CON CORTICAL ALTA Y PEGADA AL CRANEO PASARIA.
#      La de JUAN (79 cm3, cortical 56 %) esta a 2.51 mm. Con el umbral
#      en 2.0 queda como dudosa, no como aceptada, pero el margen es de
#      medio milimetro. Es el punto mas fragil de toda la calibracion.
#
# [C5] MAX_EXTENSION_CRANEO_MM = 180 es una cota anatomica para lactantes
#      y ninos pequenos, que son la poblacion de esta tesis. Si alguna
#      vez se usa en un adolescente hay que subirla. Solo afecta a quien
#      puede ser REFERENCIA, nunca elimina nada.
#
# [C6] SEIS SERIES DE CINCO PACIENTES, DE UN SOLO HOSPITAL. Los umbrales
#      son defendibles sobre estos datos y nada mas. Cada caso nuevo hay
#      que pasarlo primero por el script de diagnostico y sumarlo al CSV.
# ============================================================


# ============================================================
# CONFIG
# ============================================================

# --- Nombres de nodos (deben coincidir con la config del Corte V9) ---
NOMBRE_NODO_SEGMENTACION = "Craneo_Automatico"   # V9: NOMBRE_SEGMENTACION
NOMBRE_SEGMENTO_CRANEO   = "Craneo_Final"        # V9: NOMBRE_SEGMENTO_HUESO
NOMBRE_MODELO_CRANEO     = "Craneo_Final"        # V9: NOMBRE_MODELO_CRANEO

# Volumen a usar. None = el ultimo volumen escalar cargado en la escena.
NOMBRE_VOLUMEN = None

# Umbrales de hueso (Hounsfield).
HU_MIN = 300
HU_MAX = 3000

# --- Filtro de material (ver CONCLUSION 1) ---
# Fraccion de voxeles de la isla por encima de UMBRAL_CORTICAL_HU.
# Descarta espuma, plastico, gel y ruido. NO descarta camillas densas.
UMBRAL_CORTICAL_HU = 700.0
FRACCION_CORTICAL_MINIMA = 0.05     # 5 %  (hueso mas marginal medido: 8.6 %)

# --- Distancias (ver CONCLUSION 2) ---
# OJO con la semantica, cambio en v7.1 (ver [6]):
#   DIST_ACEPTACION_MM se mide contra la isla de REFERENCIA.
#   DIST_DUDOSA_MM se mide contra la masa ya aceptada.
DIST_ACEPTACION_MM = 2.0    # hueso craneal real medido: maximo 1.70 mm
DIST_DUDOSA_MM     = 15.0   # entre 2 y 15 -> se conserva y se marca

# Envolvencia: resolucion de la esfera de direcciones (ver [7]).
# Binning de area igual: uniforme en cos(theta) y en phi.
ENV_BINS_COSTHETA = 12
ENV_BINS_PHI      = 24

# Piso ABSOLUTO de volumen para que una isla se analice.
VOL_MINIMO_CM3 = 0.05

# Volumen oseo TOTAL minimo para considerar que lo aceptado es un craneo.
# Solo decide si hace falta el rescate de la pasada 2.
VOLUMEN_CRANEO_MINIMO_CM3 = 20.0

# Extension maxima de un craneo pediatrico. Solo limita quien puede ser
# REFERENCIA; nunca elimina. Impide que la almohadilla envolvente de JUAN
# (195.7 mm de extension) sea elegida como craneo.
MAX_EXTENSION_CRANEO_MM = 180.0

# Detector de losa de soporte (tabla de la camilla). Nunca agarro hueso
# real en los 6 casos.
ESPESOR_LOSA_MM = 30.0
FRACCION_LOSA_FOV = 0.85

# Exportacion de malla.
MINIMO_PUNTOS_MALLA = 200      # umbral ABSOLUTO de ruido de marching cubes
REDUCCION_MALLA = 0.7          # 0.0 desactiva la decimacion
SUAVIZADO_SUPERFICIE = "0.3"   # bajar a "0.1" si se disuelven placas finas

# Rendimiento.
MAX_PUNTOS_KDTREE = 40000


# ============================================================
# CODIGO
# ============================================================

import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

CRANIOPLAN_BLOQUE_A = "v7.1 (2026-08-07) consola - calibrado sobre 6 series"

_estado = {
    "volumeNode": None,
    "labels": None,
    "espaciadoZYX": None,
    "filas": [],
    "candidatas": [],
    "segNode": None,
    "colores": {},
    "modelo": None,
}

COLOR_ACEPTADA = (0.90, 0.80, 0.60)
COLOR_DUDOSA   = (1.00, 0.60, 0.10)


# ------------------------------------------------------------
# Utilidades puras (numpy)
# ------------------------------------------------------------

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
    dilatacion pediria otro array de volumen entero."""
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
    delgada en su eje mas fino."""
    if max(extMM) <= 0.0:
        return False
    return bool(min(extMM) <= ESPESOR_LOSA_MM and max(fracFOV) >= FRACCION_LOSA_FOV)


def _puntosDeSuperficie(labels, etiqueta, caja, espaciadoZYX, maxPuntos):
    """Coordenadas fisicas (mm) de los voxeles de SUPERFICIE de una isla.
    La distancia minima entre dos solidos se alcanza siempre en sus
    superficies."""
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


def analizarVolumen(volArr, espaciadoZYX):
    """Etiqueta las islas y calcula metricas. No decide nada."""
    mascara = (volArr >= HU_MIN) & (volArr <= HU_MAX)
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
        p95 = float(np.percentile(valores, 95))
        fracCort = float(np.count_nonzero(valores >= UMBRAL_CORTICAL_HU)) / float(valores.size)
        filas.append({
            "etiqueta": k, "vol": vol, "caja": caja,
            "extMM": extMM, "fracFOV": fracFOV,
            "p95": p95, "fracCort": fracCort,
            "tocaBorde": k in etiquetasBorde,
            "esLosa": _esLosaDeSoporte(extMM, fracFOV),
            "distRef": None,   # a la isla de REFERENCIA  -> decide ACEPTADA
            "distMasa": None,  # a la masa ya aceptada    -> decide DUDOSA
            "env": None,       # envolvencia, solo se mide (ver [7])
            "zona": None, "motivo": "", "referencia": False,
        })

    filas.sort(key=lambda f: f["vol"], reverse=True)
    return labels, filas, {"nIslas": nIslas, "valorRelleno": valorRelleno,
                           "totalCM3": totalCM3}


def _esMaterialOseo(f):
    """Filtro de material. Descarta espuma, plastico, gel y ruido; NO
    descarta camillas densas (para eso esta la distancia)."""
    return f["fracCort"] >= FRACCION_CORTICAL_MINIMA


def _puedeSerReferencia(f, exentarDelBorde):
    if not _esMaterialOseo(f):
        return False
    if f["esLosa"]:
        return False
    if max(f["extMM"]) > MAX_EXTENSION_CRANEO_MM:
        return False        # almohadilla envolvente (ver [C5])
    if not exentarDelBorde and f["tocaBorde"]:
        return False
    return True


def _envolvencia(ptsMM, centroMM):
    """
    Fraccion de direcciones del espacio, vistas desde el centro del craneo,
    en las que aparece esta isla.

    Mide "esta cosa me rodea?". Una almohadilla que envuelve la cabeza
    aparece en una porcion grande del cielo; una mandibula solo en un cono
    hacia adelante y abajo; una vertebra en un cono muy angosto; la camilla
    en una franja por debajo.

    Las direcciones se agrupan en bins de AREA IGUAL (uniformes en
    cos(theta) y en phi), no en una grilla lat-lon, que concentraria bins
    diminutos en los polos y sesgaria la cuenta.

    ATENCION: en v7.1 esto NO decide nada. Se imprime para juntar datos.
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


def _prepararGeometria(oseas, cache, labels, espaciadoZYX):
    """Puntos de superficie y arbol KD de cada isla, una sola vez."""
    for f in oseas:
        k = f["etiqueta"]
        if k not in cache:
            pts = _puntosDeSuperficie(labels, k, f["caja"], espaciadoZYX,
                                      MAX_PUNTOS_KDTREE)
            cache[k] = (pts, cKDTree(pts))
    return cache


def _matrizDeDistancias(oseas, cache):
    """Distancia minima superficie-superficie entre cada par de islas."""
    etiquetas = [f["etiqueta"] for f in oseas]
    D = {a: {a: 0.0} for a in etiquetas}
    for i, a in enumerate(etiquetas):
        for b in etiquetas[i + 1:]:
            d, _ = cache[b][1].query(cache[a][0], k=1)
            dmin = float(d.min())
            D[a][b] = dmin
            D[b][a] = dmin
    return D


def _aceptarDesdeReferencia(referencia, oseas, D):
    """
    v7.1: se acepta UNICAMENTE lo que esta a <= DIST_ACEPTACION_MM de la
    REFERENCIA. Sin encadenamiento (ver [6]).

    El v7 usaba crecimiento transitivo y eso hacia que una vertebra
    aceptada sirviera de trampolin para la siguiente: en JUAN subia la
    columna cervical entera, en los demas casos no. Mismo algoritmo,
    resultados distintos segun la flexion del cuello.
    """
    kRef = referencia["etiqueta"]
    aceptadas = {kRef}
    referencia["distRef"] = 0.0
    referencia["distMasa"] = 0.0
    for f in oseas:
        k = f["etiqueta"]
        if k == kRef:
            continue
        f["distRef"] = D[k][kRef]
        if f["distRef"] <= DIST_ACEPTACION_MM:
            aceptadas.add(k)
    # distancia a la masa aceptada, para ordenar y clasificar el resto
    for f in oseas:
        k = f["etiqueta"]
        f["distMasa"] = min(D[k][a] for a in aceptadas)
    return aceptadas


def decidir(labels, filas, espaciadoZYX):
    """Asigna zona a cada isla. Devuelve (referencia, aceptadas, dudosas)."""
    for f in filas:
        f["distRef"] = None
        f["distMasa"] = None
        f["env"] = None
        f["zona"] = None
        f["motivo"] = ""
        f["referencia"] = False

    # --- 1. filtro de material ---
    oseas = []
    for f in filas:
        if not _esMaterialOseo(f):
            f["zona"] = "ELIMINADA"
            f["motivo"] = ("material no oseo (cortical %.1f%% < %.0f%%)"
                           % (f["fracCort"] * 100.0, FRACCION_CORTICAL_MINIMA * 100.0))
        elif f["esLosa"]:
            f["zona"] = "ELIMINADA"
            f["motivo"] = ("losa de soporte (lado corto %.1f mm, %.0f%% del campo)"
                           % (min(f["extMM"]), max(f["fracFOV"]) * 100.0))
        else:
            oseas.append(f)

    if not oseas:
        return None, [], []

    cache = _prepararGeometria(oseas, {}, labels, espaciadoZYX)
    D = _matrizDeDistancias(oseas, cache)

    # --- 2. referencia, en dos pasadas ---
    referencia = None
    aceptadas = set()
    total = 0.0

    plausibles = [f for f in oseas if _puedeSerReferencia(f, False)]
    if plausibles:
        referencia = max(plausibles, key=lambda f: f["vol"])
        aceptadas = _aceptarDesdeReferencia(referencia, oseas, D)
        total = sum(f["vol"] for f in oseas if f["etiqueta"] in aceptadas)
        print("CranioPlan: pasada 1 (referencia sin tocar el borde): "
              "referencia %.2f cm3, %.1f cm3 aceptados." % (referencia["vol"], total))
    else:
        print("CranioPlan: pasada 1 sin candidata a referencia.")

    if referencia is None or total < VOLUMEN_CRANEO_MINIMO_CM3:
        print("CranioPlan: eso no alcanza para ser un craneo (minimo %.1f cm3). "
              "Modo RESCATE: la referencia queda exenta del borde, porque lo mas "
              "probable es que el FOV este recortando el craneo (caso TC 175)."
              % VOLUMEN_CRANEO_MINIMO_CM3)
        plausibles2 = [f for f in oseas if _puedeSerReferencia(f, True)]
        if plausibles2:
            ref2 = max(plausibles2, key=lambda f: f["vol"])
            if referencia is not None:
                referencia["referencia"] = False
            ac2 = _aceptarDesdeReferencia(ref2, oseas, D)
            total2 = sum(f["vol"] for f in oseas if f["etiqueta"] in ac2)
            if total2 > total:
                referencia, aceptadas, total = ref2, ac2, total2
                print("CranioPlan: rescate aplicado -> referencia %.2f cm3, "
                      "%.1f cm3 aceptados." % (referencia["vol"], total))
            else:
                print("CranioPlan: el rescate no mejoro nada; se vuelve a la pasada 1.")
                aceptadas = _aceptarDesdeReferencia(referencia, oseas, D)

    if referencia is None:
        return None, [], []
    referencia["referencia"] = True

    # --- 3. envolvencia (solo se mide, ver [7]) ---
    centroMM = np.mean(cache[referencia["etiqueta"]][0], axis=0)
    for f in oseas:
        f["env"] = _envolvencia(cache[f["etiqueta"]][0], centroMM)

    # --- 4. zonas ---
    listaAceptadas, listaDudosas = [], []
    for f in oseas:
        if f["etiqueta"] in aceptadas:
            f["zona"] = "ACEPTADA"
            listaAceptadas.append(f)
        elif f["distMasa"] is not None and f["distMasa"] <= DIST_DUDOSA_MM:
            f["zona"] = "DUDOSA"
            f["motivo"] = ("a %.1f mm de la referencia (limite %.1f mm)"
                           % (f["distRef"], DIST_ACEPTACION_MM))
            listaDudosas.append(f)
        else:
            f["zona"] = "ELIMINADA"
            f["motivo"] = ("lejos de la masa aceptada (%.0f mm)" % f["distMasa"])

    listaAceptadas.sort(key=lambda f: f["vol"], reverse=True)
    listaDudosas.sort(key=lambda f: f["vol"], reverse=True)
    return referencia, listaAceptadas, listaDudosas


# ------------------------------------------------------------
# Glue con Slicer
# ------------------------------------------------------------

def _borrarNodosPorNombre(nombre):
    """Con nombres fijos, un duplicado rompe el Corte V9 en silencio: si
    queda un 'Craneo_Automatico' viejo, Slicer bautiza al nuevo
    'Craneo_Automatico_1' y getNode() agarra el VIEJO."""
    nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    while nodo is not None:
        slicer.mrmlScene.RemoveNode(nodo)
        nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)


def _forzarSuperficie(segNode):
    """CreateClosedSurfaceRepresentation() no reconvierte si la
    representacion ya existe: despues de tocar el labelmap puede devolver la
    malla vieja. Hay que borrarla primero."""
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


def generar():
    """BLOQUE A v7. Segmenta, decide y deja las piezas listas para la
    revision. No fusiona nada: eso lo hace confirmar_craneo()."""
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
    nombre = volumeNode.GetName().lower()
    if ("cerebro" in nombre or "brain" in nombre) and "hueso" not in nombre:
        print("CranioPlan: AVISO - serie de kernel BLANDO. El hueso delgado "
              "pierde HU. Si el estudio tiene serie de HUESO, preferila.")

    _limpiarEscenaPrevia()

    volArr = slicer.util.arrayFromVolume(volumeNode)
    espaciadoZYX = (esp[2], esp[1], esp[0])

    labels, filas, info = analizarVolumen(volArr, espaciadoZYX)
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

    referencia, aceptadas, dudosas = decidir(labels, filas, espaciadoZYX)
    if referencia is None:
        print("CranioPlan: ninguna isla parece hueso. Revisa los umbrales HU o "
              "el volumen elegido.")
        _imprimirTabla(filas)
        return None

    _imprimirTabla(filas)

    # --- construir la segmentacion ---
    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_SEGMENTACION)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentacion = segNode.GetSegmentation()

    buf = np.zeros(labels.shape, dtype=np.uint8)
    lista = []
    for i, f in enumerate(aceptadas + dudosas, start=1):
        etiquetaZona = "OK" if f["zona"] == "ACEPTADA" else "DUDOSA"
        nombreSeg = "%02d_%s_%.2fcm3" % (i, etiquetaZona, f["vol"])
        segId = segmentacion.AddEmptySegment("", nombreSeg)
        caja = f["caja"]
        buf[caja] = (labels[caja] == f["etiqueta"]).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(buf, segNode, segId, volumeNode)
        buf[caja] = 0
        color = COLOR_ACEPTADA if f["zona"] == "ACEPTADA" else COLOR_DUDOSA
        segmentacion.GetSegment(segId).SetColor(*color)
        _estado["colores"][segId] = color
        lista.append({
            "numero": i, "segId": segId, "vol": f["vol"],
            "distRef": f["distRef"], "distMasa": f["distMasa"], "env": f["env"],
            "fracCort": f["fracCort"], "p95": f["p95"], "zona": f["zona"],
            "nota": f["motivo"], "referencia": f["referencia"],
        })
    del buf

    _forzarSuperficie(segNode)

    _estado["labels"] = labels
    _estado["espaciadoZYX"] = espaciadoZYX
    _estado["filas"] = filas
    _estado["candidatas"] = lista
    _estado["segNode"] = segNode

    print("")
    print("--- PIEZAS EN LA ESCENA ---")
    for d in lista:
        marca = " *referencia*" if d["referencia"] else ""
        extra = ("   <- " + d["nota"]) if d["nota"] else ""
        print("  [%2d] %-8s vol=%8.2f cm3  cort=%5.1f%%  env=%5.1f%%  "
              "dRef=%5.1f mm%s%s"
              % (d["numero"], d["zona"], d["vol"], d["fracCort"] * 100.0,
                 (d["env"] or 0.0) * 100.0,
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
    print("  env% = fraccion de direcciones, vistas desde el centro del craneo,")
    print("         en las que aparece la isla. SOLO SE MIDE, no decide nada.")
    print("")
    print("   vol(cm3)  p95(HU)  cort%   env%   dRef   dMasa  espMin  %FOV  borde"
          "  zona       motivo")
    for f in filas:
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "   -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "    -" if f["distMasa"] is None else "%6.1f" % f["distMasa"]
        marca = "  *ref*" if f["referencia"] else ""
        print("  %9.2f  %7.0f  %5.1f  %s  %s  %s  %6.1f  %4.0f%%  %-5s  %-9s %s%s"
              % (f["vol"], f["p95"], f["fracCort"] * 100.0, env, dr, dm,
                 min(f["extMM"]), max(f["fracFOV"]) * 100,
                 "si" if f["tocaBorde"] else "no",
                 f["zona"] or "?", f["motivo"], marca))
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
    segmentacion = segNode.GetSegmentation()
    for d in _estado["candidatas"]:
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
              "queres borrarla, revisa antes que el volumen elegido sea el "
              "correcto.")
        return
    segNode.GetSegmentation().RemoveSegment(d["segId"])
    _estado["candidatas"] = [x for x in _estado["candidatas"] if x["numero"] != numero]
    print("CranioPlan: pieza %d eliminada. Quedan %d."
          % (numero, len(_estado["candidatas"])))


def eliminar_dudosas():
    """Borra de una sola vez todas las marcadas como dudosas. Es un atajo:
    en los 6 casos relevados casi todas eran vertebras cervicales."""
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
    print("  listar_volumenes()          -> lista las series cargadas")
    print("  usar_volumen(i)             -> elige serie por indice o nombre")
    print("  generar()                   -> corre el Bloque A v7")
    print("  solo_dudosas()              -> muestra solo lo no decidido")
    print("  resaltar(n)                 -> pinta esa pieza de rojo")
    print("  eliminar(n)                 -> borra esa pieza")
    print("  eliminar_dudosas()          -> borra todas las dudosas de una")
    print("  mostrar_todas()             -> restaura colores")
    print("  confirmar_craneo()          -> fusiona en '%s'" % NOMBRE_SEGMENTO_CRANEO)
    print("  enviar_a_planner()          -> exporta el modelo '%s'" % NOMBRE_MODELO_CRANEO)
    print("  verificar_compatibilidad_v9() -> chequea antes de cortar()")


generar()
