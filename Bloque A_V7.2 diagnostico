# ============================================================
# CRANIOPLAN - BLOQUE A v7.2  ::  MODO DIAGNOSTICO
#   script de consola para 3D Slicer 5.10 - NO BORRA NADA
# ============================================================
#
# PARA QUE SIRVE
# --------------
# Relevar una tomografia nueva y sumarla a la evidencia con la que se
# calibran los criterios del Bloque A. NO es un paso obligatorio del uso
# normal: para producir el craneo de un paciente se corre directamente
# BloqueA_v72_consola.py. Esto se usa cuando se quiere que el caso quede
# medido, etiquetado y trazable en el CSV.
#
# ------------------------------------------------------------
# LO MAS IMPORTANTE DE ESTA VERSION
# ------------------------------------------------------------
# El bloque "NUCLEO DE DECISION" de este archivo es IDENTICO al del modulo
# BloqueA_v72_consola.py. No es una copia parecida: es el mismo codigo. En
# v6-DIAG y v7.1-DIAG la logica estaba reimplementada aparte, y eso significa
# que con el tiempo las dos podian derivar y el relevamiento habria dejado de
# medir lo que el modulo realmente hace. Ahora el veredicto que ves aca es,
# por construccion, el veredicto que el modulo va a dar.
#
# La unica diferencia entre los dos archivos es lo que hacen DESPUES de
# decidir: el modulo crea solo las aceptadas y las dudosas y sigue hacia el
# corte; este crea TODAS las islas, las pinta segun el veredicto, y agrega el
# etiquetado manual y la exportacion a CSV.
#
# ------------------------------------------------------------
# FLUJO POR CASO
# ------------------------------------------------------------
#   listar_volumenes()
#   usar_volumen("Hueso")            # elegi la serie
#   caso("MESKY 912")                # <-- NUEVO: nombre del caso para el CSV
#   diagnostico()
#   barrido_hu()
#   crudo()                          # opcional
#   ir_a(4) ; solo(4) ; et(4, "mandibula")   # etiquetar isla por isla
#   pendientes()
#   reporte()
#   csv()
#   limpiar()
#
# ------------------------------------------------------------
# ETIQUETAS SUGERIDAS (sin acentos, para poder agrupar)
# ------------------------------------------------------------
#   hueso_craneal  base_craneo  mandibula  maxilar  diente
#   vertebra  hueso_otro
#   camilla  almohadilla  banda  electrodo  chupete  tubo  metal  ruido
#   no_se      <- usarla sin culpa, es informacion valida
#
# Sobre almohadilla vs ruido: para la calibracion la distincion NO cambia
# nada, porque las dos cuentan como "no es hueso". Si dudas, mira la cortical
# y la extension en la tabla: espuma y gel dan cortical casi 0; un riel o una
# banda dan cortical media pero extension grande y compacidad baja. Si
# igual no te decidis, poné no_se.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS
# ------------------------------------------------------------
# - Crea como maximo MAX_SEGMENTOS islas; el resto se agrupa en
#   "99_Resto_pequenas". Con 793 islas, crearlas todas cuelga Slicer.
# - Las distancias se miden superficie-superficie con submuestreo. Error
#   tipico < 1 mm: sirve para decidir "pegada" vs "lejos", no para medir el
#   ancho de una sutura.
# - Si la referencia sale mal elegida, dRef, dMasa y env no significan nada.
#   Mira PRIMERO la tabla de candidatas; si eligio mal, usar_referencia(n).
#   Esto ya nos paso: en la TC 175 (RUIZ) el v6 eligio una vertebra como
#   referencia y hubo que descartar el caso entero de la calibracion.
# ============================================================


# ============================================================
# CONFIG
# ============================================================

NOMBRE_VOLUMEN = None

HU_MIN = 300
HU_MAX = 3000

# --- Los mismos umbrales del modulo v7.2. No tocar sin motivo: si difieren,
#     el diagnostico deja de predecir lo que el modulo hace. ---
UMBRAL_CORTICAL_HU = 700.0
FRACCION_CORTICAL_MINIMA = 0.05
COMPACIDAD_MINIMA_REFERENCIA = 16.0
FACTOR_EXTENSION_MAXIMA = 1.25
COMPACIDAD_MINIMA_PIEZA = 6.0
DIST_ACEPTACION_MM = 2.0
DIST_DUDOSA_MM = 15.0
MODO_POSTOP = False
DIST_POSTOP_ACEPTACION_MM = 60.0
VOLUMEN_CRANEO_MINIMO_CM3 = 20.0
ESPESOR_LOSA_MM = 30.0
FRACCION_LOSA_FOV = 0.85
ENV_BINS_COSTHETA = 12
ENV_BINS_PHI = 24

# Piso de listado. Mas bajo que el del modulo a proposito: queremos VER lo
# que el modulo ni analiza.
VOL_MINIMO_CM3 = 0.05

MAX_SEGMENTOS = 30
MAX_PUNTOS_KDTREE = 40000
BARRIDO_HU = [150, 200, 250, 300, 350, 400]

NOMBRE_NODO_DIAG = "DIAG_Islas"
NOMBRE_NODO_CRUDO = "DIAG_Crudo"
SUAVIZADO_SUPERFICIE = "0.3"

# CSV acumulativo. Esquema nuevo respecto del v7.1: agrega caso y compacidad.
RUTA_CSV = None   # None = ~/Documents/CranioPlan_diag_islas_v72.csv


# ============================================================
# CODIGO
# ============================================================

import os
import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

CRANIOPLAN_DIAG = "v7.2a-DIAG (2026-08-18) consola - no destructivo"

_d = {
    "volumeNode": None, "labels": None, "espaciadoZYX": None,
    "filas": [], "creadas": [], "segNode": None, "crudoNode": None,
    "colores": {}, "huMin": None, "totalCM3": None, "nIslas": None,
    "referencia": None, "caso": "",
}

COLOR = {
    "ACEPTADA": (0.90, 0.80, 0.60),   # hueso
    "DUDOSA":   (1.00, 0.60, 0.10),   # naranja
    "MATERIAL": (0.25, 0.75, 0.70),   # turquesa
    "LOSA":     (0.25, 0.45, 0.95),   # azul
    "EXTENSION":(0.95, 0.20, 0.20),   # rojo
    "DIFUSA":   (0.60, 0.60, 0.30),   # oliva: ruido de reconstruccion
    "LEJOS":    (0.70, 0.35, 0.90),   # violeta
}

_ETIQUETAS = ["hueso_craneal", "base_craneo", "mandibula", "maxilar", "diente",
              "vertebra", "hueso_otro", "camilla", "almohadilla", "banda",
              "electrodo", "chupete", "tubo", "metal", "ruido", "no_se"]
_CRANEO = {"hueso_craneal", "base_craneo", "mandibula", "maxilar", "diente"}
_BASURA = {"camilla", "almohadilla", "banda", "electrodo", "chupete",
           "tubo", "metal", "ruido"}


def _colorDe(f):
    if f["veredicto"] == "ACEPTADA":
        return COLOR["ACEPTADA"]
    if f["veredicto"] == "DUDOSA":
        return COLOR["DUDOSA"]
    return COLOR.get(f["categoria"], (0.55, 0.55, 0.55))


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
# Volumen y escena
# ------------------------------------------------------------

def listar_volumenes():
    nodos = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
    if not nodos:
        print("DIAG: no hay volumenes cargados.")
        return []
    print("DIAG: volumenes en la escena:")
    for i, n in enumerate(nodos):
        dims = n.GetImageData().GetDimensions() if n.GetImageData() else (0, 0, 0)
        esp = n.GetSpacing()
        print("  [%d] %-45s  %dx%dx%d  %.3f/%.3f/%.3f mm"
              % (i, n.GetName(), dims[0], dims[1], dims[2], esp[0], esp[1], esp[2]))
    return nodos


def usar_volumen(referencia=None):
    nodos = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
    if not nodos:
        print("DIAG: no hay volumenes cargados.")
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
        print("DIAG: no encontre ese volumen. Usa listar_volumenes().")
        return None
    _d["volumeNode"] = elegido
    print("DIAG: volumen de trabajo -> %s" % elegido.GetName())
    return elegido


def caso(nombre):
    """
    Registra el nombre del caso. IMPRESCINDIBLE antes de csv().

    Motivo: los nombres de serie NO son unicos entre pacientes. En el
    relevamiento anterior, MESKY y el postoperatorio de MORENO se llamaban
    los dos '4: Cerebro Ped 0.5 Vol Vol.' y quedaron 43 filas mezcladas en el
    CSV que solo se podian separar porque la numeracion se reinicia. Con esta
    columna eso no vuelve a pasar.
    """
    _d["caso"] = str(nombre).strip()
    print("DIAG: caso -> '%s'" % _d["caso"])
    return _d["caso"]


def _volumen():
    v = _d.get("volumeNode")
    if v is None or slicer.mrmlScene.GetNodeByID(v.GetID()) is None:
        v = usar_volumen(NOMBRE_VOLUMEN)
    return v


def _encabezado(volumeNode):
    dims = volumeNode.GetImageData().GetDimensions()
    esp = volumeNode.GetSpacing()
    print("DIAG: caso    : %s" % (_d.get("caso") or "<<< SIN NOMBRE: usa caso('...') >>>"))
    print("DIAG: volumen : %s" % volumeNode.GetName())
    print("DIAG: dims    : %d x %d x %d voxels" % dims)
    print("DIAG: spacing : %.3f x %.3f x %.3f mm" % esp)
    print("DIAG: FOV     : %.1f x %.1f x %.1f mm"
          % (dims[0] * esp[0], dims[1] * esp[1], dims[2] * esp[2]))
    n = volumeNode.GetName().lower()
    if ("cerebro" in n or "brain" in n) and "hueso" not in n:
        print("DIAG: *** AVISO - serie de KERNEL BLANDO (cerebro). El hueso")
        print("DIAG:     delgado pierde HU. Si el estudio tiene serie de HUESO,")
        print("DIAG:     corre tambien sobre esa y compara 'TOTAL en mascara'.")


def _borrarNodosPorNombre(nombre):
    nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    while nodo is not None:
        slicer.mrmlScene.RemoveNode(nodo)
        nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)


def _forzarSuperficie(segNode):
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


# ------------------------------------------------------------
# Comandos principales
# ------------------------------------------------------------

def diagnostico(referenciaForzada=None):
    """Analiza, aplica la MISMA decision que el modulo v7.2 y crea TODAS las
    islas relevantes como segmentos. No elimina ninguna."""
    volumeNode = _volumen()
    if volumeNode is None:
        return None
    _d["huMin"] = HU_MIN

    print("=" * 66)
    print("CranioPlan - Bloque A %s" % CRANIOPLAN_DIAG)
    print("DIAG: MODO NO DESTRUCTIVO - no se elimina ninguna isla.")
    if MODO_POSTOP:
        print("DIAG: *** MODO POSTOPERATORIO ACTIVO ***")
    _encabezado(volumeNode)

    volArr = slicer.util.arrayFromVolume(volumeNode)
    esp = volumeNode.GetSpacing()
    espaciadoZYX = (esp[2], esp[1], esp[0])

    labels, filas, info = analizarVolumen(volArr, espaciadoZYX, HU_MIN, HU_MAX)
    if labels is None or not filas:
        print("DIAG: no hay islas por encima de %.2f cm3 en %d-%d HU."
              % (VOL_MINIMO_CM3, HU_MIN, HU_MAX))
        return None

    _d.update({"nIslas": info["nIslas"], "totalCM3": info["totalCM3"],
               "labels": labels, "espaciadoZYX": espaciadoZYX, "filas": filas})
    print("DIAG: islas detectadas : %d  (%d por encima de %.2f cm3)"
          % (info["nIslas"], len(filas), VOL_MINIMO_CM3))
    print("DIAG: TOTAL en mascara : %.2f cm3" % info["totalCM3"])
    if info["valorRelleno"] is not None:
        print("DIAG: relleno fuera del FOV en %d HU (borde circular activo)."
              % info["valorRelleno"])
    else:
        print("DIAG: sin relleno detectable: el borde circular NO se testea.")

    referencia, aceptadas, dudosas, relajado = decidir(
        labels, filas, espaciadoZYX, referenciaForzada=referenciaForzada)
    if referencia is None:
        print("DIAG: ninguna isla puede ser referencia. Revisa HU o la serie.")
        return None
    _d["referencia"] = referencia

    print("")
    print(tablaDeCandidatas(filas, referencia))
    print("DIAG: referencia -> %.2f cm3, cortical %.1f%%, compacidad %.1f, "
          "ext %.0f mm, env %.1f%%"
          % (referencia["vol"], referencia["fracCort"] * 100.0,
             referencia["compac"], max(referencia["extMM"]),
             (referencia["env"] or 0) * 100.0))
    if relajado:
        print("DIAG: OJO - se relajo el filtro de %s. Mira la referencia con "
              "desconfianza." % relajado)

    _imprimirTabla(filas)
    _construirEscena(labels, volumeNode, filas)
    ayuda()
    return filas


def usar_referencia(numero):
    """Rehace el diagnostico forzando la referencia, sin re-etiquetar el
    volumen. Las etiquetas manuales ya puestas se conservan."""
    if not _d.get("filas"):
        print("DIAG: primero corre diagnostico().")
        return None
    guardadas = {f["etiqueta"]: (f.get("etiqueta_real", ""), f.get("nota_real", ""))
                 for f in _d["filas"]}
    labels, filas = _d["labels"], _d["filas"]
    referencia, aceptadas, dudosas, relajado = decidir(
        labels, filas, _d["espaciadoZYX"], referenciaForzada=int(numero))
    if referencia is None:
        return None
    for f in filas:
        et, nota = guardadas.get(f["etiqueta"], ("", ""))
        f["etiqueta_real"], f["nota_real"] = et, nota
    _d["referencia"] = referencia
    print("")
    print(tablaDeCandidatas(filas, referencia))
    _imprimirTabla(filas)
    _construirEscena(labels, _d["volumeNode"], filas)
    return filas


def postop(activar=True):
    """Modo postoperatorio: sube la distancia de aceptacion para recuperar los
    colgajos ya cortados. Rehace la decision sin re-segmentar."""
    global MODO_POSTOP
    MODO_POSTOP = bool(activar)
    print("DIAG: MODO POSTOPERATORIO %s (aceptacion %.0f mm)."
          % ("ACTIVADO" if MODO_POSTOP else "desactivado",
             DIST_POSTOP_ACEPTACION_MM if MODO_POSTOP else DIST_ACEPTACION_MM))
    if _d.get("filas"):
        return usar_referencia(
            _d["filas"].index(_d["referencia"]) + 1 if _d.get("referencia") else 1)
    return None


def _construirEscena(labels, volumeNode, filas):
    _borrarNodosPorNombre(NOMBRE_NODO_DIAG)
    _d["colores"] = {}
    _d["creadas"] = []
    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_DIAG)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentacion = segNode.GetSegmentation()

    aCrear, resto = filas[:MAX_SEGMENTOS], filas[MAX_SEGMENTOS:]
    buf = np.zeros(labels.shape, dtype=np.uint8)
    creadas = []
    for i, f in enumerate(aCrear, start=1):
        etq = f["veredicto"] if f["veredicto"] == "ACEPTADA" or f["veredicto"] == "DUDOSA" \
              else f["categoria"]
        segId = segmentacion.AddEmptySegment("", "%02d_%s_%.2fcm3" % (i, etq, f["vol"]))
        caja = f["caja"]
        buf[caja] = (labels[caja] == f["etiqueta"]).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(buf, segNode, segId, volumeNode)
        buf[caja] = 0
        color = _colorDe(f)
        segmentacion.GetSegment(segId).SetColor(*color)
        _d["colores"][segId] = color
        f["numero"] = i
        f["segId"] = segId
        creadas.append(f)

    if resto:
        for f in resto:
            caja = f["caja"]
            buf[caja] |= (labels[caja] == f["etiqueta"]).astype(np.uint8)
        segId = segmentacion.AddEmptySegment("", "99_Resto_pequenas_%d" % len(resto))
        slicer.util.updateSegmentBinaryLabelmapFromArray(buf, segNode, segId, volumeNode)
        segmentacion.GetSegment(segId).SetColor(0.40, 0.40, 0.40)
        buf[:] = 0
        print("DIAG: %d isla(s) fuera del tope de %d agrupadas en "
              "'99_Resto_pequenas'." % (len(resto), MAX_SEGMENTOS))
    del buf
    _forzarSuperficie(segNode)
    _d["segNode"] = segNode
    _d["creadas"] = creadas
    print("")
    print("DIAG: nodo '%s' con %d segmento(s)."
          % (NOMBRE_NODO_DIAG, segmentacion.GetNumberOfSegments()))
    print("DIAG: hueso=ACEPTADA  naranja=DUDOSA  turquesa=MATERIAL  azul=LOSA")
    print("DIAG: rojo=EXTENSION  violeta=LEJOS")
    return creadas


def _imprimirTabla(filas):
    print("")
    print("--- TODAS LAS ISLAS (>= %.2f cm3) - NINGUNA FUE ELIMINADA ---"
          % VOL_MINIMO_CM3)
    print("  compac = volumen / eje_mayor^3 x1e6. Craneos medidos 24-78;")
    print("           almohadilla envolvente 11; riel de camilla 1.")
    print("  env%   = se mide, no decide nada.")
    print("")
    print("  #   vol(cm3)   p50   p95    max  cort%  compac   env%   dRef  dMasa"
          "  espMin   ext  %FOV  borde  veredicto  motivo")
    for i, f in enumerate(filas, start=1):
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "    -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "    -" if f["distMasa"] is None else "%5.1f" % f["distMasa"]
        marca = " *ref*" if f["referencia"] else ""
        print("%3d %9.2f %5.0f %5.0f %6.0f %6.1f %7.1f %s %s %s %7.1f %5.0f %4.0f%%"
              "  %-5s  %-9s %s%s"
              % (i, f["vol"], f["p50"], f["p95"], f["max"], f["fracCort"] * 100.0,
                 f["compac"], env, dr, dm, min(f["extMM"]), max(f["extMM"]),
                 max(f["fracFOV"]) * 100, "si" if f["tocaBorde"] else "no",
                 f["veredicto"] or "?", f["motivo"], marca))
    porZona = {}
    for f in filas:
        clave = f["veredicto"] if f["veredicto"] in ("ACEPTADA", "DUDOSA") else f["categoria"]
        porZona.setdefault(clave, []).append(f["vol"])
    print("")
    print("  TOTAL en mascara : %8.2f cm3" % (_d.get("totalCM3") or 0.0))
    for z in ["ACEPTADA", "DUDOSA", "MATERIAL", "DIFUSA", "LOSA", "EXTENSION",
              "LEJOS"]:
        if z in porZona:
            print("  %-9s        : %8.2f cm3  en %d isla(s)"
                  % (z, sum(porZona[z]), len(porZona[z])))
    print("-" * 74)


def barrido_hu(umbrales=None):
    """Cuanto hueso aparece o desaparece al mover HU_MIN. No crea nada."""
    volumeNode = _volumen()
    if volumeNode is None:
        return
    lista = BARRIDO_HU if umbrales is None else umbrales
    esp = volumeNode.GetSpacing()
    volVoxCM3 = (esp[0] * esp[1] * esp[2]) / 1000.0
    volArr = slicer.util.arrayFromVolume(volumeNode)
    estructura = ndimage.generate_binary_structure(3, 3)
    print("")
    print("--- BARRIDO DE HU_MIN (HU_MAX = %d) ---" % HU_MAX)
    print("  HU_MIN   TOTAL(cm3)   n_islas   isla_mayor(cm3)   resto(cm3)")
    for h in lista:
        m = (volArr >= h) & (volArr <= HU_MAX)
        total = int(np.count_nonzero(m)) * volVoxCM3
        lab, n = ndimage.label(m, structure=estructura)
        if n == 0:
            print("  %6d   %10.2f   %7d   %15s   %10s" % (h, total, 0, "-", "-"))
            continue
        conteos = np.bincount(lab.ravel())
        conteos[0] = 0
        mayor = int(conteos.max()) * volVoxCM3
        print("  %6d   %10.2f   %7d   %15.2f   %10.2f"
              % (h, total, n, mayor, total - mayor))
        del lab, m
    print("---------------------------------------")


def crudo(huMin=None, huMax=None):
    """Un unico segmento con TODO lo que cae en el rango HU, sin ningun
    filtro. Es el ground truth visual: si falta hueso aca, el problema es el
    umbral o el kernel de la serie, no el algoritmo de seleccion."""
    volumeNode = _volumen()
    if volumeNode is None:
        return None
    hmin = HU_MIN if huMin is None else huMin
    hmax = HU_MAX if huMax is None else huMax
    volArr = slicer.util.arrayFromVolume(volumeNode)
    esp = volumeNode.GetSpacing()
    volVoxCM3 = (esp[0] * esp[1] * esp[2]) / 1000.0
    mascara = ((volArr >= hmin) & (volArr <= hmax)).astype(np.uint8)
    nVox = int(np.count_nonzero(mascara))
    print("DIAG CRUDO: %d voxeles = %.2f cm3 en %d-%d HU."
          % (nVox, nVox * volVoxCM3, hmin, hmax))
    if nVox == 0:
        return None
    _borrarNodosPorNombre(NOMBRE_NODO_CRUDO)
    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_CRUDO)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segId = segNode.GetSegmentation().AddEmptySegment("", "Crudo_%dHU" % hmin)
    slicer.util.updateSegmentBinaryLabelmapFromArray(mascara, segNode, segId, volumeNode)
    segNode.GetSegmentation().GetSegment(segId).SetColor(0.90, 0.80, 0.60)
    _forzarSuperficie(segNode)
    _d["crudoNode"] = segNode
    print("DIAG: nodo '%s' creado. Oculta 'DIAG_Islas' para verlo."
          % NOMBRE_NODO_CRUDO)
    return segNode


# ------------------------------------------------------------
# Inspeccion
# ------------------------------------------------------------

def _porNumero(numero):
    return next((f for f in _d.get("creadas", []) if f.get("numero") == numero), None)


def resaltar(numero):
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    segmentacion = segNode.GetSegmentation()
    for f in _d["creadas"]:
        seg = segmentacion.GetSegment(f["segId"])
        if seg is None:
            continue
        seg.SetColor(1.0, 0.05, 0.05) if f["numero"] == numero \
            else seg.SetColor(*_d["colores"].get(f["segId"], (0.6, 0.6, 0.6)))


def solo(numero):
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    disp = segNode.GetDisplayNode()
    for f in _d["creadas"]:
        disp.SetSegmentVisibility(f["segId"], f["numero"] == numero)


def mostrar_todas():
    segNode = _d.get("segNode")
    if segNode is None:
        return
    disp = segNode.GetDisplayNode()
    segmentacion = segNode.GetSegmentation()
    for f in _d["creadas"]:
        disp.SetSegmentVisibility(f["segId"], True)
        seg = segmentacion.GetSegment(f["segId"])
        if seg is not None:
            seg.SetColor(*_d["colores"].get(f["segId"], (0.6, 0.6, 0.6)))


def ocultar_aceptadas():
    """Deja a la vista todo lo que el v7.2 NO acepta (dudosas incluidas)."""
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    disp = segNode.GetDisplayNode()
    n = 0
    for f in _d["creadas"]:
        v = f["veredicto"] != "ACEPTADA"
        disp.SetSegmentVisibility(f["segId"], v)
        n += 1 if v else 0
    print("DIAG: %d isla(s) NO aceptadas a la vista." % n)


def ver_veredicto(nombre):
    """Filtra por veredicto o por categoria: ACEPTADA, DUDOSA, MATERIAL,
    LOSA, EXTENSION, LEJOS."""
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    disp = segNode.GetDisplayNode()
    clave = nombre.upper()
    n = 0
    for f in _d["creadas"]:
        v = (f["veredicto"] == clave) or (f["categoria"].upper() == clave)
        disp.SetSegmentVisibility(f["segId"], v)
        n += 1 if v else 0
    print("DIAG: %d isla(s) con veredicto/categoria %s." % (n, clave))


def ir_a(numero):
    """Centra las tres vistas 2D en el centroide de esa isla."""
    volumeNode = _d.get("volumeNode")
    f = _porNumero(numero)
    if f is None or volumeNode is None:
        print("DIAG: no existe la isla %s." % numero)
        return
    m = vtk.vtkMatrix4x4()
    volumeNode.GetIJKToRASMatrix(m)
    cx, cy, cz = f["ijk"]
    ras = [0.0, 0.0, 0.0, 0.0]
    m.MultiplyPoint([cx, cy, cz, 1.0], ras)
    slicer.vtkMRMLSliceNode.JumpAllSlices(
        slicer.mrmlScene, ras[0], ras[1], ras[2],
        slicer.vtkMRMLSliceNode.CenteredJumpSlice)
    print("DIAG: isla %d -> RAS (%.1f, %.1f, %.1f) | %.2f cm3, cort %.1f%%, "
          "compac %.1f, ext %.0f mm, dRef %.1f mm -> %s (%s)"
          % (numero, ras[0], ras[1], ras[2], f["vol"], f["fracCort"] * 100.0,
             f["compac"], max(f["extMM"]),
             f["distRef"] if f["distRef"] is not None else -1,
             f["veredicto"], f["categoria"]))


# ------------------------------------------------------------
# Etiquetado manual
# ------------------------------------------------------------

def etiquetas():
    print("DIAG: etiquetas sugeridas:")
    for e in _ETIQUETAS:
        marca = "   (cuenta como craneo)" if e in _CRANEO else (
                "   (cuenta como basura)" if e in _BASURA else "")
        print("   %s%s" % (e, marca))


def et(numero, etiqueta, nota=""):
    """Registra que es la isla. Alias corto de etiquetar()."""
    return etiquetar(numero, etiqueta, nota)


def etiquetar(numero, etiqueta, nota=""):
    f = _porNumero(numero)
    if f is None:
        print("DIAG: no existe la isla %s en la escena." % numero)
        return
    clave = etiqueta.strip().lower().replace(" ", "_")
    if clave not in _ETIQUETAS:
        print("DIAG: '%s' no esta en la lista sugerida. Se guarda igual, pero si "
              "es un tipeo conviene corregirlo (etiquetas())." % clave)
    f["etiqueta_real"] = clave
    f["nota_real"] = nota
    esperado = "ACEPTADA" if clave in _CRANEO else "no-ACEPTADA"
    real = "ACEPTADA" if f["veredicto"] == "ACEPTADA" else "no-ACEPTADA"
    flecha = "" if (clave == "no_se" or esperado == real) else "   <-- DISCREPANCIA"
    print("DIAG: isla %d = %s  (v7.2 dijo %s | %.2f cm3, cort %.1f%%, compac %.1f, "
          "dRef %.1f mm)%s"
          % (numero, clave, f["veredicto"], f["vol"], f["fracCort"] * 100.0,
             f["compac"], f["distRef"] if f["distRef"] is not None else -1, flecha))


def pendientes():
    faltan = [f for f in _d.get("creadas", []) if not f.get("etiqueta_real")]
    if not faltan:
        print("DIAG: todas las islas estan etiquetadas.")
        return []
    print("DIAG: faltan etiquetar %d isla(s):" % len(faltan))
    for f in faltan:
        print("   %2d  %8.2f cm3  cort %5.1f%%  compac %6.1f  dRef %6.1f  -> %s"
              % (f["numero"], f["vol"], f["fracCort"] * 100.0, f["compac"],
                 f["distRef"] if f["distRef"] is not None else -1, f["veredicto"]))
    return faltan


def discrepancias():
    """Las unicas filas que importan para calibrar. Un falso negativo (hueso
    craneal no aceptado) es grave; un falso positivo (basura aceptada) es
    molesto pero lo corrige el cirujano en un click."""
    fn, fp = [], []
    for f in _d.get("creadas", []):
        e = f.get("etiqueta_real")
        if not e or e == "no_se":
            continue
        if e in _CRANEO and f["veredicto"] != "ACEPTADA":
            fn.append(f)
        elif e in _BASURA and f["veredicto"] == "ACEPTADA":
            fp.append(f)
    print("")
    print("--- DISCREPANCIAS v7.2 vs realidad ---")
    print("  FALSOS NEGATIVOS (hueso craneal NO aceptado): %d" % len(fn))
    for f in fn:
        print("     %2d  %-14s %7.2f cm3  cort %5.1f%%  compac %6.1f  dRef %6.1f mm"
              "  -> %s (%s)"
              % (f["numero"], f["etiqueta_real"], f["vol"], f["fracCort"] * 100.0,
                 f["compac"], f["distRef"] if f["distRef"] is not None else -1,
                 f["veredicto"], f["categoria"]))
    print("  FALSOS POSITIVOS (basura aceptada): %d" % len(fp))
    for f in fp:
        print("     %2d  %-14s %7.2f cm3  cort %5.1f%%  compac %6.1f  dRef %6.1f mm"
              % (f["numero"], f["etiqueta_real"], f["vol"], f["fracCort"] * 100.0,
                 f["compac"], f["distRef"] if f["distRef"] is not None else -1))
    print("  (vertebras y 'no_se' no cuentan en ninguna de las dos listas)")
    print("--------------------------------------")
    return fn, fp


def metricas_por_etiqueta():
    """Resume las metricas de decision agrupadas por lo que la isla ES en
    realidad. Es la tabla que sostiene o mueve los umbrales."""
    grupos = {}
    for f in _d.get("creadas", []):
        e = f.get("etiqueta_real")
        if not e:
            continue
        grupos.setdefault(e, []).append(f)
    if not grupos:
        print("DIAG: todavia no hay islas etiquetadas.")
        return
    print("")
    print("--- METRICAS POR TIPO DE OBJETO (este caso) ---")
    print("  etiqueta          n   cort% min/max    compac min/max    dRef min/max"
          "    env% max")
    for e in sorted(grupos, key=lambda k: -max(x["vol"] for x in grupos[k])):
        g = grupos[e]
        co = [x["fracCort"] * 100.0 for x in g]
        cp = [x["compac"] for x in g]
        dr = [x["distRef"] for x in g if x["distRef"] is not None] or [-1]
        en = [(x["env"] or 0) * 100.0 for x in g]
        print("  %-14s %4d   %5.1f / %5.1f    %6.1f / %6.1f    %5.1f / %5.1f"
              "    %5.1f"
              % (e, len(g), min(co), max(co), min(cp), max(cp),
                 min(dr), max(dr), max(en)))
    print("-----------------------------------------------")


def reporte():
    """Bloque compacto, listo para copiar y pasar."""
    vol = _d.get("volumeNode")
    creadas = _d.get("creadas", [])
    ref = _d.get("referencia")
    if vol is None or not creadas:
        print("DIAG: primero corre diagnostico().")
        return
    dims = vol.GetImageData().GetDimensions()
    esp = vol.GetSpacing()
    print("")
    print("========== REPORTE CRANIOPLAN - BLOQUE A v7.2 DIAG ==========")
    print("caso         : %s" % (_d.get("caso") or "<<< SIN NOMBRE >>>"))
    print("serie        : %s" % vol.GetName())
    print("dims         : %d x %d x %d" % dims)
    print("spacing      : %.3f x %.3f x %.3f mm" % esp)
    print("FOV          : %.1f x %.1f x %.1f mm"
          % (dims[0] * esp[0], dims[1] * esp[1], dims[2] * esp[2]))
    print("HU_MIN       : %s   modo: %s"
          % (_d.get("huMin"), "POSTOP" if MODO_POSTOP else "preoperatorio"))
    print("islas        : %s" % _d.get("nIslas"))
    print("total mascara: %.2f cm3" % (_d.get("totalCM3") or 0.0))
    if ref is not None:
        print("referencia   : %.2f cm3, cort %.1f%%, compac %.1f, ext %.0f mm, "
              "toca borde: %s"
              % (ref["vol"], ref["fracCort"] * 100.0, ref["compac"],
                 max(ref["extMM"]), "si" if ref["tocaBorde"] else "no"))
    print("")
    print(tablaDeCandidatas(_d["filas"], ref))
    print("")
    print("  #  vol_cm3   p50   p95    max  cort% compac   env%   dRef  dMasa"
          " espMin   ext  %FOV borde  veredicto  categoria  ES_EN_REALIDAD")
    for f in sorted(creadas, key=lambda x: x["numero"]):
        e = f.get("etiqueta_real") or "SIN_ETIQUETAR"
        nota = f.get("nota_real", "")
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "    -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "    -" if f["distMasa"] is None else "%5.1f" % f["distMasa"]
        print("%3d %8.2f %5.0f %5.0f %6.0f %6.1f %6.1f %s %s %s %6.1f %5.0f %4.0f%%"
              "  %-4s  %-9s  %-10s %s%s"
              % (f["numero"], f["vol"], f["p50"], f["p95"], f["max"],
                 f["fracCort"] * 100.0, f["compac"], env, dr, dm,
                 min(f["extMM"]), max(f["extMM"]), max(f["fracFOV"]) * 100,
                 "si" if f["tocaBorde"] else "no", f["veredicto"],
                 f["categoria"], e, ("  (%s)" % nota) if nota else ""))
    print("")
    discrepancias()
    metricas_por_etiqueta()
    faltan = [f for f in creadas if not f.get("etiqueta_real")]
    if faltan:
        print("ATENCION: quedan %d isla(s) sin etiquetar (pendientes())."
              % len(faltan))
    if not _d.get("caso"):
        print("ATENCION: el caso no tiene nombre. Corre caso('...') antes de csv().")
    print("=============================================================")


def csv(ruta=None):
    """Exporta y ACUMULA. Un solo archivo con todas las TC es lo que sirve
    para calibrar; no lo pises entre pacientes."""
    creadas = _d.get("creadas", [])
    if not creadas:
        print("DIAG: primero corre diagnostico().")
        return
    if not _d.get("caso"):
        print("DIAG: falta el nombre del caso. Corre caso('MESKY 912') y repeti, "
              "o las filas quedan sin poder separarse de otro paciente.")
        return
    if ruta is None:
        ruta = RUTA_CSV or os.path.join(
            os.path.expanduser("~"), "Documents", "CranioPlan_diag_islas_v72.csv")
    vol = _d.get("volumeNode")
    nuevo = not os.path.exists(ruta)
    with open(ruta, "a") as fh:
        if nuevo:
            fh.write("caso;serie;modo;hu_min;n;vol_cm3;p50;p95;p99;max;"
                     "frac_cortical;compacidad;env;esp_min_mm;ext_max_mm;pct_fov;"
                     "toca_borde;es_losa;dist_ref_mm;dist_masa_mm;es_referencia;"
                     "veredicto_v72;categoria;etiqueta_real;nota\n")
        for f in sorted(creadas, key=lambda x: x["numero"]):
            fh.write("%s;%s;%s;%s;%d;%.3f;%.0f;%.0f;%.0f;%.0f;%.4f;%.2f;%s;%.2f;"
                     "%.2f;%.1f;%s;%s;%s;%s;%s;%s;%s;%s;%s\n"
                     % (_d["caso"], vol.GetName() if vol else "?",
                        "postop" if MODO_POSTOP else "preop", _d.get("huMin"),
                        f["numero"], f["vol"], f["p50"], f["p95"], f["p99"],
                        f["max"], f["fracCort"], f["compac"],
                        "" if f["env"] is None else "%.4f" % f["env"],
                        min(f["extMM"]), max(f["extMM"]), max(f["fracFOV"]) * 100,
                        "si" if f["tocaBorde"] else "no",
                        "si" if f["esLosa"] else "no",
                        "" if f["distRef"] is None else "%.2f" % f["distRef"],
                        "" if f["distMasa"] is None else "%.2f" % f["distMasa"],
                        "si" if f["referencia"] else "no",
                        f["veredicto"], f["categoria"],
                        f.get("etiqueta_real", ""), f.get("nota_real", "")))
    print("DIAG: %s -> %s" % ("creado" if nuevo else "agregado a", ruta))
    print("DIAG: esquema nuevo respecto del v7.1 (agrega caso, modo, compacidad,")
    print("      es_referencia y categoria), por eso el archivo cambia de nombre.")


def limpiar():
    _borrarNodosPorNombre(NOMBRE_NODO_DIAG)
    _borrarNodosPorNombre(NOMBRE_NODO_CRUDO)
    _d["segNode"] = None
    _d["crudoNode"] = None
    _d["creadas"] = []
    _d["caso"] = ""
    print("DIAG: escena limpia. Acordate de caso('...') para el proximo estudio.")


def ayuda():
    print("")
    print("Analisis:")
    print("  listar_volumenes() / usar_volumen(i)")
    print("  caso('MESKY 912')        -> nombre del caso, HACE FALTA para csv()")
    print("  diagnostico()            -> todas las islas, ninguna borrada")
    print("  usar_referencia(n)       -> forzar la referencia y volver a decidir")
    print("  postop()                 -> modo postoperatorio")
    print("  barrido_hu()             -> cuanto hueso aparece al bajar HU_MIN")
    print("  crudo()                  -> mascara HU sin ningun filtro")
    print("Inspeccion:")
    print("  ir_a(n) / solo(n) / resaltar(n) / mostrar_todas()")
    print("  ocultar_aceptadas()      -> ver todo lo que el v7.2 no acepta")
    print("  ver_veredicto('DUDOSA')  -> tambien MATERIAL, DIFUSA, LOSA,")
    print("                              EXTENSION, LEJOS")
    print("Etiquetado:")
    print("  etiquetas() / et(n, 'mandibula') / pendientes()")
    print("  discrepancias()          -> donde el v7.2 se equivoca")
    print("  metricas_por_etiqueta()  -> cort, compac, dRef y env por tipo de objeto")
    print("  reporte()                -> bloque compacto para pasar")
    print("  csv()                    -> exporta y ACUMULA todos los casos")
    print("  limpiar()                -> borrar los nodos DIAG_*")


diagnostico()
