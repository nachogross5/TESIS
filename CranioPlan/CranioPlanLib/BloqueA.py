# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - BLOQUE A v7.2  (version libreria, para el modulo)
#              PREPARACION DEL CRANEO
# ============================================================
# Este archivo es el script de consola BloqueA_v72_consola.py convertido en
# libreria. El NUCLEO DE DECISION (todo lo que esta entre las lineas de
# "NUCLEO" mas abajo) es el mismo codigo: mismas funciones, mismos umbrales,
# mismas reglas. Lo unico que cambio es el envoltorio:
#
#   - el diccionario global _estado paso a ser una instancia de BloqueA, asi
#     el modulo puede tener el estado de un paciente sin pisarlo con globals;
#   - print() paso a ser self.log(), para que los avisos lleguen tambien a la
#     interfaz y no solo a la consola de Python, que el medico no mira;
#   - MODO_POSTOP dejo de ser una variable global y es un argumento, porque
#     una variable global de modo es exactamente el tipo de estado invisible
#     que despues explica los "a mi ayer me funcionaba".
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
# LAS REGLAS, Y POR QUE
# ------------------------------------------------------------
# [1] LA REFERENCIA SE ELIGE POR CORTICAL + COMPACIDAD + VOLUMEN, sin cota
#     fija de tamano. Una cota en milimetros se rompe con la edad del
#     paciente; la compacidad (volumen / eje_mayor^3) es adimensional y no.
#     Craneos medidos: 24 a 78. Impostores: almohadilla de JUAN 11, riel de
#     camilla de JUSTINA601 1. El hueco entre 11 y 24 es el margen.
#
# [2] TOCAR EL BORDE DEL CAMPO NO ELIMINA NADA. Solo genera una NOTA DE
#     CALIDAD: el craneo puede venir recortado. En la TC 175 (RUIZ) el craneo
#     SI toca el borde, asi que vetarlo dejaba al craneo afuera.
#
# [3] EL FILTRO DE MATERIAL NO ES BORRADO DURO. En MORENO aparecio hueso
#     craneal real con 6.0 % de cortical contra un umbral de 5 %, y el chupete
#     de JUAN tiene 3.88 %. Los rangos se tocan: no existe umbral que deje el
#     hueso adentro y el chupete afuera. Una isla de cortical baja pegada a la
#     referencia pasa a DUDOSA en vez de eliminarse. Cuesta un click y no
#     puede perder una lamina de hueso.
#
# [4] ELIMINACION POR EXTENSION RELATIVA. Un fragmento no puede ser mas grande
#     que el craneo del que salio. Limite relativo (1.25x la referencia), asi
#     que escala con el paciente. Peor hueso no-referencia medido: 0.50x.
#
# [5] REFERENCIA VISIBLE Y FORZABLE. La interfaz muestra la tabla de
#     candidatas y permite corregir la eleccion.
#
# [6] MODO POSTOPERATORIO. En un craneo ya operado los colgajos estan
#     genuinamente separados (MORENO 501: 2.2 a 16.7 mm). Ningun umbral
#     preoperatorio razonable los recupera, y no deberia. El modo sube la
#     distancia de aceptacion y avisa que la revision manual es obligatoria.
#
# [7] COMPACIDAD MINIMA DE UNA PIEZA. El ruido de reconstruccion queda pegado
#     al craneo, asi que la regla blanda de material [3] lo conservaba. Se
#     separa por FORMA: ruido y soportes 0.5-2.1, hueso real no-referencia en
#     preoperatorio 14.3-157. Desactivada en postoperatorio (alli el colgajo
#     mas flaco de MORENO da 8.5).
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS (documentados a proposito, no ignorados)
# ------------------------------------------------------------
# [C1] LA COMPACIDAD ESTA CALIBRADA CON DOS IMPOSTORES. Siete craneos de un
#      lado (24 a 78) y dos objetos del otro (11 y 1). Mitigacion: la eleccion
#      de referencia es visible y forzable desde la interfaz.
# [C2] EL CHUPETE NO ES SEPARABLE POR MATERIAL. Queda DUDOSO. Es un click, a
#      cambio de no poder perder una lamina fina de hueso.
# [C3] LAS VERTEBRAS SON HUESO REAL. Ningun criterio las distingue de la
#      mandibula. Entran como aceptadas o dudosas y se sacan a mano. NO se
#      eliminan solas a proposito.
# [C4] SI LA BOVEDA VIENE PARTIDA EN DOS CASCARAS FINAS y ninguna llega a
#      compacidad 16, se relaja el filtro y se avisa. Nunca se observo.
# [C5] NUEVE SERIES DE SIETE PACIENTES, UN SOLO HOSPITAL.
# [C6] LA ENVOLVENCIA NO DECIDE NADA. Quedo descartada con datos.
# ============================================================

import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

from .Comun import (NOMBRE_NODO_SEGMENTACION, NOMBRE_SEGMENTO_CRANEO,
                    NOMBRE_MODELO_CRANEO)

BLOQUE_A_VERSION = "v7.2 (calibrado sobre 9 series del Garrahan)"

# ============================================================
# CONFIG
# ============================================================
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
#   ruido difuso medido ............ 1.5 - 2.1
#   rieles y almohadillas .......... 0.5 - 4.8
#   hueso real no-referencia PREOP . 14.3 - 157
COMPACIDAD_MINIMA_PIEZA = 6.0

# --- Distancias ---
DIST_ACEPTACION_MM = 2.0     # a la REFERENCIA. Hueso craneal medido: max 1.70
DIST_DUDOSA_MM     = 15.0    # a la masa aceptada. Mas alla, se elimina.

# --- Modo postoperatorio (ver [6]) ---
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

COLOR_ACEPTADA = (0.90, 0.80, 0.60)
COLOR_DUDOSA   = (1.00, 0.60, 0.10)
COLOR_RESALTE  = (1.00, 0.15, 0.15)


# ============================================================
# ============================================================
# NUCLEO DE DECISION
#
# Este bloque es el mismo codigo que corre el script de diagnostico con el
# que se hizo el relevamiento de los 9 casos. Es a proposito: si el
# diagnostico usara una copia parecida pero no igual, con el tiempo las dos
# derivarian y el relevamiento dejaria de medir lo que el modulo realmente
# hace. Si hay que tocar un umbral, se toca aca y en ningun otro lado.
# ============================================================
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
    fallar por memoria a la version 5.

    OJO: tocar el borde NO elimina ni condiciona nada. Se mide y se informa
    porque es un dato de CALIDAD DEL ESTUDIO util para el cirujano (el craneo
    puede venir recortado), no un criterio de segmentacion."""
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

    ES EL CRITERIO DE FORMA QUE REEMPLAZA A LA COTA FIJA DE TAMANO.

    Por que sirve: un craneo es una masa COMPACTA, mucho volumen metido en su
    propio tamano. Una almohadilla o un riel de camilla son ALARGADOS: se
    estiran mucho y ocupan poco. Medido sobre 7 craneos y 2 impostores:

        riel de camilla (JUSTINA601) ....  1
        almohadilla envolvente (JUAN) ... 11
        craneo de JUAN .................. 24   <- el mas flaco de los craneos
        craneo de MESKY ................. 46
        craneo de JUSTINA ............... 78

    Y lo importante: la magnitud es ADIMENSIONAL, o sea invariante de escala.
    El craneo de un chico de 4 anios es mas grande pero igual de compacto; la
    almohadilla sigue siendo alargada. Por eso este criterio no se rompe con
    la edad, que era el problema de la cota fija en milimetros.
    """
    e = max(extMM)
    if e <= 0.0:
        return 0.0
    return float(volCM3) / (e ** 3) * 1.0e6


def _envolvencia(ptsMM, centroMM):
    """
    Fraccion de direcciones del espacio, vistas desde el centro del craneo, en
    las que aparece esta isla. Bins de AREA IGUAL (uniformes en cos(theta) y
    en phi), no grilla lat-lon, que concentraria bins diminutos en los polos.

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
    mismo resultado que la transformada de distancia sobre el volumen completo
    con dos ordenes de magnitud menos de puntos."""
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


def analizarVolumen(volArr, espaciadoZYX, huMin=HU_MIN, huMax=HU_MAX):
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
        })
    filas.sort(key=lambda f: f["vol"], reverse=True)
    return labels, filas, {"nIslas": nIslas, "valorRelleno": valorRelleno,
                           "totalCM3": totalCM3}


def candidatasAReferencia(filas):
    """
    Admisibilidad para competir por REFERENCIA. Dos filtros, cada uno tapa el
    punto ciego del otro:

      cortical >= 5 %   -> saca espuma, gel, plastico, ruido y el riel de la
                           camilla de JUSTINA601.
      compacidad >= 16  -> saca la almohadilla envolvente de JUAN, que es el
                           unico impostor del dataset que le GANA en volumen
                           al craneo (79.5 vs 58.7 cm3). Con el criterio de
                           "mayor volumen" a secas, la referencia de JUAN
                           habria sido la almohadilla.

    Entre las admisibles gana la de mayor volumen.
    """
    for f in filas:
        f["admisible"] = (f["fracCort"] >= FRACCION_CORTICAL_MINIMA
                          and f["compac"] >= COMPACIDAD_MINIMA_REFERENCIA)
    return [f for f in filas if f["admisible"]]


def _elegirReferencia(filas, avisar):
    """Devuelve (referencia, relajado). 'relajado' dice si hubo que soltar
    algun filtro, y en ese caso hay que mirar la eleccion con desconfianza."""
    cand = candidatasAReferencia(filas)
    if cand:
        return max(cand, key=lambda f: f["vol"]), None

    # Salvaguarda 1 (ver [C4]): ningun candidato compacto. Puede pasar si la
    # boveda esta partida en dos cascaras finas y ninguna sola es compacta.
    porCortical = [f for f in filas if f["fracCort"] >= FRACCION_CORTICAL_MINIMA]
    if porCortical:
        avisar("Ninguna pieza llega a compacidad %.0f. Se relaja ese filtro y "
               "se elige por volumen entre las de cortical suficiente. REVISA "
               "cual quedo como pieza principal." % COMPACIDAD_MINIMA_REFERENCIA)
        for f in porCortical:
            f["admisible"] = True
        return max(porCortical, key=lambda f: f["vol"]), "compacidad"

    # Salvaguarda 2: hueso muy poco mineralizado o serie de kernel muy blando.
    if filas:
        avisar("Ninguna pieza llega a cortical %.0f%%. Se relajan TODOS los "
               "filtros. Es muy probable que la serie elegida no sirva "
               "(kernel blando?). REVISA todo."
               % (FRACCION_CORTICAL_MINIMA * 100.0))
        for f in filas:
            f["admisible"] = True
        return max(filas, key=lambda f: f["vol"]), "cortical"
    return None, None


def decidir(labels, filas, espaciadoZYX, referenciaForzada=None,
            modoPostop=False, avisar=print):
    """
    Asigna veredicto a cada isla.

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
            avisar("No existe la pieza %s." % referenciaForzada)
            return None, [], [], None
        candidatasAReferencia(filas)          # deja marcada la admisibilidad
        referencia = filas[referenciaForzada - 1]
        avisar("Pieza principal FORZADA a mano -> pieza %d (%.2f cm3)."
               % (referenciaForzada, referencia["vol"]))
    else:
        referencia, relajado = _elegirReferencia(filas, avisar=avisar)
        if referencia is None:
            return None, [], [], None
    referencia["referencia"] = True

    distAceptacion = (DIST_POSTOP_ACEPTACION_MM if modoPostop
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
            f["categoria"] = "pieza principal"
            aceptadas.append(f)
            continue

        if max(f["extMM"]) > limiteExt:
            # Un fragmento de craneo NO puede ser mas grande que el craneo del
            # que salio. Limite RELATIVO a la referencia: escala con el
            # paciente, asi que no se rompe con la edad.
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "EXTENSION"
            f["motivo"] = ("mide %.0f mm, mas de %.2fx la pieza principal (%.0f mm)"
                           % (max(f["extMM"]), FACTOR_EXTENSION_MAXIMA, limiteExt))
            if f["vol"] >= 5.0 and f["fracCort"] >= 0.30:
                avisar("Se descarto por tamano una pieza de %.1f cm3 con %.0f%% "
                       "de cortical. Si eso era hueso, la PIEZA PRINCIPAL esta "
                       "mal elegida: cambiala desde la lista."
                       % (f["vol"], f["fracCort"] * 100.0))
            continue

        if (not modoPostop) and f["compac"] < COMPACIDAD_MINIMA_PIEZA:
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
            f["motivo"] = ("camilla o soporte: lado corto %.1f mm y %.0f%% del campo"
                           % (min(f["extMM"]), max(f["fracFOV"]) * 100.0))
            continue

        if f["fracCort"] < FRACCION_CORTICAL_MINIMA:
            # Ver [3]: el filtro de material no es borrado duro. Lo que esta
            # pegado al craneo no se borra, se manda a revision.
            if f["distRef"] is not None and f["distRef"] <= distAceptacion:
                f["veredicto"] = "DUDOSA"
                f["categoria"] = "MATERIAL_PEGADO"
                f["motivo"] = ("poco denso (%.1f%% de cortical) pero pegado al "
                               "craneo, a %.1f mm"
                               % (f["fracCort"] * 100.0, f["distRef"]))
            else:
                f["veredicto"] = "ELIMINADA"
                f["categoria"] = "MATERIAL"
                f["motivo"] = ("poco denso (%.1f%% < %.0f%%) y lejos, a %.0f mm"
                               % (f["fracCort"] * 100.0,
                                  FRACCION_CORTICAL_MINIMA * 100.0, f["distRef"]))
            continue

        if f["distRef"] is not None and f["distRef"] <= distAceptacion:
            f["veredicto"] = "ACEPTADA"
            f["categoria"] = "pegada al craneo"
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
            f["motivo"] = ("separada del craneo por %.1f mm (el limite para "
                           "aceptarla sola es %.1f mm)"
                           % (f["distRef"], distAceptacion))
            dudosas.append(f)
        else:
            f["veredicto"] = "ELIMINADA"
            f["categoria"] = "LEJOS"
            f["motivo"] = "a %.0f mm del craneo" % f["distMasa"]

    aceptadas.sort(key=lambda f: f["vol"], reverse=True)
    dudosas.sort(key=lambda f: f["vol"], reverse=True)

    totalAceptado = sum(f["vol"] for f in aceptadas)
    if totalAceptado < VOLUMEN_CRANEO_MINIMO_CM3:
        avisar("Lo aceptado suma %.1f cm3, menos que el minimo anatomico de "
               "%.1f cm3. La pieza principal puede estar mal elegida, o la "
               "serie no sirve." % (totalAceptado, VOLUMEN_CRANEO_MINIMO_CM3))

    if referencia["tocaBorde"]:
        avisar("NOTA DE CALIDAD: el craneo llega al borde del campo de vision "
               "de la tomografia. Puede venir RECORTADO y el modelo quedar "
               "incompleto. Esto no cambia la segmentacion: es un dato para el "
               "cirujano.")

    return referencia, aceptadas, dudosas, relajado


def tablaDeCandidatas(filas, referencia, cuantas=3):
    """Las mejores candidatas a pieza principal, para que la decision sea
    visible y no solo su resultado."""
    orden = sorted(filas, key=lambda f: f["vol"], reverse=True)[:max(cuantas, 3)]
    lineas = ["  Candidatas a pieza principal (por volumen):",
              "    #   vol(cm3)   ext(mm)  cort%   compac   estado"]
    for i, f in enumerate(orden, start=1):
        if referencia is not None and f is referencia:
            estado = "<= ELEGIDA"
        elif not f["admisible"]:
            motivos = []
            if f["fracCort"] < FRACCION_CORTICAL_MINIMA:
                motivos.append("poco densa")
            if f["compac"] < COMPACIDAD_MINIMA_REFERENCIA:
                motivos.append("alargada")
            estado = "excluida: " + " y ".join(motivos)
        else:
            estado = "admisible"
        lineas.append("   %2d  %9.2f  %8.1f  %5.1f  %7.1f   %s"
                      % (i, f["vol"], max(f["extMM"]), f["fracCort"] * 100.0,
                         f["compac"], estado))
    return "\n".join(lineas)

# ============================================================
# FIN DEL NUCLEO DE DECISION
# ============================================================


# ------------------------------------------------------------
# Utilidades de escena
# ------------------------------------------------------------
def _borrarNodosPorNombre(nombre):
    """Con nombres fijos, un duplicado rompe el Bloque F en silencio: si queda
    un 'Craneo_Automatico' viejo, Slicer bautiza al nuevo
    'Craneo_Automatico_1' y la busqueda por nombre agarra el VIEJO."""
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


def _contarRegiones(polyData):
    if polyData is None or polyData.GetNumberOfPoints() == 0:
        return 0
    c = vtk.vtkPolyDataConnectivityFilter()
    c.SetInputData(polyData)
    c.SetExtractionModeToAllRegions()
    c.Update()
    return int(c.GetNumberOfExtractedRegions())


def _limpiarRuidoMalla(polyData, minimoPuntos, log=print):
    """Umbral ABSOLUTO en cantidad de puntos.

    Con un umbral RELATIVO, cualquier placa craneal desconectada por una
    sutura abierta que mida menos del 10 % de la boveda se borra sola: por eso
    'Craneo_Final' aparecia sin mandibula en el caso 403."""
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
    log("Bloque A: limpieza de malla: %d region(es) -> %d conservada(s) "
        "(umbral ABSOLUTO %d puntos)." % (nRegiones, len(conservar), minimoPuntos))
    if descartadas and descartadas[0] >= minimoPuntos * 0.5:
        log("Bloque A:   ATENCION - la mayor region descartada esta cerca del "
            "umbral. Si falta hueso, bajar MINIMO_PUNTOS_MALLA.")
    return salida


def _decimarMalla(polyData, reduccion, log=print):
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
        log("Bloque A: la decimacion al %.0f%% perdio regiones; reintento."
            % (intento * 100))
        intento /= 2.0
    log("Bloque A: se devuelve la malla SIN decimar para no perder piezas.")
    return polyData


# ============================================================
# LA CLASE QUE USA EL MODULO
# ============================================================
class BloqueA:
    """
    Estado de UNA corrida del Bloque A sobre UN paciente.

    Por que una clase y no funciones sueltas con globals (que es como estaba
    en el script de consola): el modulo puede cerrar la escena, cargar otro
    estudio y volver a correr. Con globals, los restos del paciente anterior
    sobreviven a eso y aparecen mezclados. Con una instancia, se tira y se
    hace una nueva.
    """

    def __init__(self, log=print):
        self.log = log
        self.volumeNode = None
        self.labels = None
        self.espaciadoZYX = None
        self.filas = []
        self.piezas = []          # lo que ve la interfaz
        self.segNode = None
        self.colores = {}
        self.modelo = None
        self.referencia = None
        self.modoPostop = False
        self.avisos = []

    # --------------------------------------------------------
    def _avisar(self, mensaje):
        self.avisos.append(mensaje)
        self.log("Bloque A: " + str(mensaje))

    def _limpiarEscenaPrevia(self):
        for atributo in ("segNode", "modelo"):
            nodo = getattr(self, atributo, None)
            if nodo is not None:
                try:
                    slicer.mrmlScene.RemoveNode(nodo)
                except Exception:
                    pass
            setattr(self, atributo, None)
        _borrarNodosPorNombre(NOMBRE_NODO_SEGMENTACION)
        _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)
        self.colores = {}
        self.piezas = []

    # --------------------------------------------------------
    # PASO PRINCIPAL
    # --------------------------------------------------------
    def generar(self, volumeNode, modoPostop=False):
        """
        Segmenta, decide y deja las piezas listas para la revision.
        NO fusiona nada: eso lo hace confirmarCraneo().

        Devuelve la lista de piezas (dicts listos para la tabla) o None.
        """
        if volumeNode is None:
            self.log("Bloque A: no hay volumen cargado.")
            return None

        self.avisos = []
        self.volumeNode = volumeNode
        self.modoPostop = bool(modoPostop)

        self.log("=" * 66)
        self.log("Bloque A %s%s" % (BLOQUE_A_VERSION,
                                    "  *** MODO POSTOPERATORIO ***"
                                    if self.modoPostop else ""))
        dims = volumeNode.GetImageData().GetDimensions()
        esp = volumeNode.GetSpacing()
        self.log("Bloque A: volumen  %s" % volumeNode.GetName())
        self.log("Bloque A: %d x %d x %d voxels, %.3f x %.3f x %.3f mm"
                 % (dims + esp))
        nombre = volumeNode.GetName().lower()
        if ("cerebro" in nombre or "brain" in nombre) and "hueso" not in nombre:
            self._avisar("La serie cargada es de ventana BLANDA (cerebro). El "
                         "hueso delgado pierde densidad. Si el estudio tiene "
                         "serie de hueso, conviene usar esa.")

        self._limpiarEscenaPrevia()

        volArr = slicer.util.arrayFromVolume(volumeNode)
        self.espaciadoZYX = (esp[2], esp[1], esp[0])

        labels, filas, info = analizarVolumen(volArr, self.espaciadoZYX)
        if labels is None or not filas:
            self.log("Bloque A: no se encontro ninguna pieza util entre %d y "
                     "%d HU." % (HU_MIN, HU_MAX))
            return None

        self.log("Bloque A: %d region(es) detectada(s), %d por encima de "
                 "%.2f cm3." % (info["nIslas"], len(filas), VOL_MINIMO_CM3))
        self.log("Bloque A: total de hueso en la mascara: %.2f cm3" % info["totalCM3"])

        self.labels = labels
        self.filas = filas

        referencia, aceptadas, dudosas, relajado = decidir(
            labels, filas, self.espaciadoZYX, modoPostop=self.modoPostop,
            avisar=self._avisar)
        if referencia is None:
            self.log("Bloque A: ninguna pieza puede ser la principal. Revisa "
                     "la serie elegida.")
            return None

        self.log("")
        self.log(tablaDeCandidatas(filas, referencia))
        if relajado:
            self._avisar("Se relajo el filtro de %s para poder elegir la pieza "
                         "principal. Mirala con desconfianza." % relajado)

        return self._construirEscena(referencia, aceptadas, dudosas)

    def redecidir(self, referenciaForzada=None, modoPostop=None):
        """Rehace la decision con el analisis ya en memoria (la parte lenta,
        etiquetar el volumen, no se repite)."""
        if self.labels is None or not self.filas or self.volumeNode is None:
            self.log("Bloque A: no hay un analisis previo en memoria.")
            return None
        if modoPostop is not None:
            self.modoPostop = bool(modoPostop)
        self.avisos = []
        referencia, aceptadas, dudosas, _ = decidir(
            self.labels, self.filas, self.espaciadoZYX,
            referenciaForzada=referenciaForzada, modoPostop=self.modoPostop,
            avisar=self._avisar)
        if referencia is None:
            return None
        self.log("")
        self.log(tablaDeCandidatas(self.filas, referencia))
        return self._construirEscena(referencia, aceptadas, dudosas)

    # --------------------------------------------------------
    def _construirEscena(self, referencia, aceptadas, dudosas):
        self._limpiarEscenaPrevia()
        segNode = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLSegmentationNode', NOMBRE_NODO_SEGMENTACION)
        segNode.CreateDefaultDisplayNodes()
        segNode.SetReferenceImageGeometryParameterFromVolumeNode(self.volumeNode)
        segmentacion = segNode.GetSegmentation()

        buf = np.zeros(self.labels.shape, dtype=np.uint8)
        lista = []
        for i, f in enumerate(aceptadas + dudosas, start=1):
            etiquetaZona = "OK" if f["veredicto"] == "ACEPTADA" else "DUDOSA"
            segId = segmentacion.AddEmptySegment(
                "", "%02d_%s_%.2fcm3" % (i, etiquetaZona, f["vol"]))
            caja = f["caja"]
            buf[caja] = (self.labels[caja] == f["etiqueta"]).astype(np.uint8)
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                buf, segNode, segId, self.volumeNode)
            buf[caja] = 0
            color = COLOR_ACEPTADA if f["veredicto"] == "ACEPTADA" else COLOR_DUDOSA
            segmentacion.GetSegment(segId).SetColor(*color)
            self.colores[segId] = color
            lista.append({
                "numero": i, "segId": segId, "vol": f["vol"],
                "distRef": f["distRef"], "fracCort": f["fracCort"],
                "compac": f["compac"], "zona": f["veredicto"],
                "categoria": f["categoria"], "nota": f["motivo"],
                "referencia": f["referencia"],
            })
        del buf

        _forzarSuperficie(segNode)
        self.segNode = segNode
        self.piezas = lista
        self.referencia = referencia

        self.log("")
        self.log("--- PIEZAS EN LA ESCENA ---")
        for d in lista:
            marca = "  *pieza principal*" if d["referencia"] else ""
            extra = ("   <- " + d["nota"]) if d["nota"] else ""
            self.log("  [%2d] %-8s vol=%8.2f cm3  cort=%5.1f%%  dist=%5.1f mm%s%s"
                     % (d["numero"], d["zona"], d["vol"], d["fracCort"] * 100.0,
                        d["distRef"] if d["distRef"] is not None else -1,
                        marca, extra))
        nd = sum(1 for d in lista if d["zona"] == "DUDOSA")
        self.log("Bloque A: %d aceptada(s) + %d dudosa(s)."
                 % (len(lista) - nd, nd))
        return lista

    # --------------------------------------------------------
    # REVISION
    # --------------------------------------------------------
    def resaltar(self, numero):
        """Pinta esa pieza de rojo; el resto vuelve a su color."""
        if self.segNode is None:
            return
        segmentacion = self.segNode.GetSegmentation()
        for d in self.piezas:
            seg = segmentacion.GetSegment(d["segId"])
            if seg is None:
                continue
            if d["numero"] == numero:
                seg.SetColor(*COLOR_RESALTE)
            else:
                seg.SetColor(*self.colores.get(d["segId"], (0.6, 0.6, 0.6)))

    def mostrarTodas(self):
        if self.segNode is None:
            return
        disp = self.segNode.GetDisplayNode()
        segmentacion = self.segNode.GetSegmentation()
        for d in self.piezas:
            disp.SetSegmentVisibility(d["segId"], True)
            seg = segmentacion.GetSegment(d["segId"])
            if seg is not None:
                seg.SetColor(*self.colores.get(d["segId"], (0.6, 0.6, 0.6)))

    def soloDudosas(self):
        """Deja a la vista unicamente lo que el algoritmo no supo decidir."""
        if self.segNode is None:
            return 0
        disp = self.segNode.GetDisplayNode()
        n = 0
        for d in self.piezas:
            v = d["zona"] == "DUDOSA"
            disp.SetSegmentVisibility(d["segId"], v)
            n += 1 if v else 0
        return n

    def eliminar(self, numero):
        if self.segNode is None:
            return False
        d = next((x for x in self.piezas if x["numero"] == numero), None)
        if d is None:
            return False
        if d["referencia"]:
            self.log("Bloque A: esa es la pieza principal del craneo. Si de "
                     "verdad hay que cambiarla, se elige otra como principal "
                     "en vez de borrarla.")
            return False
        self.segNode.GetSegmentation().RemoveSegment(d["segId"])
        self.piezas = [x for x in self.piezas if x["numero"] != numero]
        return True

    def eliminarDudosas(self):
        """Saca de una sola vez todas las marcadas como dudosas. En los 9
        casos relevados casi todas eran vertebras cervicales."""
        dudosas = [d["numero"] for d in self.piezas if d["zona"] == "DUDOSA"]
        for n in sorted(dudosas, reverse=True):
            self.eliminar(n)
        return len(dudosas)

    # --------------------------------------------------------
    # CONFIRMACION Y EXPORTACION
    # --------------------------------------------------------
    def confirmarCraneo(self):
        """
        Fusiona las piezas que quedaron en un unico segmento con el nombre que
        espera el Bloque F.

        Se hace con un OR de mascaras en numpy y UNA sola escritura, no con el
        efecto 'Logical operators' del Segment Editor: ese efecto arrastra las
        reglas de masking del editor (estado oculto) y dejaba viva la
        representacion de superficie vieja de cada segmento.
        """
        if self.segNode is None or self.volumeNode is None:
            self.log("Bloque A: primero hay que generar el craneo.")
            return None
        if not self.piezas:
            self.log("Bloque A: no queda ninguna pieza.")
            return None

        segmentacion = self.segNode.GetSegmentation()
        esp = self.volumeNode.GetSpacing()
        volVoxCM3 = (esp[0] * esp[1] * esp[2]) / 1000.0

        union = None
        sumaPiezas = 0.0
        leidas = 0
        ids = [d["segId"] for d in self.piezas
               if segmentacion.GetSegment(d["segId"]) is not None]
        for segId in ids:
            arr = slicer.util.arrayFromSegmentBinaryLabelmap(
                self.segNode, segId, self.volumeNode)
            if arr is None:
                self.log("Bloque A: ADVERTENCIA - no pude leer %s, queda afuera."
                         % segId)
                continue
            arr = arr.astype(bool)
            sumaPiezas += int(np.count_nonzero(arr)) * volVoxCM3
            union = arr.copy() if union is None else np.logical_or(union, arr)
            leidas += 1

        if union is None or not union.any():
            self.log("Bloque A: la union quedo vacia.")
            return None

        volUnion = int(np.count_nonzero(union)) * volVoxCM3
        base = ids[0]
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            union.astype(np.uint8), self.segNode, base, self.volumeNode)
        for i in range(segmentacion.GetNumberOfSegments() - 1, -1, -1):
            segId = segmentacion.GetNthSegmentID(i)
            if segId != base:
                segmentacion.RemoveSegment(segId)
        segmentacion.GetSegment(base).SetName(NOMBRE_SEGMENTO_CRANEO)
        segmentacion.GetSegment(base).SetColor(*COLOR_ACEPTADA)
        _forzarSuperficie(self.segNode)

        estructura = ndimage.generate_binary_structure(3, 3)
        _, piezasVoxel = ndimage.label(union, structure=estructura)
        poly = vtk.vtkPolyData()
        self.segNode.GetClosedSurfaceRepresentation(base, poly)
        regiones = _contarRegiones(poly)

        self.log("")
        self.log("--- CONFIRMACION ---")
        self.log("  piezas fusionadas       : %d de %d" % (leidas, len(ids)))
        self.log("  suma de las piezas      : %.2f cm3" % sumaPiezas)
        self.log("  union (sin solapamiento): %.2f cm3" % volUnion)
        self.log("  piezas conexas (voxels) : %d" % piezasVoxel)
        self.log("  regiones en la malla    : %d" % regiones)
        if regiones < piezasVoxel:
            self._avisar("La malla tiene menos regiones que los voxeles: alguna "
                         "placa fina se disolvio en el suavizado. Bajar "
                         "SUAVIZADO_SUPERFICIE a '0.1' y repetir.")

        return {"volumenCM3": volUnion, "piezasFusionadas": leidas,
                "piezasConexas": piezasVoxel, "regionesMalla": regiones}

    def enviarAPlanner(self):
        """Exporta el segmento a un Model node con el nombre que espera el
        Bloque F."""
        if self.segNode is None:
            self.log("Bloque A: primero generar y confirmar el craneo.")
            return None
        segmentacion = self.segNode.GetSegmentation()
        segId = segmentacion.GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_CRANEO)
        if not segId:
            self.log("Bloque A: todavia no existe el segmento '%s'."
                     % NOMBRE_SEGMENTO_CRANEO)
            return None

        _forzarSuperficie(self.segNode)
        bruto = vtk.vtkPolyData()
        self.segNode.GetClosedSurfaceRepresentation(segId, bruto)
        if bruto.GetNumberOfPoints() == 0:
            self.log("Bloque A: la conversion a superficie no devolvio geometria.")
            return None

        regBruto, celBruto = _contarRegiones(bruto), bruto.GetNumberOfCells()
        malla = _limpiarRuidoMalla(bruto, MINIMO_PUNTOS_MALLA, log=self.log)
        regLimpio, celLimpio = _contarRegiones(malla), malla.GetNumberOfCells()
        malla = _decimarMalla(malla, REDUCCION_MALLA, log=self.log)

        if self.modelo is not None:
            try:
                slicer.mrmlScene.RemoveNode(self.modelo)
            except Exception:
                pass
        _borrarNodosPorNombre(NOMBRE_MODELO_CRANEO)

        modelo = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLModelNode', NOMBRE_MODELO_CRANEO)
        modelo.SetAndObservePolyData(malla)
        modelo.CreateDefaultDisplayNodes()
        modelo.GetDisplayNode().SetColor(*COLOR_ACEPTADA)
        modelo.GetDisplayNode().SetScalarVisibility(False)
        self.modelo = modelo

        display = self.segNode.GetDisplayNode()
        if display is not None:
            # A partir de aca se trabaja sobre el MODELO. Si quedan los dos
            # visibles, Slicer renderiza el craneo dos veces en cada
            # movimiento del mouse: es una de las causas del lag al rotar.
            display.SetVisibility(False)

        self.log("")
        self.log("--- MODELO 3D ---")
        self.log("  bruto  : %7d triangulos, %3d region(es)" % (celBruto, regBruto))
        self.log("  limpio : %7d triangulos, %3d region(es)" % (celLimpio, regLimpio))
        self.log("  final  : %7d triangulos, %3d region(es)"
                 % (malla.GetNumberOfCells(), _contarRegiones(malla)))

        if modelo.GetName() != NOMBRE_MODELO_CRANEO:
            self._avisar("El modelo quedo como '%s' porque habia otro nodo con "
                         "ese nombre. El corte va a usar el equivocado."
                         % modelo.GetName())
        return modelo

    # --------------------------------------------------------
    def verificarCompatibilidad(self):
        """Chequea que la escena tenga EXACTAMENTE lo que el Bloque F busca.
        Devuelve (ok, lineas)."""
        lineas = []
        ok = True

        segs = [n for n in slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
                if n.GetName() == NOMBRE_NODO_SEGMENTACION]
        if len(segs) == 1:
            lineas.append("[OK]    segmentacion '%s'" % NOMBRE_NODO_SEGMENTACION)
            sid = segs[0].GetSegmentation().GetSegmentIdBySegmentName(
                NOMBRE_SEGMENTO_CRANEO)
            if sid:
                lineas.append("[OK]    segmento '%s'" % NOMBRE_SEGMENTO_CRANEO)
            else:
                lineas.append("[FALTA] segmento '%s': confirmar el craneo"
                              % NOMBRE_SEGMENTO_CRANEO)
                ok = False
        elif not segs:
            lineas.append("[FALTA] segmentacion: generar el craneo 3D")
            ok = False
        else:
            lineas.append("[ERROR] hay %d nodos llamados '%s'. Borrar los sobrantes."
                          % (len(segs), NOMBRE_NODO_SEGMENTACION))
            ok = False

        mods = [n for n in slicer.util.getNodesByClass("vtkMRMLModelNode")
                if n.GetName() == NOMBRE_MODELO_CRANEO]
        if len(mods) == 1:
            lineas.append("[OK]    modelo 3D '%s'" % NOMBRE_MODELO_CRANEO)
        elif not mods:
            lineas.append("[FALTA] modelo 3D: confirmar el craneo")
            ok = False
        else:
            lineas.append("[ERROR] hay %d modelos llamados '%s'."
                          % (len(mods), NOMBRE_MODELO_CRANEO))
            ok = False

        return ok, lineas
