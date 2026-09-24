# ============================================================
# CRANIOPLAN - BLOQUE A v7.1  ::  MODO DIAGNOSTICO
#   script de consola para 3D Slicer 5.10 - NO BORRA NADA
#   incluye el etiquetado manual (antes era un archivo aparte)
# ============================================================
#
# PARA QUE SIRVE
# --------------
# Relevar cada tomografia nueva ANTES de correr el Bloque A v7.1, para que
# los umbrales sigan saliendo de datos y no de intuiciones. Reemplaza a
# BloqueA_v6_DIAGNOSTICO.py + su complemento: ahora es un solo pegado.
#
# QUE HACE
# --------
# 1. Calcula las mismas metricas que el v7.1 y aplica el MISMO veredicto,
#    pero NO elimina nada: crea todas las islas como segmentos y las pinta
#    segun lo que el v7.1 habria hecho con ellas.
# 2. Mide la ENVOLVENCIA (env%) de todas las islas, incluidas las que el
#    v7.1 elimina. Es la metrica candidata para detectar almohadillas de
#    inmovilizacion y hoy NO decide nada: se releva para calibrarla.
# 3. Deja registrar a mano que es cada isla y exporta todo a un CSV que se
#    ACUMULA entre pacientes.
#
# ------------------------------------------------------------
# QUE CAMBIO RESPECTO DEL DIAGNOSTICO ANTERIOR (v6-DIAG)
# ------------------------------------------------------------
# [A] LA REFERENCIA SE ELIGE CON LA LOGICA DEL v7.1, no con la del v6.
#     Esto importa mucho mas de lo que parece: en la TC 175 (RUIZ) el v6
#     eliminaba la calota por tocar el borde, la referencia terminaba
#     siendo una vertebra de 2.38 cm3, y TODA la columna de distancias de
#     ese caso quedo medida contra la vertebra en vez de contra el craneo.
#     Por eso las mandibulas figuraban a 19-20 mm y hubo que excluir el
#     caso entero de la calibracion. Con la referencia correcta eso no
#     vuelve a pasar.
#
# [B] DOS COLUMNAS DE DISTANCIA en vez de una, porque son cosas distintas:
#       dRef  = distancia a la isla de REFERENCIA  -> decide ACEPTADA
#       dMasa = distancia a la masa ya aceptada    -> decide DUDOSA
#     En el v6 estaban mezcladas en un solo numero.
#
# [C] COLUMNA env% nueva (ver arriba).
#
# [D] EL VEREDICTO QUE SE COMPARA ES EL DEL v7.1, asi que discrepancias()
#     mide los errores del algoritmo que vamos a usar, no del anterior.
#
# ------------------------------------------------------------
# FLUJO POR CASO
# ------------------------------------------------------------
#   listar_volumenes()
#   usar_volumen("Hueso")        # elegi la serie y anotala
#   diagnostico()
#   barrido_hu()
#   crudo()                      # opcional: mascara sin ningun filtro
#   ir_a(3) ; solo(3) ; et(3, "mandibula")      # etiquetar isla por isla
#   pendientes()                 # que falta
#   reporte()                    # bloque compacto para pasar
#   csv()                        # exporta y ACUMULA
#
# ------------------------------------------------------------
# ETIQUETAS SUGERIDAS (sin acentos, para poder agrupar)
# ------------------------------------------------------------
#   hueso_craneal  base_craneo  mandibula  maxilar  diente
#   vertebra  hueso_otro
#   camilla  almohadilla  banda  electrodo  chupete  tubo  metal  ruido
#   no_se      <- usarla sin culpa, es informacion valida
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS
# ------------------------------------------------------------
# - Crea como maximo MAX_SEGMENTOS islas; el resto se agrupa en
#   "99_Resto_pequenas". Con 793 islas, crearlas todas cuelga Slicer.
# - Las distancias se miden superficie-superficie con submuestreo. Error
#   tipico < 1 mm: sirve para decidir "pegada" vs "lejos", no para medir
#   el ancho de una sutura.
# - env% depende de donde este el centro de la referencia. Si la
#   referencia sale mal elegida, env% no significa nada. Mirá primero que
#   la linea "referencia ->" sea el craneo.
# ============================================================


# ============================================================
# CONFIG
# ============================================================

NOMBRE_VOLUMEN = None          # None = ultimo volumen cargado

HU_MIN = 300
HU_MAX = 3000

# --- umbrales del v7.1, replicados SOLO para etiquetar el veredicto ---
FRACCION_CORTICAL_MINIMA = 0.05
UMBRAL_CORTICAL_HU       = 700.0
DIST_ACEPTACION_MM       = 2.0     # contra la REFERENCIA
DIST_DUDOSA_MM           = 15.0    # contra la masa aceptada
ESPESOR_LOSA_MM          = 30.0
FRACCION_LOSA_FOV        = 0.85
MAX_EXTENSION_CRANEO_MM  = 180.0
VOLUMEN_CRANEO_MINIMO_CM3 = 20.0

# Envolvencia: bins de area igual, uniformes en cos(theta) y phi.
ENV_BINS_COSTHETA = 12
ENV_BINS_PHI      = 24

# Piso para que una isla se liste. Mas bajo que el del v7.1 a proposito:
# queremos VER lo que el v7.1 descarta sin analizar.
PISO_LISTADO_CM3 = 0.05

MAX_SEGMENTOS = 30
MAX_PUNTOS_KDTREE = 40000

BARRIDO_HU = [150, 200, 250, 300, 350, 400]

NOMBRE_NODO_DIAG  = "DIAG_Islas"
NOMBRE_NODO_CRUDO = "DIAG_Crudo"
SUAVIZADO_SUPERFICIE = "0.3"

# CSV acumulativo. Esquema nuevo (agrega env, dRef, dMasa), por eso el
# nombre cambia: no se puede mezclar con el CSV del v6-DIAG, que tenia
# menos columnas. Los 143 registros viejos siguen sirviendo aparte.
RUTA_CSV = None   # None = ~/Documents/CranioPlan_diag_islas_v71.csv


# ============================================================
# CODIGO
# ============================================================

import os
import numpy as np
import vtk
import slicer
from scipy import ndimage
from scipy.spatial import cKDTree

CRANIOPLAN_DIAG = "v7.1-DIAG (2026-08-08) consola - no destructivo"

_d = {
    "volumeNode": None, "labels": None, "espaciadoZYX": None,
    "filas": [], "creadas": [], "segNode": None, "crudoNode": None,
    "colores": {}, "huMin": None, "totalCM3": None, "nIslas": None,
}

COLOR = {
    "ACEPTADA":  (0.90, 0.80, 0.60),   # hueso
    "DUDOSA":    (1.00, 0.60, 0.10),   # naranja
    "MATERIAL":  (0.25, 0.75, 0.70),   # turquesa: cortical insuficiente
    "LOSA":      (0.25, 0.45, 0.95),   # azul: tabla de camilla
    "LEJOS":     (0.70, 0.35, 0.90),   # violeta
    "CHICA":     (0.55, 0.55, 0.55),   # gris
}

_ETIQUETAS = ["hueso_craneal", "base_craneo", "mandibula", "maxilar", "diente",
              "vertebra", "hueso_otro", "camilla", "almohadilla", "banda",
              "electrodo", "chupete", "tubo", "metal", "ruido", "no_se"]
_CRANEO = {"hueso_craneal", "base_craneo", "mandibula", "maxilar", "diente"}
_BASURA = {"camilla", "almohadilla", "banda", "electrodo", "chupete",
           "tubo", "metal", "ruido"}


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


def _contarRegiones(polyData):
    if polyData is None or polyData.GetNumberOfPoints() == 0:
        return 0
    c = vtk.vtkPolyDataConnectivityFilter()
    c.SetInputData(polyData)
    c.SetExtractionModeToAllRegions()
    c.Update()
    return int(c.GetNumberOfExtractedRegions())


# ------------------------------------------------------------
# Seleccion de volumen
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


def _volumen():
    v = _d.get("volumeNode")
    if v is None or slicer.mrmlScene.GetNodeByID(v.GetID()) is None:
        v = usar_volumen(NOMBRE_VOLUMEN)
    return v


def _encabezado(volumeNode):
    dims = volumeNode.GetImageData().GetDimensions()
    esp = volumeNode.GetSpacing()
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
# Analisis
# ------------------------------------------------------------

def _analizar(volArr, espaciadoZYX, huMin, huMax):
    mascara = (volArr >= huMin) & (volArr <= huMax)
    estructura = ndimage.generate_binary_structure(3, 3)
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
        if vol < PISO_LISTADO_CM3:
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
            "tocaBorde": k in etiquetasBorde,
            "esLosa": _esLosaDeSoporte(extMM, fracFOV),
            "ijk": (idx[:, 2].mean() + caja[2].start,
                    idx[:, 1].mean() + caja[1].start,
                    idx[:, 0].mean() + caja[0].start),
            "distRef": None, "distMasa": None, "env": None,
            "veredicto": None, "motivo": "", "referencia": False,
            "etiqueta_real": "", "nota_real": "",
        })
    filas.sort(key=lambda f: f["vol"], reverse=True)
    return labels, filas, {"nIslas": nIslas, "valorRelleno": valorRelleno,
                           "totalCM3": totalCM3}


def _veredictoV71(labels, filas, espaciadoZYX):
    """Aplica EXACTAMENTE la logica del v7.1, pero solo para etiquetar.
    No elimina nada."""
    for f in filas:
        f.update({"distRef": None, "distMasa": None, "env": None,
                  "veredicto": None, "motivo": "", "referencia": False})

    oseas = []
    for f in filas:
        if f["fracCort"] < FRACCION_CORTICAL_MINIMA:
            f["veredicto"] = "MATERIAL"
            f["motivo"] = ("cortical %.1f%% < %.0f%%"
                           % (f["fracCort"] * 100.0, FRACCION_CORTICAL_MINIMA * 100.0))
        elif f["esLosa"]:
            f["veredicto"] = "LOSA"
            f["motivo"] = ("lado corto %.1f mm y %.0f%% del campo"
                           % (min(f["extMM"]), max(f["fracFOV"]) * 100.0))
        else:
            oseas.append(f)
    if not oseas:
        return None

    cache = {}
    for f in oseas:
        k = f["etiqueta"]
        cache[k] = None
        pts = _puntosDeSuperficie(labels, k, f["caja"], espaciadoZYX,
                                  MAX_PUNTOS_KDTREE)
        cache[k] = (pts, cKDTree(pts))
    D = _matrizDeDistancias(oseas, cache)

    def puedeSerReferencia(f, exentarBorde):
        if f["esLosa"] or f["fracCort"] < FRACCION_CORTICAL_MINIMA:
            return False
        if max(f["extMM"]) > MAX_EXTENSION_CRANEO_MM:
            return False
        if not exentarBorde and f["tocaBorde"]:
            return False
        return True

    def aceptarDesde(ref):
        kRef = ref["etiqueta"]
        ac = {kRef}
        ref["distRef"] = 0.0
        for f in oseas:
            k = f["etiqueta"]
            if k == kRef:
                continue
            f["distRef"] = D[k][kRef]
            if f["distRef"] <= DIST_ACEPTACION_MM:
                ac.add(k)
        for f in oseas:
            f["distMasa"] = min(D[f["etiqueta"]][a] for a in ac)
        return ac

    referencia, aceptadas, total = None, set(), 0.0
    plausibles = [f for f in oseas if puedeSerReferencia(f, False)]
    if plausibles:
        referencia = max(plausibles, key=lambda f: f["vol"])
        aceptadas = aceptarDesde(referencia)
        total = sum(f["vol"] for f in oseas if f["etiqueta"] in aceptadas)
        print("DIAG: pasada 1 (referencia sin tocar el borde): referencia "
              "%.2f cm3, %.1f cm3 aceptados." % (referencia["vol"], total))
    else:
        print("DIAG: pasada 1 sin candidata a referencia.")

    if referencia is None or total < VOLUMEN_CRANEO_MINIMO_CM3:
        print("DIAG: no alcanza para ser un craneo (minimo %.1f cm3) -> modo "
              "RESCATE, la referencia queda exenta del borde."
              % VOLUMEN_CRANEO_MINIMO_CM3)
        p2 = [f for f in oseas if puedeSerReferencia(f, True)]
        if p2:
            ref2 = max(p2, key=lambda f: f["vol"])
            ac2 = aceptarDesde(ref2)
            total2 = sum(f["vol"] for f in oseas if f["etiqueta"] in ac2)
            if total2 > total:
                referencia, aceptadas, total = ref2, ac2, total2
                print("DIAG: rescate aplicado -> referencia %.2f cm3, %.1f cm3."
                      % (referencia["vol"], total))
            elif referencia is not None:
                aceptadas = aceptarDesde(referencia)

    if referencia is None:
        return None
    referencia["referencia"] = True

    centroMM = np.mean(cache[referencia["etiqueta"]][0], axis=0)
    for f in oseas:
        f["env"] = _envolvencia(cache[f["etiqueta"]][0], centroMM)

    for f in oseas:
        if f["etiqueta"] in aceptadas:
            f["veredicto"] = "ACEPTADA"
        elif f["distMasa"] is not None and f["distMasa"] <= DIST_DUDOSA_MM:
            f["veredicto"] = "DUDOSA"
            f["motivo"] = ("dRef %.1f mm > %.1f mm"
                           % (f["distRef"], DIST_ACEPTACION_MM))
        else:
            f["veredicto"] = "LEJOS"
            f["motivo"] = "dMasa %.0f mm" % f["distMasa"]
    return referencia


# ------------------------------------------------------------
# Comandos principales
# ------------------------------------------------------------

def diagnostico(huMin=None, huMax=None):
    """Analiza, aplica el veredicto del v7.1 y crea TODAS las islas
    relevantes como segmentos. No elimina ninguna."""
    volumeNode = _volumen()
    if volumeNode is None:
        return None
    hmin = HU_MIN if huMin is None else huMin
    hmax = HU_MAX if huMax is None else huMax
    _d["huMin"] = hmin

    print("=" * 66)
    print("CranioPlan - Bloque A %s" % CRANIOPLAN_DIAG)
    print("DIAG: MODO NO DESTRUCTIVO - no se elimina ninguna isla.")
    _encabezado(volumeNode)

    _borrarNodosPorNombre(NOMBRE_NODO_DIAG)
    _d["colores"] = {}
    _d["creadas"] = []

    volArr = slicer.util.arrayFromVolume(volumeNode)
    esp = volumeNode.GetSpacing()
    espaciadoZYX = (esp[2], esp[1], esp[0])

    labels, filas, info = _analizar(volArr, espaciadoZYX, hmin, hmax)
    if labels is None or not filas:
        print("DIAG: no hay islas por encima de %.2f cm3 en %d-%d HU."
              % (PISO_LISTADO_CM3, hmin, hmax))
        return None

    _d["nIslas"] = info["nIslas"]
    _d["totalCM3"] = info["totalCM3"]
    print("DIAG: islas detectadas : %d  (%d por encima de %.2f cm3)"
          % (info["nIslas"], len(filas), PISO_LISTADO_CM3))
    print("DIAG: TOTAL en mascara : %.2f cm3" % info["totalCM3"])
    if info["valorRelleno"] is not None:
        print("DIAG: relleno fuera del FOV en %d HU (borde circular activo)."
              % info["valorRelleno"])
    else:
        print("DIAG: sin relleno detectable: el borde circular NO se testea.")

    referencia = _veredictoV71(labels, filas, espaciadoZYX)
    if referencia is None:
        print("DIAG: ninguna isla pasa el filtro de material. Revisa HU o serie.")
        _imprimirTabla(filas)
        return None
    print("DIAG: referencia -> %.2f cm3, cortical %.1f%%, env %.1f%%"
          % (referencia["vol"], referencia["fracCort"] * 100.0,
             (referencia["env"] or 0) * 100.0))

    _imprimirTabla(filas)

    # --- crear los segmentos ---
    segNode = slicer.mrmlScene.AddNewNodeByClass(
        'vtkMRMLSegmentationNode', NOMBRE_NODO_DIAG)
    segNode.CreateDefaultDisplayNodes()
    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segmentacion = segNode.GetSegmentation()

    aCrear, resto = filas[:MAX_SEGMENTOS], filas[MAX_SEGMENTOS:]
    buf = np.zeros(labels.shape, dtype=np.uint8)
    creadas = []
    for i, f in enumerate(aCrear, start=1):
        segId = segmentacion.AddEmptySegment(
            "", "%02d_%s_%.2fcm3" % (i, f["veredicto"], f["vol"]))
        caja = f["caja"]
        buf[caja] = (labels[caja] == f["etiqueta"]).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(buf, segNode, segId, volumeNode)
        buf[caja] = 0
        color = COLOR.get(f["veredicto"], (0.6, 0.6, 0.6))
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

    _d.update({"labels": labels, "espaciadoZYX": espaciadoZYX,
               "filas": filas, "creadas": creadas, "segNode": segNode})

    print("")
    print("DIAG: nodo '%s' con %d segmento(s)."
          % (NOMBRE_NODO_DIAG, segmentacion.GetNumberOfSegments()))
    print("DIAG: colores -> hueso=ACEPTADA  naranja=DUDOSA  turquesa=MATERIAL")
    print("DIAG:            azul=LOSA  violeta=LEJOS")
    ayuda()
    return filas


def _imprimirTabla(filas):
    print("")
    print("--- TODAS LAS ISLAS (>= %.2f cm3) - NINGUNA FUE ELIMINADA ---"
          % PISO_LISTADO_CM3)
    print("  env% = fraccion de direcciones, desde el centro del craneo, en las")
    print("         que aparece la isla. SOLO SE MIDE, todavia no decide nada.")
    print("")
    print("  #   vol(cm3)   p50   p95   p99    max  cort%   env%   dRef   dMasa"
          "  espMin  %FOV  borde  veredicto  motivo")
    for i, f in enumerate(filas, start=1):
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "    -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "     -" if f["distMasa"] is None else "%6.1f" % f["distMasa"]
        marca = " *ref*" if f["referencia"] else ""
        print("%3d %9.2f %5.0f %5.0f %5.0f %6.0f %6.1f %s %s %s %7.1f %4.0f%%"
              "  %-5s  %-9s %s%s"
              % (i, f["vol"], f["p50"], f["p95"], f["p99"], f["max"],
                 f["fracCort"] * 100.0, env, dr, dm, min(f["extMM"]),
                 max(f["fracFOV"]) * 100, "si" if f["tocaBorde"] else "no",
                 f["veredicto"] or "?", f["motivo"], marca))

    porZona = {}
    for f in filas:
        porZona.setdefault(f["veredicto"], []).append(f["vol"])
    print("")
    print("  TOTAL en mascara : %8.2f cm3" % (_d.get("totalCM3") or 0.0))
    for z in ["ACEPTADA", "DUDOSA", "MATERIAL", "LOSA", "LEJOS"]:
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
    """Pinta esa isla de rojo; el resto vuelve a su color de veredicto."""
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    segmentacion = segNode.GetSegmentation()
    for f in _d["creadas"]:
        seg = segmentacion.GetSegment(f["segId"])
        if seg is None:
            continue
        if f["numero"] == numero:
            seg.SetColor(1.0, 0.05, 0.05)
        else:
            seg.SetColor(*_d["colores"].get(f["segId"], (0.6, 0.6, 0.6)))


def solo(numero):
    """Muestra unicamente esa isla."""
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
    """Deja a la vista SOLO lo que el v7.1 no acepta: la vista que contesta
    'se esta yendo hueso del craneo?'."""
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
    """Filtra por veredicto: ACEPTADA, DUDOSA, MATERIAL, LOSA, LEJOS."""
    segNode = _d.get("segNode")
    if segNode is None:
        print("DIAG: primero corre diagnostico().")
        return
    disp = segNode.GetDisplayNode()
    clave = nombre.upper()
    n = 0
    for f in _d["creadas"]:
        v = f["veredicto"] == clave
        disp.SetSegmentVisibility(f["segId"], v)
        n += 1 if v else 0
    print("DIAG: %d isla(s) con veredicto %s." % (n, clave))


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
          "env %.1f%%, dRef %.1f mm -> %s"
          % (numero, ras[0], ras[1], ras[2], f["vol"], f["fracCort"] * 100.0,
             (f["env"] or 0) * 100.0,
             f["distRef"] if f["distRef"] is not None else -1, f["veredicto"]))


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
        print("DIAG: '%s' no esta en la lista sugerida. Se guarda igual, pero "
              "si es un tipeo conviene corregirlo (etiquetas())." % clave)
    f["etiqueta_real"] = clave
    f["nota_real"] = nota
    esperado = "ACEPTADA" if clave in _CRANEO else "no-ACEPTADA"
    real = "ACEPTADA" if f["veredicto"] == "ACEPTADA" else "no-ACEPTADA"
    flecha = "" if (clave == "no_se" or esperado == real) else "   <-- DISCREPANCIA"
    print("DIAG: isla %d = %s  (v7.1 dijo %s | %.2f cm3, cort %.1f%%, "
          "env %.1f%%, dRef %.1f mm)%s"
          % (numero, clave, f["veredicto"], f["vol"], f["fracCort"] * 100.0,
             (f["env"] or 0) * 100.0,
             f["distRef"] if f["distRef"] is not None else -1, flecha))


def pendientes():
    faltan = [f for f in _d.get("creadas", []) if not f.get("etiqueta_real")]
    if not faltan:
        print("DIAG: todas las islas estan etiquetadas.")
        return []
    print("DIAG: faltan etiquetar %d isla(s):" % len(faltan))
    for f in faltan:
        print("   %2d  %8.2f cm3  cort %5.1f%%  env %5.1f%%  dRef %6.1f  -> %s"
              % (f["numero"], f["vol"], f["fracCort"] * 100.0,
                 (f["env"] or 0) * 100.0,
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
    print("--- DISCREPANCIAS v7.1 vs realidad ---")
    print("  FALSOS NEGATIVOS (hueso craneal NO aceptado): %d" % len(fn))
    for f in fn:
        print("     %2d  %-14s %7.2f cm3  cort %5.1f%%  env %5.1f%%  "
              "dRef %6.1f mm  -> %s"
              % (f["numero"], f["etiqueta_real"], f["vol"],
                 f["fracCort"] * 100.0, (f["env"] or 0) * 100.0,
                 f["distRef"] if f["distRef"] is not None else -1, f["veredicto"]))
    print("  FALSOS POSITIVOS (basura aceptada): %d" % len(fp))
    for f in fp:
        print("     %2d  %-14s %7.2f cm3  cort %5.1f%%  env %5.1f%%  "
              "dRef %6.1f mm"
              % (f["numero"], f["etiqueta_real"], f["vol"],
                 f["fracCort"] * 100.0, (f["env"] or 0) * 100.0,
                 f["distRef"] if f["distRef"] is not None else -1))
    print("  (vertebras y 'no_se' no cuentan en ninguna de las dos listas)")
    print("--------------------------------------")
    return fn, fp


def envolvencia_por_etiqueta():
    """Resume env% agrupado por lo que la isla ES en realidad. Esta es la
    tabla que va a decidir si la envolvencia se convierte en criterio de
    eliminacion automatica o se descarta."""
    grupos = {}
    for f in _d.get("creadas", []):
        e = f.get("etiqueta_real")
        if not e or f["env"] is None:
            continue
        grupos.setdefault(e, []).append(f["env"] * 100.0)
    if not grupos:
        print("DIAG: todavia no hay islas etiquetadas con env% medido.")
        return
    print("")
    print("--- env%% POR TIPO DE OBJETO (este caso) ---")
    print("  etiqueta          n    min    med    max")
    for e in sorted(grupos, key=lambda k: -max(grupos[k])):
        v = sorted(grupos[e])
        med = v[len(v) // 2]
        print("  %-14s %4d  %5.1f  %5.1f  %5.1f" % (e, len(v), v[0], med, v[-1]))
    print("-------------------------------------------")


def reporte():
    """Bloque compacto, listo para copiar y pasar."""
    vol = _d.get("volumeNode")
    creadas = _d.get("creadas", [])
    if vol is None or not creadas:
        print("DIAG: primero corre diagnostico().")
        return
    dims = vol.GetImageData().GetDimensions()
    esp = vol.GetSpacing()
    print("")
    print("========== REPORTE CRANIOPLAN - BLOQUE A v7.1 DIAG ==========")
    print("caso         : <<< COMPLETAR: nombre/numero de TC >>>")
    print("serie        : %s" % vol.GetName())
    print("dims         : %d x %d x %d" % dims)
    print("spacing      : %.3f x %.3f x %.3f mm" % esp)
    print("FOV          : %.1f x %.1f x %.1f mm"
          % (dims[0] * esp[0], dims[1] * esp[1], dims[2] * esp[2]))
    print("HU_MIN       : %s" % _d.get("huMin"))
    print("islas        : %s" % _d.get("nIslas"))
    print("total mascara: %.2f cm3" % (_d.get("totalCM3") or 0.0))
    print("")
    print("  #  vol_cm3   p50   p95   p99    max  cort%   env%   dRef   dMasa"
          " espMin  %FOV borde  veredicto  ES_EN_REALIDAD")
    for f in sorted(creadas, key=lambda x: x["numero"]):
        e = f.get("etiqueta_real") or "SIN_ETIQUETAR"
        nota = f.get("nota_real", "")
        env = "  n/d" if f["env"] is None else "%5.1f" % (f["env"] * 100.0)
        dr = "    -" if f["distRef"] is None else "%5.1f" % f["distRef"]
        dm = "     -" if f["distMasa"] is None else "%6.1f" % f["distMasa"]
        print("%3d %8.2f %5.0f %5.0f %5.0f %6.0f %6.1f %s %s %s %6.1f %4.0f%%"
              "  %-4s  %-9s  %s%s"
              % (f["numero"], f["vol"], f["p50"], f["p95"], f["p99"], f["max"],
                 f["fracCort"] * 100.0, env, dr, dm, min(f["extMM"]),
                 max(f["fracFOV"]) * 100, "si" if f["tocaBorde"] else "no",
                 f["veredicto"], e, ("  (%s)" % nota) if nota else ""))
    print("")
    discrepancias()
    envolvencia_por_etiqueta()
    faltan = [f for f in creadas if not f.get("etiqueta_real")]
    if faltan:
        print("ATENCION: quedan %d isla(s) sin etiquetar (pendientes())."
              % len(faltan))
    print("=============================================================")


def csv(ruta=None):
    """Exporta y ACUMULA. Un solo archivo con todas las TC es lo que
    necesitamos para calibrar; no lo pises entre pacientes."""
    creadas = _d.get("creadas", [])
    if not creadas:
        print("DIAG: primero corre diagnostico().")
        return
    if ruta is None:
        ruta = RUTA_CSV or os.path.join(
            os.path.expanduser("~"), "Documents", "CranioPlan_diag_islas_v71.csv")
    vol = _d.get("volumeNode")
    nuevo = not os.path.exists(ruta)
    with open(ruta, "a") as fh:
        if nuevo:
            fh.write("serie;hu_min;n;vol_cm3;p50;p95;p99;max;frac_cortical;"
                     "env;esp_min_mm;ext_max_mm;pct_fov;toca_borde;es_losa;"
                     "dist_ref_mm;dist_masa_mm;veredicto_v71;etiqueta_real;nota\n")
        for f in sorted(creadas, key=lambda x: x["numero"]):
            fh.write("%s;%s;%d;%.3f;%.0f;%.0f;%.0f;%.0f;%.4f;%s;%.2f;%.2f;%.1f;"
                     "%s;%s;%s;%s;%s;%s;%s\n"
                     % (vol.GetName() if vol else "?", _d.get("huMin"),
                        f["numero"], f["vol"], f["p50"], f["p95"], f["p99"],
                        f["max"], f["fracCort"],
                        "" if f["env"] is None else "%.4f" % f["env"],
                        min(f["extMM"]), max(f["extMM"]),
                        max(f["fracFOV"]) * 100,
                        "si" if f["tocaBorde"] else "no",
                        "si" if f["esLosa"] else "no",
                        "" if f["distRef"] is None else "%.2f" % f["distRef"],
                        "" if f["distMasa"] is None else "%.2f" % f["distMasa"],
                        f["veredicto"], f.get("etiqueta_real", ""),
                        f.get("nota_real", "")))
    print("DIAG: %s -> %s" % ("creado" if nuevo else "agregado a", ruta))
    print("DIAG: el archivo acumula casos. Esquema nuevo respecto del CSV del")
    print("      v6-DIAG (agrega env, ext_max, es_losa, dRef y dMasa), por eso")
    print("      el nombre cambia y no se mezclan.")


def limpiar():
    _borrarNodosPorNombre(NOMBRE_NODO_DIAG)
    _borrarNodosPorNombre(NOMBRE_NODO_CRUDO)
    _d["segNode"] = None
    _d["crudoNode"] = None
    _d["creadas"] = []
    print("DIAG: escena limpia.")


def ayuda():
    print("")
    print("Analisis:")
    print("  listar_volumenes() / usar_volumen(i)")
    print("  diagnostico()            -> todas las islas, ninguna borrada")
    print("  barrido_hu()             -> cuanto hueso aparece al bajar HU_MIN")
    print("  crudo()                  -> mascara HU sin ningun filtro")
    print("Inspeccion:")
    print("  ir_a(n) / solo(n) / resaltar(n) / mostrar_todas()")
    print("  ocultar_aceptadas()      -> ver solo lo que el v7.1 no acepta")
    print("  ver_veredicto('DUDOSA')  -> filtrar por veredicto")
    print("Etiquetado:")
    print("  etiquetas()              -> lista de etiquetas sugeridas")
    print("  et(n, 'mandibula')       -> registrar que es")
    print("  pendientes()             -> que falta etiquetar")
    print("  discrepancias()          -> donde el v7.1 se equivoca")
    print("  envolvencia_por_etiqueta() -> env% agrupado por tipo de objeto")
    print("  reporte()                -> bloque compacto para pasar")
    print("  csv()                    -> exporta y ACUMULA todos los casos")
    print("  limpiar()                -> borrar los nodos DIAG_*")


diagnostico()
