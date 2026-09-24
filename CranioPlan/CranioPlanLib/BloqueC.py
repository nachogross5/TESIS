# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - BLOQUE C v2.1
#        ELECCION Y CARGA DE LA SERIE DE TRABAJO
# ============================================================
#
# ESTE ARCHIVO ES EL BLOQUE C v2 QUE YA ESTABA EN EL MODULO, sacado a
# libreria. El criterio de clasificacion NO se toco: esta validado contra 6
# estudios reales del Garrahan y acierta en los 6. Lo unico que se agrego es
# lo de v2.1, al final de este encabezado.
#
# ------------------------------------------------------------
# POR QUE SE REESCRIBIO EN SU MOMENTO (evidencia de los 6 estudios)
# ------------------------------------------------------------
# (1) La regla vieja descartaba toda serie cuya DESCRIPCION contuviera
#     "axial", "coronal" o "sagittal". En el caso 403 la serie buena se llama
#     "Hueso 0.5 Axial VOLUME VOL SCAN": la descartaba junto con todas las
#     demas del estudio y terminaba eligiendo el localizador de 2 imagenes.
#     Una adquisicion volumetrica ES axial: la palabra no distingue nada. Lo
#     que distingue a un reformateo es que es DERIVED.
#
# (2) Miraba solo el primer paciente y el primer estudio (patientUIDs[0],
#     studyUIDs[0]). Un zip del portal puede traer mas de un estudio del
#     mismo paciente, y las series del resto no se miraban nunca.
#
# (3) Ordenaba por SliceThickness, que es el espesor NOMINAL del corte. En
#     estos equipos varias series empatan en 0.5 siendo muy distintas: el
#     intervalo REAL de reconstruccion va de 0.25 a 0.5 mm. Ahora se calcula
#     con ImagePositionPatient del primer y el ultimo corte.
#
# ------------------------------------------------------------
# CRITERIOS, todos objetivos (tags DICOM, no nombres)
# ------------------------------------------------------------
#   - Se recorren TODOS los estudios y series de la carpeta, recursivamente.
#   - Se descarta por SOPClassUID lo que no es imagen de CT (los VR 3D y los
#     informes SUMMARY son Secondary Capture).
#   - Se descarta por ImageType: LOCALIZER (scout), DERIVED (reformateos MPR)
#     y SECONDARY (KEY_IMAGES). En los 6 estudios, la serie correcta es
#     SIEMPRE ORIGINAL\PRIMARY y ninguna otra lo es.
#   - Se descarta el plano coronal/sagital por ImageOrientationPatient
#     (geometria), no por la palabra en la descripcion.
#   - Se exige que sea volumetrica: >= 60 imagenes, cortes <= 2 mm y
#     cobertura >= 40 mm en Z.
#   - Entre las que quedan, el orden es: HUESO primero, despues CEREBRO,
#     despues cualquier otra volumetrica; y dentro de cada grupo, menor
#     espaciado real y despues mas imagenes.
#
# ------------------------------------------------------------
# QUE AGREGA v2.1
# ------------------------------------------------------------
# [1] EL CEREBRO ES UN ESCALON PROPIO. Antes habia dos grupos: hueso y "todo
#     lo demas". Ahora son tres. Sobre los casos vistos da el mismo resultado
#     (cuando no hay serie de hueso, la de cerebro es igual la de menor
#     espaciado), pero deja el criterio ESCRITO en vez de que dependa de que
#     el cerebro casualmente sea la mas fina. Es lo que se pidio: si no hay
#     hueso, la de cerebro con mas cortes.
#
# [2] SE DETECTA EL EMPATE. Si las dos primeras candidatas tienen la misma
#     descripcion y practicamente la misma cantidad de imagenes, son dos
#     reconstrucciones del mismo protocolo y el algoritmo NO puede saber cual
#     quiere el medico. Se marcan como ambiguas para que la interfaz lo diga
#     en vez de elegir uno callado.
#
# [3] EL ZIP SE MANEJA ACA. Descomprimir el .zip del portal era parte del
#     widget; paso a la libreria, que es donde puede reusarlo cualquiera.
#
# ------------------------------------------------------------
# LIMITE CONOCIDO
# ------------------------------------------------------------
# La lista de kernels oseos cubre Canon/Toshiba (FC30...), Siemens
# (H60/H70/B60/B70), GE (BONE...) y Philips (YA/YB...). Un equipo con otra
# convencion caeria en el fallback volumetrico, que sigue siendo una serie
# valida. Por eso el Paso 1 muestra el ranking completo y permite cambiar la
# eleccion a mano.
# ============================================================

import os
import shutil
import tempfile
import zipfile
import re

import slicer
from DICOMLib import DICOMUtils

try:
    import pydicom
except ImportError:      # Slicer siempre lo trae, pero no dependemos de eso
    pydicom = None


CT_IMAGE_STORAGE = "1.2.840.10008.5.1.4.1.1.2"
MIN_IMAGENES_VOLUMETRICA = 60
MAX_SPACING_Z_MM = 2.0
MIN_EXTENSION_Z_MM = 40.0

KERNELS_HUESO_EXACTOS = {
    "FC30", "FC31", "FC35", "FC81", "FC82", "FC83",   # Canon / Toshiba
    "BONE", "BONEPLUS", "DETAIL", "EDGE",             # GE
    "YA", "YB", "YC", "YD", "EB", "EC",               # Philips
}
KERNELS_HUESO_PREFIJOS = ("H60", "H70", "B60", "B70", "BR6", "BR7", "HR6", "HR7")
PALABRAS_HUESO = ("hueso", "bone", "knochen", "osso")
PALABRAS_CEREBRO = ("cerebro", "brain", "neuro", "hirn")

# Dos series "empatan" si se llaman igual y la diferencia de imagenes es menor
# que esto. Dos reconstrucciones del mismo protocolo difieren en unas pocas
# imagenes; series genuinamente distintas difieren en decenas o cientos.
EMPATE_DIFERENCIA_ABSOLUTA = 10       # imagenes
EMPATE_DIFERENCIA_RELATIVA = 0.03     # 3 % de la mayor


# ============================================================
# CLASIFICACION  (nucleo: no depende de Slicer, se puede probar aparte)
# ============================================================
def _txt(v):
    return "" if v is None else str(v)


def _normalizar(texto):
    """Descripcion en minusculas, sin puntos y con un solo espacio, para poder
    comparar 'Hueso 0.5 Vol Vol.' con 'HUESO 0.5 VOL VOL' sin sorpresas."""
    t = _txt(texto).lower().replace(".", " ")
    return re.sub(r"\s+", " ", t).strip()


def esKernelDeHueso(kernel):
    k = _txt(kernel).upper().replace("\\", " ").replace("/", " ")
    for token in k.split():
        t = token.strip(" .,-")
        if t in KERNELS_HUESO_EXACTOS:
            return True
        if t.startswith(KERNELS_HUESO_PREFIJOS):
            return True
        if "BONE" in t or "SHARP" in t:
            return True
    return False


def esDescripcionDeHueso(descripcion):
    d = _txt(descripcion).lower()
    return any(p in d for p in PALABRAS_HUESO)


def esDescripcionDeCerebro(descripcion):
    d = _txt(descripcion).lower()
    return any(p in d for p in PALABRAS_CEREBRO)


def planoDesdeOrientacion(iop):
    """
    Plano de adquisicion a partir de ImageOrientationPatient.

    IOP son los vectores de fila y columna del plano de imagen; su producto
    vectorial es la normal, y el eje dominante de la normal da el plano. Es
    geometria real: detecta un reformateo coronal aunque se llame "VOL", y NO
    descarta una volumetrica por decir "Axial".
    """
    if not iop or len(iop) < 6:
        return "?"
    try:
        f = [float(v) for v in iop]
    except (TypeError, ValueError):
        return "?"
    fila, col = f[0:3], f[3:6]
    normal = [fila[1] * col[2] - fila[2] * col[1],
              fila[2] * col[0] - fila[0] * col[2],
              fila[0] * col[1] - fila[1] * col[0]]
    ejes = ["SAGITAL", "CORONAL", "AXIAL"]
    absn = [abs(v) for v in normal]
    i = absn.index(max(absn))
    return ejes[i] if absn[i] > 0.95 else "OBLICUO"


def motivoDeExclusion(s):
    """Motivo por el que una serie no puede ser la de trabajo, o None."""
    tipo = _txt(s.get("image_type")).upper()
    partes = [p.strip() for p in tipo.replace("/", "\\").split("\\") if p.strip()]

    if _txt(s.get("modality")).upper() not in ("CT", ""):
        return "modalidad %s (no es CT)" % s.get("modality")
    sop = _txt(s.get("sop_class"))
    if sop and sop != CT_IMAGE_STORAGE:
        return "no es imagen de CT (captura secundaria: VR, SUMMARY, informe)"
    if "LOCALIZER" in partes:
        return "localizador / scout"
    if partes and partes[0] == "DERIVED":
        return "reconstruccion derivada (MPR), no la adquisicion original"
    if len(partes) > 1 and partes[1] == "SECONDARY":
        return "captura secundaria (KEY_IMAGES, anotaciones)"

    plano = _txt(s.get("plano")).upper()
    if plano in ("CORONAL", "SAGITAL"):
        return "adquisicion en plano %s (no axial)" % plano

    n = s.get("n_imagenes") or 0
    if n < MIN_IMAGENES_VOLUMETRICA:
        return "solo %d imagen(es), no es volumetrica" % n
    sp = s.get("spacing_z_real")
    if sp is not None and sp > MAX_SPACING_Z_MM:
        return "cortes cada %.2f mm (demasiado gruesa)" % sp
    ext = s.get("extension_z_mm")
    if ext is not None and ext < MIN_EXTENSION_Z_MM:
        return "cobertura de solo %.0f mm en Z" % ext
    return None


def clasificarSeries(series):
    """
    Devuelve (ordenadas, descartadas), de mejor a peor para hueso.

    El orden es por tres claves, en este orden:
      1. TIPO DE VENTANA: hueso (0) > cerebro (1) > cualquier otra (2).
         El hueso primero porque es la ventana donde se ve la cortical, que
         es lo que necesita la segmentacion. El cerebro segundo porque cuando
         no hay serie de hueso es la unica volumetrica fina que mandan.
      2. MENOR ESPACIADO REAL. Es lo mismo que decir "la que mas cortes tiene"
         para una misma cobertura, pero medido en milimetros y calculado con
         las posiciones reales de los cortes, no con el tag nominal.
      3. MAS IMAGENES, para desempatar.
    """
    candidatas, descartadas = [], []
    for s in series:
        motivo = motivoDeExclusion(s)
        d = dict(s)
        if motivo:
            d["motivo"] = motivo
            descartadas.append(d)
            continue

        porKernel = esKernelDeHueso(s.get("kernel"))
        porTexto = esDescripcionDeHueso(s.get("series_desc"))
        esCerebro = esDescripcionDeCerebro(s.get("series_desc"))
        d["esHueso"] = porKernel or porTexto
        d["esCerebro"] = bool(esCerebro and not d["esHueso"])

        if porKernel and porTexto:
            d["razon"] = "hueso (descripcion + kernel %s)" % s.get("kernel")
        elif porKernel:
            d["razon"] = "hueso (kernel %s)" % s.get("kernel")
        elif porTexto:
            d["razon"] = "hueso (descripcion)"
        elif d["esCerebro"]:
            d["razon"] = "ventana de cerebro (sirve si no hay serie de hueso)"
        else:
            d["razon"] = "volumetrica de tejido blando"

        d["grupo"] = 0 if d["esHueso"] else (1 if d["esCerebro"] else 2)
        candidatas.append(d)

    candidatas.sort(key=lambda d: (
        d["grupo"],
        d.get("spacing_z_real") if d.get("spacing_z_real") is not None else 9e9,
        -(d.get("n_imagenes") or 0),
    ))
    return candidatas, descartadas


def detectarEmpate(ordenadas):
    """
    [2] Series indistinguibles entre si: mismo grupo, misma descripcion y
    practicamente la misma cantidad de imagenes.

    Cuando esto pasa, el estudio trae dos reconstrucciones del mismo
    protocolo y solo el medico sabe cual quiere. El algoritmo se declara
    incompetente a proposito en vez de tirar una moneda.

    Devuelve la lista de series empatadas (>= 2) o [].
    """
    if len(ordenadas) < 2:
        return []
    primera = ordenadas[0]
    empatadas = [primera]
    for s in ordenadas[1:]:
        if s.get("grupo") != primera.get("grupo"):
            break
        if _normalizar(s.get("series_desc")) != _normalizar(primera.get("series_desc")):
            break
        n1 = primera.get("n_imagenes") or 0
        n2 = s.get("n_imagenes") or 0
        mayor = max(n1, n2, 1)
        if abs(n1 - n2) > max(EMPATE_DIFERENCIA_ABSOLUTA,
                              EMPATE_DIFERENCIA_RELATIVA * mayor):
            break
        empatadas.append(s)
    return empatadas if len(empatadas) > 1 else []


# ============================================================
# LECTURA DE LA CARPETA
# ============================================================
def _leerCabecera(ruta):
    if pydicom is None:
        return None
    try:
        return pydicom.dcmread(ruta, stop_before_pixels=True, force=True)
    except Exception:
        return None


def prepararCarpeta(ruta):
    """
    [3] Devuelve (carpeta, temporal).

    Si `ruta` es un .zip, lo descomprime en una carpeta temporal y devuelve
    esa carpeta junto con su ruta, para que quien la creo la pueda borrar
    despues. Si es una carpeta, la devuelve tal cual y temporal = None.
    """
    if not str(ruta).lower().endswith(".zip"):
        return ruta, None
    temporal = tempfile.mkdtemp(prefix="cranioplan_zip_")
    try:
        with zipfile.ZipFile(ruta) as z:
            z.extractall(temporal)
    except Exception:
        shutil.rmtree(temporal, ignore_errors=True)
        raise
    return temporal, temporal


def borrarTemporal(temporal):
    if temporal:
        shutil.rmtree(temporal, ignore_errors=True)


def analizarEstudio(carpeta, log=print):
    """
    Recorre la carpeta RECURSIVAMENTE con pydicom y devuelve el ranking de
    series.

    Se lee con pydicom en vez de la base DICOM de Slicer porque hacen falta
    tags (ImageType, ImageOrientationPatient, ImagePositionPatient del primer
    y ultimo corte) y asi se evita ademas depender del estado de la base. El
    recorrido recursivo cubre las carpetas anidadas del portal (exam/XXXX/).

    Devuelve {"ordenadas", "descartadas", "empatadas", "nEstudios",
              "nPacientes", "nArchivos"} o None.
    """
    log("Bloque C v2.1: analizando %s" % carpeta)
    if pydicom is None:
        log("Bloque C: pydicom no esta disponible.")
        return None

    porSerie = {}
    nArchivos = 0
    for raiz, _dirs, archivos in os.walk(carpeta):
        for nombre in archivos:
            if nombre.lower() in ("dicomdir", "version"):
                continue
            ds = _leerCabecera(os.path.join(raiz, nombre))
            nArchivos += 1
            if ds is None or not hasattr(ds, "SOPClassUID"):
                continue
            uid = str(getattr(ds, "SeriesInstanceUID", "sin_uid"))
            porSerie.setdefault(uid, []).append(ds)

    log("Bloque C: %d archivo(s) recorrido(s), %d serie(s)."
        % (nArchivos, len(porSerie)))
    if not porSerie:
        return None

    series = []
    for uid, lista in porSerie.items():
        primero = lista[0]
        iop = getattr(primero, "ImageOrientationPatient", None)
        plano = planoDesdeOrientacion(iop)
        eje = {"AXIAL": 2, "CORONAL": 1, "SAGITAL": 0}.get(plano, 2)

        posiciones = []
        for ds in lista:
            ipp = getattr(ds, "ImagePositionPatient", None)
            if ipp is not None and len(ipp) == 3:
                try:
                    posiciones.append(float(ipp[eje]))
                except (TypeError, ValueError):
                    pass
        spacing = extension = None
        if len(posiciones) >= 2:
            posiciones.sort()
            extension = posiciones[-1] - posiciones[0]
            if extension > 0:
                spacing = extension / (len(posiciones) - 1)

        series.append({
            "series_uid": uid,
            "series_num": str(getattr(primero, "SeriesNumber", "")),
            "series_desc": str(getattr(primero, "SeriesDescription", "")),
            "modality": str(getattr(primero, "Modality", "")),
            "sop_class": str(getattr(primero, "SOPClassUID", "")),
            "image_type": "\\".join(str(x) for x in
                                    (getattr(primero, "ImageType", []) or [])),
            "kernel": str(getattr(primero, "ConvolutionKernel", "")),
            "plano": plano,
            "n_imagenes": len(lista),
            "spacing_z_real": spacing,
            "extension_z_mm": extension,
            "study_uid": str(getattr(primero, "StudyInstanceUID", "")),
            "study_date": str(getattr(primero, "StudyDate", "")),
            "paciente": str(getattr(primero, "PatientName", "")),
        })

    ordenadas, descartadas = clasificarSeries(series)
    empatadas = detectarEmpate(ordenadas)

    log("Bloque C: --------- SERIES DEL ESTUDIO ---------")
    for i, s in enumerate(ordenadas):
        sp = "%.3f" % s["spacing_z_real"] if s["spacing_z_real"] else "?"
        log("Bloque C:   %s n#%-5s %-40s %4d img  %s mm  %s"
            % ("ELEGIDA " if i == 0 else "        ", s["series_num"],
               (s["series_desc"] or "")[:40], s["n_imagenes"], sp, s["razon"]))
    for s in descartadas:
        log("Bloque C:   descartada n#%-5s %-40s <- %s"
            % (s["series_num"], (s["series_desc"] or "")[:40], s["motivo"]))
    log("Bloque C: --------------------------------------")

    if empatadas:
        log("Bloque C: hay %d series con el mismo nombre y casi la misma "
            "cantidad de cortes. No se puede saber cual quiere el medico: "
            "hay que preguntar." % len(empatadas))

    return {
        "ordenadas": ordenadas,
        "descartadas": descartadas,
        "empatadas": empatadas,
        "nEstudios": len(set(s["study_uid"] for s in series if s["study_uid"])),
        "nPacientes": len(set(s["paciente"] for s in series if s["paciente"])),
        "nArchivos": nArchivos,
    }


def cargarSerie(carpeta, seriesUID, log=print):
    """Importa la carpeta a una base temporal y carga esa serie en la escena.
    Devuelve el vtkMRMLScalarVolumeNode o None."""
    try:
        with DICOMUtils.TemporaryDICOMDatabase() as db:
            DICOMUtils.importDicom(carpeta, db)
            ids = DICOMUtils.loadSeriesByUID([seriesUID])
            for nodeID in ids:
                node = slicer.mrmlScene.GetNodeByID(nodeID)
                if node is not None and node.IsA("vtkMRMLScalarVolumeNode"):
                    log("Bloque C: serie cargada -> %s" % node.GetName())
                    return node
    except Exception as e:
        log("Bloque C: error al cargar la serie: %s" % e)
    return None


def etiquetaDeSerie(s, recomendada=False):
    """Texto de una serie para el desplegable del Paso 1."""
    sp = "%.2f mm" % s["spacing_z_real"] if s.get("spacing_z_real") else "? mm"
    return "%s n#%s - %s  (%d cortes, %s, %s)" % (
        "*" if recomendada else "  ",
        s.get("series_num", "?"),
        s.get("series_desc") or "(sin descripcion)",
        s.get("n_imagenes", 0), sp, s.get("razon", ""))
