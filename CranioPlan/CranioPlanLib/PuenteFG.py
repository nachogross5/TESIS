# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - PUENTE F -> G  (version libreria, para el modulo)
# ============================================================
# Es el script desarrollo/valentino/PuenteF_G.py convertido en libreria. La
# logica NO se toco. Los unicos cambios son de envoltorio:
#
#   - print() -> _p(), que escribe en la consola de Python Y en la interfaz;
#   - puente() ademas de imprimir DEVUELVE un diccionario con el resultado
#     (renombrados, piezas por tipo, preexistentes, desconocidas);
#   - los nombres de nodo y de carpeta vienen de Comun.py;
#   - no imprime nada al importarse.
#
# La logica del modulo lo llama automaticamente despues de cortar(), sin
# boton: el medico no tiene que saber que existe.
#
# ------------------------------------------------------------
# QUE PROBLEMA RESUELVE
# ------------------------------------------------------------
# Hasta corte_V17.py, cortar() hacia AddEmptySegment(nombre). Con un solo
# argumento Slicer toma el nombre como ID del segmento y le autogenera un
# NOMBRE 'Segment_N'; ExportSegmentsToModels bautiza los modelos con el
# nombre, asi que las piezas llegaban al Bloque G como 'Segment_5',
# 'Segment_17', etc. El nombre correcto estaba guardado en el ID del segmento:
# este puente lo recupera.
#
# CranioPlanLib/BloqueF.py ya tiene el arreglo AddEmptySegment(nombre,
# nombre). Con eso el puente detecta que "los nombres ya estan bien", no
# renombra nada y solo escribe los atributos. Se deja igual porque:
#   - protege si el modulo abre una escena cortada con una version vieja;
#   - escribe atributos legibles por maquina en cada modelo, que sobreviven a
#     un renombrado a mano en el panel Data y al guardado de la escena:
#
#       CranioPlan.origen   = "corte" | "preexistente"
#       CranioPlan.tipo     = "tapa" | "resto" | "fragmento" | "hueso"
#       CranioPlan.curva    = nombre de la curva (solo en tapas)
#
#     El Bloque G v1.1 decide la movilidad por geometria y todavia NO lee
#     estos atributos.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS
# ------------------------------------------------------------
# - Si un corte grabo pero NO separo, la isla madre no se dividio y esa pieza
#   sale clasificada como 'Hueso_N' (preexistente) aunque el corte la haya
#   tocado. El Bloque F lo detecta (diag["separo"]) pero no lo propaga al
#   nombre. El Paso 3 avisa cuando alguna linea dio "SEPARACION: NO".
# - Necesita el nodo de trabajo 'CranioPlan_Corte' en la escena. Si se borro
#   (limpiar() del Bloque F), ya no se puede recuperar el mapa.
# - Si se corta dos veces en la misma escena, los modelos del primer corte
#   siguen en la carpeta y tambien reciben atributos (identificado por
#   revision de codigo, no observado todavia en casos reales).
# ============================================================

import vtk
import slicer

from .Comun import NOMBRE_SEG_TRABAJO, NOMBRE_CARPETA_SH

PUENTE_FG_VERSION = "v1"

PREFIJOS_TIPO = (
    ("Tapa_",      "tapa",       "corte"),
    ("Resto_",     "resto",      "corte"),
    ("Fragmento_", "fragmento",  "corte"),
    ("Hueso_",     "hueso",      "preexistente"),
)


# ------------------------------------------------------------
# Salida de texto
# ------------------------------------------------------------
_SALIDA = {"log": print}


def _p(mensaje):
    _SALIDA["log"](mensaje)


def configurar(log=None):
    """Enchufa el log del modulo (None = print)."""
    _SALIDA["log"] = log if log is not None else print


# ============================================================
# LOGICA (identica a PuenteF_G.py)
# ============================================================
def _sh():
    return slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)


def _carpeta():
    shNode = _sh()
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(shNode.GetSceneItemID(), hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if shNode.GetItemName(itemId) == NOMBRE_CARPETA_SH:
            return shNode, itemId
    return shNode, None


def _modelos(shNode, itemId):
    salida = []

    def recorrer(parent):
        hijos = vtk.vtkIdList()
        shNode.GetItemChildren(parent, hijos)
        for i in range(hijos.GetNumberOfIds()):
            hid = hijos.GetId(i)
            nodo = shNode.GetItemDataNode(hid)
            if nodo is not None and nodo.IsA("vtkMRMLModelNode"):
                salida.append(nodo)
            elif nodo is None:
                recorrer(hid)

    recorrer(itemId)
    return salida


def _clasificar(nombreReal):
    """Devuelve (tipo, origen, curva) a partir del nombre real de la pieza."""
    for prefijo, tipo, origen in PREFIJOS_TIPO:
        if nombreReal.startswith(prefijo):
            curva = ""
            if tipo == "tapa":
                curva = nombreReal[len(prefijo):]
                # las tapas multiples llevan sufijo _a, _b...
                if len(curva) > 2 and curva[-2] == "_" and curva[-1].isalpha():
                    curva = curva[:-2]
            return tipo, origen, curva
    return "desconocido", "desconocido", ""


def puente():
    """
    Recupera los nombres reales y escribe los atributos CranioPlan.

    Devuelve un dict para la interfaz, o None si no pudo correr:
        {"renombrados": int, "porTipo": {tipo: [nombres]},
         "preexistentes": [nombres], "desconocidas": [nombres],
         "mapaUsado": bool}
    """
    try:
        segTrabajo = slicer.util.getNode(NOMBRE_SEG_TRABAJO)
    except Exception:
        _p(f"ERROR: no encuentro '{NOMBRE_SEG_TRABAJO}'.")
        _p("  Si corriste limpiar() del Bloque F, hay que volver a cortar().")
        return None

    seg = segTrabajo.GetSegmentation()
    n = seg.GetNumberOfSegments()
    if n == 0:
        _p("ERROR: el nodo de trabajo no tiene segmentos.")
        return None

    # nombre autogenerado -> nombre real (el ID guarda el nombre de Nacho)
    mapa = {}
    for i in range(n):
        sid = seg.GetNthSegmentID(i)
        s = seg.GetNthSegment(i)
        nombreSeg = s.GetName() if s is not None else ""
        if nombreSeg and nombreSeg != sid:
            mapa[nombreSeg] = sid

    shNode, itemId = _carpeta()
    if itemId is None:
        _p(f"ERROR: no encuentro la carpeta '{NOMBRE_CARPETA_SH}'.")
        return None
    modelos = _modelos(shNode, itemId)
    if not modelos:
        _p("ERROR: la carpeta no tiene modelos.")
        return None

    _p("--- PUENTE F -> G ---")
    if mapa:
        _p(f"Mapa recuperado desde los IDs de segmento: {len(mapa)} entrada(s).")
    else:
        _p("Los nombres ya estan bien (el Bloque F fue arreglado).")
        _p("Solo escribo los atributos.")

    _p("")
    _p(f"{'ANTES':<16}{'AHORA':<22}{'TIPO':<12}{'ORIGEN':<15}CURVA")
    _p("-" * 78)

    renombrados = 0
    usados = set()
    for m in modelos:
        antes = m.GetName()
        real = mapa.get(antes, antes)

        if real != antes:
            candidato = real
            k = 1
            while candidato in usados:
                candidato = f"{real}_{k}"
                k += 1
            m.SetName(candidato)
            real = candidato
            renombrados += 1
        usados.add(real)

        tipo, origen, curva = _clasificar(real)
        m.SetAttribute("CranioPlan.origen", origen)
        m.SetAttribute("CranioPlan.tipo", tipo)
        if curva:
            m.SetAttribute("CranioPlan.curva", curva)

        _p(f"{antes:<16}{real:<22}{tipo:<12}{origen:<15}{curva}")

    _p("-" * 78)
    _p(f"{renombrados} modelo(s) renombrado(s), {len(modelos)} con atributos.")

    desconocidos = [m.GetName() for m in modelos
                    if m.GetAttribute("CranioPlan.tipo") == "desconocido"]
    if desconocidos:
        _p("")
        _p(f"AVISO: {len(desconocidos)} pieza(s) sin clasificar: {desconocidos}")
        _p("  No matchean ningun prefijo conocido. En el Bloque G van a")
        _p("  decidirse por proximidad a las curvas.")

    preex = [m.GetName() for m in modelos
             if m.GetAttribute("CranioPlan.origen") == "preexistente"]
    _p("")
    _p(f"Piezas PREEXISTENTES (deberian ser mandibula, vertebras, tubuladuras): "
       f"{preex if preex else 'ninguna'}")
    _p("Si alguna pieza craneal aparece aca, es que su corte grabo pero NO")
    _p("separo: revisa el 'SEPARACION: NO' en la salida de cortar().")
    _p("")
    _p("Listo. Ahora corre preparar() del Bloque G.")

    # --- envoltorio del modulo: resumen para la interfaz ---
    porTipo = {}
    for m in modelos:
        porTipo.setdefault(m.GetAttribute("CranioPlan.tipo"), []).append(m.GetName())
    return {"renombrados": renombrados, "porTipo": porTipo,
            "preexistentes": preex, "desconocidas": desconocidos,
            "mapaUsado": bool(mapa)}
