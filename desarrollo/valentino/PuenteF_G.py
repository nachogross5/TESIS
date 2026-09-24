# ============================================================
# CRANIOPLAN - PUENTE F -> G
#
# Repara el contrato entre el Bloque F y el Bloque G SIN tocar el archivo de
# Nacho. Se corre despues de cortar() y antes de preparar().
#
# ------------------------------------------------------------
# QUE PROBLEMA RESUELVE
# ------------------------------------------------------------
# En cortar() (v12 hasta v17) esta linea:
#
#     sid = segTrabajo.GetSegmentation().AddEmptySegment(nombre)
#
# La firma real de vtkSegmentation es:
#
#     AddEmptySegment(segmentId, segmentName, color)
#
# Con UN solo argumento, Slicer lo interpreta como el ID, no como el nombre.
# Resultado: el segmento queda con ID 'Fragmento_1' y NOMBRE vacio, y Slicer
# le autogenera 'Segment_5'. Despues ExportSegmentsToModels bautiza los
# modelos con el NOMBRE, no con el ID -> los fragmentos llegan al Bloque G
# llamandose 'Segment_5', 'Segment_17', etc.
#
# Consecuencia en el Bloque G: la regla "Hueso_ = isla preexistente = fija"
# no matchea nunca, y la mandibula, las vertebras y las tubuladuras quedan
# movibles.
#
# LA BUENA NOTICIA: la clasificacion de Nacho nunca estuvo rota. El nombre
# correcto esta guardado en el ID del segmento. Este script lo recupera.
#
# ARREGLO DEFINITIVO en el Bloque F (una sola linea, hacerlo igual):
#     sid = segTrabajo.GetSegmentation().AddEmptySegment(nombre, nombre)
#
# ------------------------------------------------------------
# QUE HACE ESTE SCRIPT
# ------------------------------------------------------------
# 1. Recorre el nodo de segmentacion de trabajo y arma el mapa
#    nombre_autogenerado (Segment_N)  ->  nombre_real (ID del segmento).
# 2. Renombra los modelos de la carpeta CranioPlan_Fragmentos.
# 3. Escribe atributos legibles por maquina en cada modelo, para que el
#    Bloque G no dependa de parsear strings:
#
#       CranioPlan.origen   = "corte" | "preexistente"
#       CranioPlan.tipo     = "tapa" | "resto" | "fragmento" | "hueso"
#       CranioPlan.curva    = nombre de la curva (solo en tapas)
#
#    Apoyar la semantica solo en un prefijo de nombre es fragil: cualquiera
#    renombra una pieza en el panel Data y se pierde. Los atributos MRML
#    sobreviven al guardado de la escena.
#
# 4. Imprime la tabla resultante para verificar a ojo.
#
# ------------------------------------------------------------
# LIMITES CONOCIDOS
# ------------------------------------------------------------
# - Si el Bloque F ya fue arreglado, los nombres ya estan bien: el script lo
#   detecta, no renombra nada, y solo escribe los atributos.
# - Si un corte grabo pero NO separo, la isla madre no se dividio y esa pieza
#   sale clasificada como 'Hueso_N' (preexistente) aunque el corte la haya
#   tocado. El Bloque F lo detecta (diag["separo"]) pero no lo propaga al
#   nombre. Revisa la salida de cortar(): si dice "SEPARACION: NO", esa pieza
#   va a llegar a G marcada como fija y hay que liberarla a mano.
# - El script asume que el nodo de trabajo 'CranioPlan_Corte' sigue en la
#   escena. Si corriste limpiar() del Bloque F, ya no se puede recuperar el
#   mapa y hay que volver a cortar.
# ============================================================

NOMBRE_SEG_TRABAJO   = "CranioPlan_Corte"
NOMBRE_CARPETA_SH    = "CranioPlan_Fragmentos"

PREFIJOS_TIPO = (
    ("Tapa_",      "tapa",       "corte"),
    ("Resto_",     "resto",      "corte"),
    ("Fragmento_", "fragmento",  "corte"),
    ("Hueso_",     "hueso",      "preexistente"),
)


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
    """Recupera los nombres reales y escribe los atributos CranioPlan."""
    try:
        segTrabajo = slicer.util.getNode(NOMBRE_SEG_TRABAJO)
    except Exception:
        print(f"ERROR: no encuentro '{NOMBRE_SEG_TRABAJO}'.")
        print("  Si corriste limpiar() del Bloque F, hay que volver a cortar().")
        return

    seg = segTrabajo.GetSegmentation()
    n = seg.GetNumberOfSegments()
    if n == 0:
        print("ERROR: el nodo de trabajo no tiene segmentos.")
        return

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
        print(f"ERROR: no encuentro la carpeta '{NOMBRE_CARPETA_SH}'.")
        return
    modelos = _modelos(shNode, itemId)
    if not modelos:
        print("ERROR: la carpeta no tiene modelos.")
        return

    print("--- PUENTE F -> G ---")
    if mapa:
        print(f"Mapa recuperado desde los IDs de segmento: {len(mapa)} entrada(s).")
    else:
        print("Los nombres ya estan bien (el Bloque F fue arreglado).")
        print("Solo escribo los atributos.")

    print("")
    print(f"{'ANTES':<16}{'AHORA':<22}{'TIPO':<12}{'ORIGEN':<15}CURVA")
    print("-" * 78)

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

        print(f"{antes:<16}{real:<22}{tipo:<12}{origen:<15}{curva}")

    print("-" * 78)
    print(f"{renombrados} modelo(s) renombrado(s), {len(modelos)} con atributos.")

    desconocidos = [m.GetName() for m in modelos
                    if m.GetAttribute("CranioPlan.tipo") == "desconocido"]
    if desconocidos:
        print("")
        print(f"AVISO: {len(desconocidos)} pieza(s) sin clasificar: {desconocidos}")
        print("  No matchean ningun prefijo conocido. En el Bloque G van a")
        print("  decidirse por proximidad a las curvas.")

    preex = [m.GetName() for m in modelos
             if m.GetAttribute("CranioPlan.origen") == "preexistente"]
    print("")
    print(f"Piezas PREEXISTENTES (deberian ser mandibula, vertebras, tubuladuras): "
          f"{preex if preex else 'ninguna'}")
    print("Si alguna pieza craneal aparece aca, es que su corte grabo pero NO")
    print("separo: revisa el 'SEPARACION: NO' en la salida de cortar().")
    print("")
    print("Listo. Ahora corre preparar() del Bloque G.")


print("Puente F -> G cargado.")
print("  puente()   -> recupera nombres reales y escribe atributos CranioPlan")
print("")
print("ARREGLO DEFINITIVO en el Bloque F, cambiar:")
print("    AddEmptySegment(nombre)        ->    AddEmptySegment(nombre, nombre)")
