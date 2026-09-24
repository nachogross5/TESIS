# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - COMUN
#
# Lo que comparten TODOS los bloques: los nombres de los nodos de la escena
# y la funcion de log.
#
# POR QUE ESTE ARCHIVO EXISTE
# ---------------------------
# Los bloques se comunican entre si por NOMBRE DE NODO, no por variables:
# el Bloque A deja un nodo llamado "Craneo_Final" y el Bloque F lo busca por
# ese nombre. Mientras cada bloque era un script suelto, cada uno tenia su
# propia copia de esos nombres al principio del archivo, y bastaba con que
# alguien renombrara uno para que el siguiente bloque fallara con un
# "no encuentro el nodo..." que no explicaba nada.
#
# Con los nombres aca, hay una sola definicion. Si manana el Bloque F cambia
# a que espera, se cambia en un solo lugar y los tres bloques siguen
# hablandose.
# ============================================================

# --- Nodos que produce el Bloque A y consume el Bloque F ---
NOMBRE_NODO_SEGMENTACION = "Craneo_Automatico"   # vtkMRMLSegmentationNode
NOMBRE_SEGMENTO_CRANEO   = "Craneo_Final"        # segmento adentro de ese nodo
NOMBRE_MODELO_CRANEO     = "Craneo_Final"        # vtkMRMLModelNode (malla)

# --- Nodos que produce el Bloque F y consume el Bloque G ---
NOMBRE_SEG_TRABAJO     = "CranioPlan_Corte"
NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"
NOMBRE_MODELO_PICOS    = "CranioPlan_Picos"
NOMBRE_FID_SIN_PASAR   = "CranioPlan_NoAtraviesa"

# --- Nodos propios del Bloque G ---
NOMBRE_MOLDE           = "CranioPlan_Molde"
NOMBRE_LANDMARKS       = "CranioPlan_LineaMedia"

# --- Prefijo de las curvas de osteotomia que crea la interfaz ---
PREFIJO_CURVA_CORTE    = "Corte_"


class Registro:
    """
    Logger minimo con dos salidas a la vez: la consola de Python de Slicer
    (para nosotros, con todo el detalle) y un callback opcional (para la
    interfaz, que muestra un resumen al medico).

    Los bloques NUNCA llaman a print() directamente. Motivo: cuando el bloque
    corre desde el modulo, el medico no ve la consola de Python, asi que un
    aviso importante impreso ahi es un aviso que no existe. Con el callback,
    el mismo mensaje puede terminar en un cartel de la interfaz sin duplicar
    el codigo del bloque.
    """

    def __init__(self, callback=None, prefijo=""):
        self.callback = callback
        self.prefijo = prefijo
        self.avisos = []      # mensajes marcados como importantes
        self.lineas = []      # historial completo, para el reporte

    def __call__(self, mensaje, importante=False):
        texto = "%s%s" % (self.prefijo, mensaje)
        print(texto)
        self.lineas.append(texto)
        if importante:
            self.avisos.append(str(mensaje))
        if self.callback is not None:
            try:
                self.callback(str(mensaje), importante)
            except Exception:
                pass

    def aviso(self, mensaje):
        self(mensaje, importante=True)

    def limpiar(self):
        self.avisos = []
        self.lineas = []
