# -*- coding: utf-8 -*-
# ============================================================
# CRANIOPLAN - LIBRERIA DE BLOQUES
#
# Cada bloque es un archivo aparte para que se pueda trabajar en paralelo sin
# pisarse: Valentino toca la interfaz (CranioPlan.py) y Nacho toca los
# algoritmos (BloqueF.py, BloqueA.py) sin que un cambio de uno obligue al otro
# a resolver un conflicto dentro del mismo archivo de 4000 lineas.
#
#   Comun.py   - nombres de nodos compartidos y logger
#   BloqueC.py - eleccion y carga de la serie DICOM correcta
#   BloqueA.py - preparacion del craneo (segmentacion + revision de piezas)
#   BloqueF.py - corte / osteotomias
#   BloqueG.py - reacomodamiento de las piezas post-corte
#
# Los bloques se comunican por NOMBRE DE NODO en la escena de Slicer, no por
# variables de Python. Esa es la razon por la que cada uno sigue funcionando
# tambien pegado suelto en la consola: el estado vive en la escena.
# ============================================================

from . import Comun          # noqa: F401
from . import BloqueC        # noqa: F401
from . import BloqueA        # noqa: F401
from . import BloqueF        # noqa: F401
from . import BloqueG        # noqa: F401


def recargar():
    """
    Vuelve a leer los cinco archivos desde el disco.

    Sirve durante el desarrollo: Slicer cachea los modulos Python importados,
    asi que despues de editar BloqueF.py el boton 'Reload' del modulo recarga
    CranioPlan.py pero NO la libreria, y uno queda probando el codigo viejo sin
    darse cuenta. El widget llama a esto en cada Reload.
    """
    import importlib
    for modulo in (Comun, BloqueC, BloqueA, BloqueF, BloqueG):
        importlib.reload(modulo)
