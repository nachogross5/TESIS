# ============================================================
# CRANIOPLAN - BLOQUE F v17
#
# ############################################################
# VUELVE LA BUSQUEDA DEL CAMINO MAS CORTO, PERO COHERENTE
# ############################################################
#
# v14 y v15 buscaban, en un abanico de direcciones, cual atravesaba el hueso
# por el camino mas corto. La idea era y es CORRECTA: en la ceja, cortar hacia
# arriba por encima del techo orbitario es un camino muchisimo mas corto que
# arrastrarse a lo largo de la placa. v16 la elimino porque producia los
# "tajos random".
#
# Pero el problema no era la idea, era COMO estaba implementada:
# CADA PUNTO ELEGIA SU ANGULO POR SU CUENTA, sin mirar a los vecinos. Si en un
# punto el mejor camino era +40 grados y en el de al lado era -40 (dos minimos
# casi empatados, cosa comun en geometria complicada), quedaban dos puntos
# VECINOS apuntando en direcciones opuestas: eso es literalmente un tajo
# abriendose en abanico. Y el suavizado que habia encima no ayudaba, porque
# promediar +40 y -40 da 0, que no es ninguna de las dos y suele ser
# justamente la direccion rasante.
#
# ---- [C1] LA ELECCION AHORA ES DE LA CURVA ENTERA, NO DE CADA PUNTO ----
# Se elige la SECUENCIA COMPLETA de angulos que minimiza el recorrido total,
# con la restriccion dura de que de un punto al siguiente el angulo NO PUEDA
# GIRAR mas de MAX_GIRO_POR_PASO_DEG. Eso se resuelve exacto con programacion
# dinamica (Viterbi) sobre la grilla de angulos.
#
# Consecuencias, todas buenas:
#   - Es IMPOSIBLE que dos puntos vecinos apunten a lados opuestos: el tajo
#     random no puede ocurrir por construccion, no por suavizado a posteriori.
#   - La curva "se compromete" con una solucion global. Si el conjunto de la
#     ceja conviene resolverlo yendo hacia arriba, todos los puntos van hacia
#     arriba; ninguno se escapa hacia la orbita porque ahi seria un poquito
#     mas corto, porque cambiar de bando cuesta caro en la suma total.
#   - Donde ninguna direccion atraviesa, ese punto recibe un costo alto pero
#     FINITO, asi la secuencia puede pasar por ahi sin romperse.
#
# ---- [C2] UNIONES QUE DE VERDAD CIERRAN ----
# Bug encontrado en v16: la regla "si la punta esta cerca de otra curva, es una
# union, no la extiendo" (puesta en v13 para evitar espolones) dejaba el hueco
# ABIERTO cuando las puntas quedaban a 0.5-4mm una de otra, que es exactamente
# la distancia que sale cuando uno dibuja tratando de que dos lineas se junten.
# En un craneo CERRADO, dos lineas abiertas solo separan si forman un lazo; un
# hueco de 2mm en una sien deja el colgajo colgado y no se separa NADA.
# Ahora la punta se extiende EXACTAMENTE hasta tocar la otra curva, mas
# SOLAPE_UNION_MM de garantia. Ni menos (hueco) ni de mas (espolon).
#
# ---- [C3] LO QUE DE VERDAD DECIDE LA SEPARACION ----
# El reporte de v16 decia "atraviesa el 95.7%", que suena bien y es engañoso:
# la separacion no es un promedio, es una CADENA. 16 puntos seguidos que no
# atraviesan son ~6mm de hueso macizo que sostienen todo el colgajo. Ahora se
# reporta la TIRADA CONTINUA MAS LARGA sin atravesar, en mm, que es el numero
# que hay que llevar a cero.
#
# ############################################################
# SE CONSERVA DE v16
# ############################################################
# Verificacion rayo por rayo (atraviesa / no atraviesa), tope en seco para los
# rayos que no atraviesan (no pueden hacer daño lejos del sitio), salida por
# espacio abierto 3D + margen, kerf = GROSOR/2, supersampleo de continuidad,
# filtro angular del vecindario de normales, guarda anti-contaminacion,
# verificacion de separacion por curva, higiene de escena.
#
# v16 hace MENOS pero lo hace SIEMPRE IGUAL:
#
#   1. Corta SIEMPRE PERPENDICULAR a la superficie donde dibujaste la curva.
#      Predecible, sin sorpresas, sin tajos en diagonal.
#
#   2. VERIFICA, rayo por rayo, si ese corte logro ATRAVESAR el hueso, y te lo
#      muestra. Donde no pudo, NO improvisa: se detiene, lo marca, y te avisa.
#
# La zona retro-orbitaria probablemente siga sin separarse, porque ahi la
# perpendicular a la superficie externa corre a lo largo de la placa y eso es
# geometria, no un bug. La diferencia es que ahora vas a VER exactamente en que
# puntos no atraviesa, en vez de descubrirlo despues mirando el modelo. La
# solucion en esa zona es quirurgica, no algoritmica: una segunda curva
# dibujada sobre la cara por la que si se puede entrar. Es, de hecho, lo que
# hace el cirujano: la osteotomia del techo orbitario es un corte aparte.
#
# ############################################################
# [V1] LA VERIFICACION POR RAYO (lo nuevo)
# ############################################################
#
# Por cada punto de la curva, el barrido avanza y se pregunta: sali del hueso
# de verdad, o me quede adentro? Se acepta la salida solo si se cumplen DOS
# condiciones a la vez:
#
#   (a) TEST 3D DE BURBUJA. El hueco tiene que ser un ESPACIO ABIERTO: aire que
#       sobrevive a una erosion de UMBRAL_ESPACIO_ABIERTO_MM. Un poro o un
#       espacio medular del diploe no sobrevive; la cavidad craneal, la orbita
#       y el exterior si. Este test mide el hueco EN VOLUMEN, no a lo largo de
#       una linea, y por eso es dificil de enganar.
#
#   (b) MARGEN A LO LARGO DEL RAYO. Ademas, desde ese punto tiene que haber
#       MARGEN_SALIDA_MM de aire continuo. Es la confirmacion de que no caimos
#       en una burbuja.
#
# El margen de 2mm solo NO alcanza para el diploe supraorbitario, donde los
# espacios medulares llegan a 3-4mm; por eso van los dos juntos.
#
# ############################################################
# [V2] SIN SALIDA => NO SE IMPROVISA (esto mata los tajos random)
# ############################################################
#
# Si un rayo NO logra atravesar, v14/v15 seguian barriendo hasta el tope de
# 25mm, y si ademas la direccion venia inclinada, eso dejaba una lengua de
# corte de varios centimetros en diagonal: el tajo random.
#
# Ahora un rayo que no atraviesa se CORTA EN SECO a
# TOPE_SIN_SALIDA_FACTOR veces el espesor tipico de esa misma curva. No puede
# hacer dano lejos del sitio. Y el punto queda marcado.
#
# ############################################################
# QUE SE CONSERVA (arreglos que no tocan la direccion)
# ############################################################
# [B1] kerf = GROSOR/2. v13 y v14 usaban (GROSOR - spacing)/2, que con spacing
#      0.45 no incluia ningun voxel vecino: la ranura quedaba de UN voxel
#      (0.45mm), mas fina que la de v12.
# [B2] supersampleo de continuidad: se insertan rayos interpolados hasta que
#      dos vecinos no puedan separarse mas de medio voxel a ninguna
#      profundidad. Sin esto la cortina sale perforada y quedan puentes.
# [B4] salida por espacio abierto 3D (ahora parte de [V1]).
# v13: filtro angular del vecindario de normales (solo entran vertices de la
#      MISMA cara: ni la tabla interna, ni los poros, ni la placa vecina),
#      guarda anti-contaminacion en la propagacion de signos, ventana de
#      suavizado en mm, no extender puntas que se unen a otra curva,
#      muestreo derivado del spacing, higiene de escena.
# v15: verificacion de separacion curva por curva.
#
# FLUJO:
#   Bloque A+B -> confirmar_craneo() -> enviar_a_planner()
#   -> dibujar curvas -> cargar este archivo
#   -> diagnosticar()   <-- MIRA LOS PICOS: los rojos no atraviesan
#   -> cortar()
#
# REQUISITO: scipy.
# ============================================================
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy import ndimage

# ==================== CONFIGURACION ====================

# --- Geometria del corte ---
GROSOR_CORTE_MM        = 1.2    # ancho REAL de la ranura. [B1] El radio de
#   dilatacion es GROSOR/2. Con spacing 0.45 da 3 voxeles (1.35mm). El ancho
#   real siempre cae en un numero impar de voxeles; cortar() imprime cual salio.

# --- Direccion: SIEMPRE perpendicular a la superficie ---
MODO_VECINDAD          = "angular"   # "angular" | "geodesico"
RADIO_NORMAL_MM        = 2.0
ANGULO_MAX_VECINDAD    = 60.0   # solo entran al promedio vertices de la MISMA
#   cara. Bloquea la tabla interna (~180 grados), las trabeculas de los poros y
#   la placa vecina cruzando una sutura.
MIN_VECINOS            = 3
REFINAR_VECINDAD       = True
ANGULO_SALTO_MAX       = 45.0   # un salto mayor entre normales consecutivas es
#   ruido: no se propaga ni contamina el resto de la curva.
VENTANA_SUAVIZADO_MM   = 1.5

# --- [V1] Cuando se considera que el rayo SALIO del hueso ---
UMBRAL_ESPACIO_ABIERTO_MM = 3.0 # test 3D de burbuja. El hueco tiene que
#   sobrevivir a una erosion de este radio para contar como salida. Los poros y
#   los espacios medulares del diploe (2-4mm) no sobreviven; la cavidad
#   craneal, la orbita y el exterior si.
#   Subir a 4-5 si el corte frena dentro de hueso muy poroso.
#   Bajar a 2 si el corte se pasa de largo hacia estructuras vecinas.
MARGEN_SALIDA_MM       = 2.0    # ademas, desde el punto de salida tiene que
#   haber esta cantidad de aire continuo a lo largo del rayo. Confirmacion de
#   que no caimos en una burbuja de hueso.
GRACIA_INICIO_MM       = 1.0    # tramo inicial exento (la malla y el borde de
#   los voxeles no coinciden exactamente: sin esto el rayo "saldria" en el
#   paso cero).
PROFUNDIDAD_MAX_MM     = 25.0   # tope duro de busqueda por lado.

# --- [V2] Que hacer con un rayo que NO atraviesa ---
TOPE_SIN_SALIDA_FACTOR = 1.5    # se lo corta en seco a este factor por el
#   espesor tipico de la misma curva. Con 1.5 no puede hacer dano lejos del
#   sitio. Subir a 2-3 solo si sabes que ahi el hueso es genuinamente mas
#   grueso que en el resto de la curva.
TOPE_SIN_SALIDA_MIN_MM = 6.0    # piso, para curvas sobre hueso muy fino.
AVISO_SI_NO_ATRAVIESA  = 0.02   # avisar si mas de este porcentaje de puntos de
#   una curva no atraviesa (0.02 = 2%).

# --- [C1] BUSQUEDA DEL CAMINO MAS CORTO (coherente por programacion dinamica) ---
BUSCAR_CAMINO_MAS_CORTO = True  # False = cortar siempre perpendicular (v16)
ANGULO_BUSQUEDA_MAX    = 45.0   # cuanto puede ladearse el corte respecto de la
#   perpendicular, EN EL PLANO PERPENDICULAR A LA CURVA (la unica libertad es
#   inclinarse hacia un lado u otro en la seccion transversal de la placa, que
#   es el mismo grado de libertad que tiene la sierra en la mano).
N_ANGULOS_BUSQUEDA     = 25     # resolucion: con +-45 y 25 angulos, 3.75 grados
MAX_GIRO_POR_PASO_DEG  = 4.0    # <<< EL PARAMETRO QUE MATA LOS TAJOS RANDOM.
#   Giro maximo del corte de un punto de la curva al siguiente (~0.37mm). Con
#   4 grados por paso la direccion puede recorrer todo el rango a lo largo de
#   la curva, pero NO puede pegar un salto. Bajarlo a 2 = cortes mas rigidos y
#   seguros; subirlo a 8-10 = mas libertad y mas riesgo de abanico.
PESO_REGULARIZACION    = 0.35   # entre dos direcciones de recorrido parecido
#   gana la mas perpendicular. Subir = mas conservador.
PESO_SUAVIDAD          = 0.05   # costo (mm) por cada grado de giro. Ademas del
#   limite duro, penaliza girar sin necesidad.
ESPESOR_MIN_VALIDO     = 1.0    # una direccion que "atraviesa" menos que esto
#   esta pasando de refilon por el borde, no cruzando la placa.
PENAL_SIN_SALIDA_MM    = 60.0   # costo de un punto donde ninguna direccion
#   atraviesa. Alto pero FINITO: la secuencia puede pasar por ahi sin romperse.
PASO_BUSQUEDA_VOXELES  = 1.0    # la busqueda usa paso grueso (el corte final
#   sigue usando medio voxel).
BLOQUE_BUSQUEDA        = 512

# --- [B2] Continuidad de la cortina ---
OBJETIVO_CONTINUIDAD   = 0.5    # separacion maxima entre rayos vecinos, en
#   voxeles, a CUALQUIER profundidad. Con 0.5 la cortina no puede perforarse.
MAX_FACTOR_SUPERSAMPLEO = 12

# --- Muestreo (derivado del spacing) ---
FACTOR_MUESTREO        = 0.9
FACTOR_PASO_BARRIDO    = 0.5
DIST_MUESTREO_MM       = 0.0    # 0 = derivar del spacing

# --- Uniones y extensiones de puntas ---
EXTENSION_EXTREMOS_MM  = 3.0    # extension de una punta LIBRE (que no toca
#   ninguna otra curva).
UNION_TOLERANCIA_MM    = 6.0    # [C2] si la punta esta a menos de esto de otra
#   curva, se considera que ahi va una UNION y se extiende EXACTAMENTE hasta
#   tocarla (no se suprime la extension, que es lo que dejaba el hueco).
SOLAPE_UNION_MM        = 1.0    # [C2] cuanto se pasa de largo en la union para
#   garantizar que cierre. Chico a proposito: es un solape, no un espolon.
PROYECTAR_EXTENSION    = True

# --- Verificacion ---
VERIFICAR_SEPARACION   = True   # comprobar, curva por curva, que el corte
#   realmente partio el hueso. Cuesta un etiquetado extra por curva.
ESPESOR_MAX_ESPERADO   = 14.0   # aviso de recorrido sospechosamente largo.

# --- Higiene de escena ---
PREFIJO_CURVAS         = ""

# --- Nodos ---
NOMBRE_MODELO_CRANEO   = "Craneo_Final"
NOMBRE_SEGMENTACION    = "Craneo_Automatico"
NOMBRE_SEGMENTO_HUESO  = "Craneo_Final"
NOMBRE_SEG_TRABAJO     = "CranioPlan_Corte"
NOMBRE_CARPETA_SH      = "CranioPlan_Fragmentos"
NOMBRE_MODELO_PICOS    = "CranioPlan_Picos"
NOMBRE_FID_SIN_PASAR   = "CranioPlan_NoAtraviesa"

# --- Clasificacion y filtrado ---
RUIDO_VOXELES          = 30
VOLUMEN_MINIMO_CM3     = 0.10
FRACCION_INTERIOR      = 0.5
MAX_CAVIDAD_MM3        = 2000.0
SUAVIZADO_EXPORT       = "0.0"
UMBRAL_LIMPIEZA_MALLA  = 0.10


# ============================================================
# PARTE 0a - EL VOLUMEN DE HUESO
# ============================================================
class VolumenHueso:
    """
    Mascara de hueso + transformada RAS->IJK + mascara de ESPACIO ABIERTO.

    ESPACIO ABIERTO [V1a]: aire que sobrevive a una erosion de
    UMBRAL_ESPACIO_ABIERTO_MM. Es el test de "burbuja" en 3D. La pregunta no es
    cuanto aire hay a lo largo del rayo (eso no distingue un espacio medular de
    3mm de la cavidad craneal, porque los dos son aire) sino QUE TAN GRANDE es
    el hueco en volumen.

    Se usa erosion binaria y no transformada de distancia porque la EDT de un
    volumen craneal completo son cientos de MB en float64; la erosion trabaja
    solo con arrays booleanos.
    """

    def __init__(self, boneMask, volumeNode):
        self.mask = boneMask
        self.dims = boneMask.shape                      # (z, y, x)
        self.spacing = tuple(volumeNode.GetSpacing())   # (x, y, z)
        self.sampling = (self.spacing[2], self.spacing[1], self.spacing[0])

        m = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(m)
        self.M = np.array([[m.GetElement(r, c) for c in range(4)]
                           for r in range(4)], dtype=float)

        self.voxel_mm3 = self.spacing[0] * self.spacing[1] * self.spacing[2]
        self.paso = max(min(self.spacing) * FACTOR_PASO_BARRIDO, 0.05)
        self.abierto = self._espacio_abierto()

    def _espacio_abierto(self):
        aire = ~self.mask
        n = int(round(UMBRAL_ESPACIO_ABIERTO_MM / min(self.spacing)))
        if n < 1:
            return aire
        cruz = ndimage.generate_binary_structure(3, 1)
        # border_value=1: el aire de AFUERA del volumen cuenta como abierto,
        # asi el exterior del craneo no se erosiona desde el borde del array.
        return ndimage.binary_erosion(aire, structure=cruz, iterations=n,
                                      border_value=1)

    def muestrear(self, pts):
        h = np.empty((pts.shape[0], 4), dtype=float)
        h[:, :3] = pts
        h[:, 3] = 1.0
        ijk = np.rint(h @ self.M.T)[:, :3].astype(np.int64)
        ii, jj, kk = ijk[:, 0], ijk[:, 1], ijk[:, 2]
        d = self.dims
        dentro = ((kk >= 0) & (kk < d[0]) &
                  (jj >= 0) & (jj < d[1]) &
                  (ii >= 0) & (ii < d[2]))
        esHueso = np.zeros(pts.shape[0], dtype=bool)
        esAbierto = np.ones(pts.shape[0], dtype=bool)   # fuera del volumen = abierto
        esHueso[dentro] = self.mask[kk[dentro], jj[dentro], ii[dentro]]
        esAbierto[dentro] = self.abierto[kk[dentro], jj[dentro], ii[dentro]]
        return ijk, dentro, esHueso, esAbierto


# ============================================================
# PARTE 0b - CAMPO DE NORMALES (la superficie, y SOLO la superficie)
# ============================================================
class CampoNormales:
    """
    El vecindario del que se promedia tiene que estar sobre LA MISMA CARA de
    hueso donde dibujaste la curva: ni la tabla interna, ni los poros, ni la
    placa de al lado. Eso lo garantiza el filtro angular.
    """

    def __init__(self, mallaHueso):
        nf = vtk.vtkPolyDataNormals()
        nf.SetInputData(mallaHueso)
        nf.ComputePointNormalsOn()
        nf.ComputeCellNormalsOff()
        nf.SplittingOff()
        nf.ConsistencyOn()
        nf.AutoOrientNormalsOn()
        nf.Update()
        self.malla = nf.GetOutput()

        arr = self.malla.GetPointData().GetNormals()
        if arr is None:
            raise RuntimeError("la malla no tiene normales calculables")
        self.normales = vtk_to_numpy(arr).astype(float)
        norm = np.linalg.norm(self.normales, axis=1)
        norm[norm < 1e-12] = 1.0
        self.normales /= norm[:, None]

        self.puntos = vtk_to_numpy(self.malla.GetPoints().GetData()).astype(float)

        self.locator = vtk.vtkPointLocator()
        self.locator.SetDataSet(self.malla)
        self.locator.BuildLocator()

        self.cellLocator = vtk.vtkCellLocator()
        self.cellLocator.SetDataSet(self.malla)
        self.cellLocator.BuildLocator()

        self._adyacencia = None
        self.cosVecindad = float(np.cos(np.radians(ANGULO_MAX_VECINDAD)))

    def _construir_adyacencia(self):
        polys = self.malla.GetPolys()
        if polys is None or polys.GetNumberOfCells() == 0:
            self._adyacencia = None
            return
        data = vtk_to_numpy(polys.GetData())
        if data.size % 4 != 0 or not np.all(data[0::4] == 3):
            print("  AVISO: la malla no es puramente de triangulos; "
                  "uso MODO_VECINDAD='angular'.")
            self._adyacencia = None
            return
        tri = data.reshape(-1, 4)[:, 1:].astype(np.int64)
        e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        e = np.vstack([e, e[:, ::-1]])
        nPts = self.puntos.shape[0]
        e = e[np.argsort(e[:, 0], kind='stable')]
        conteo = np.bincount(e[:, 0], minlength=nPts)
        indptr = np.zeros(nPts + 1, np.int64)
        np.cumsum(conteo, out=indptr[1:])
        self._adyacencia = (indptr, e[:, 1].copy())

    def _vecinos_geodesicos(self, idSemilla, radio):
        import heapq
        if self._adyacencia is None:
            self._construir_adyacencia()
        if self._adyacencia is None:
            return None
        indptr, indices = self._adyacencia
        dist = {int(idSemilla): 0.0}
        cola = [(0.0, int(idSemilla))]
        while cola:
            d, v = heapq.heappop(cola)
            if d > dist.get(v, np.inf) + 1e-12:
                continue
            for k in range(indptr[v], indptr[v + 1]):
                w = int(indices[k])
                nd = d + float(np.linalg.norm(self.puntos[w] - self.puntos[v]))
                if nd <= radio and nd < dist.get(w, np.inf):
                    dist[w] = nd
                    heapq.heappush(cola, (nd, w))
        return np.fromiter(dist.keys(), dtype=np.int64, count=len(dist))

    def proyectar(self, punto):
        cerrado = [0.0, 0.0, 0.0]
        cellId = vtk.reference(0)
        subId = vtk.reference(0)
        d2 = vtk.reference(0.0)
        self.cellLocator.FindClosestPoint(punto.tolist(), cerrado, cellId, subId, d2)
        return np.array(cerrado, dtype=float)

    def normal_en(self, punto):
        """Devuelve (normal, calidad, nUsados)."""
        idCerca = self.locator.FindClosestPoint(punto.tolist())
        if idCerca < 0:
            return None, 0.0, 0
        nRef = self.normales[idCerca]

        if MODO_VECINDAD == "geodesico":
            idx = self._vecinos_geodesicos(idCerca, RADIO_NORMAL_MM)
            if idx is None or idx.size == 0:
                idx = np.array([idCerca], dtype=np.int64)
        else:
            ids = vtk.vtkIdList()
            self.locator.FindPointsWithinRadius(RADIO_NORMAL_MM, punto.tolist(), ids)
            m = ids.GetNumberOfIds()
            if m == 0:
                idx = np.array([idCerca], dtype=np.int64)
            else:
                idx = np.fromiter((ids.GetId(j) for j in range(m)),
                                  dtype=np.int64, count=m)

        cand = self.normales[idx]
        sel = (cand @ nRef) >= self.cosVecindad
        usados = cand[sel] if np.any(sel) else nRef[None, :]

        acum = usados.sum(axis=0)
        na = np.linalg.norm(acum)
        if na < 1e-9:
            return nRef.copy(), 0.0, int(usados.shape[0])
        prom = acum / na

        if REFINAR_VECINDAD and usados.shape[0] >= MIN_VECINOS:
            sel2 = (cand @ prom) >= self.cosVecindad
            if np.any(sel2):
                usados = cand[sel2]
                acum = usados.sum(axis=0)
                na = np.linalg.norm(acum)
                if na > 1e-9:
                    prom = acum / na

        nUsados = int(usados.shape[0])
        return prom, float(na / max(nUsados, 1)), nUsados


# ============================================================
# PARTE 1 - GEOMETRIA DE LA CURVA
# ============================================================
def _resamplear_puntos(puntosVtk, distanciaMM):
    """
    Re-distribuye la polilinea densa de Slicer en pasos parejos. Tiene que ser
    fino porque la cortina se estampa punto por punto.
    """
    n = puntosVtk.GetNumberOfPoints()
    if n < 2:
        return []
    original = [np.array(puntosVtk.GetPoint(i)) for i in range(n)]
    resampleados = [original[0]]
    acumulada = 0.0
    anterior = original[0]
    for i in range(1, n):
        actual = original[i]
        segmento = actual - anterior
        largo = np.linalg.norm(segmento)
        while acumulada + largo >= distanciaMM:
            falta = distanciaMM - acumulada
            direccion = segmento / largo if largo > 1e-9 else segmento
            nuevo = anterior + direccion * falta
            resampleados.append(nuevo)
            anterior = nuevo
            segmento = actual - anterior
            largo = np.linalg.norm(segmento)
            acumulada = 0.0
        acumulada += largo
        anterior = actual
    if np.linalg.norm(resampleados[-1] - original[-1]) > 1e-6:
        resampleados.append(original[-1])
    return resampleados


def _quitar_coincidentes(posiciones, esCerrada, toleranciaMM=0.2):
    limpias = []
    for p in posiciones:
        if not limpias or np.linalg.norm(p - limpias[-1]) > toleranciaMM:
            limpias.append(p)
    if esCerrada and len(limpias) > 2:
        if np.linalg.norm(limpias[-1] - limpias[0]) <= toleranciaMM:
            limpias.pop()
    return limpias


def _extender_extremos(posiciones, largoIni, largoFin, pasoMM, campo):
    """
    [C2] Prolonga cada punta el largo que le corresponde, proyectando sobre la
    superficie para que la extension SIGA el hueso en vez de volar por la
    tangente.

    El largo ya viene calculado por _largo_extension_puntas(): para una punta
    LIBRE es EXTENSION_EXTREMOS_MM; para una punta que va a UNIRSE con otra
    curva es exactamente la distancia que falta para tocarla, mas un solape
    chico.

    v16 hacia otra cosa y estaba mal: cuando detectaba una union NO extendia
    nada, dando por sentado que las curvas ya se tocaban. Si quedaban a 2mm,
    esos 2mm de hueso jamas se cortaban, y en un craneo cerrado un solo puente
    de 2mm en una sien deja el colgajo colgado y no se separa nada.
    """
    if len(posiciones) < 2:
        return list(posiciones), np.zeros(len(posiciones), dtype=bool)

    def _muestras(desde, direccion, largo):
        out = []
        t = pasoMM
        while t <= largo + 1e-9:
            p = desde + direccion * t
            if PROYECTAR_EXTENSION and campo is not None:
                p = campo.proyectar(p)
            out.append(p)
            t += pasoMM
        return out

    extIni, extFin = [], []
    if largoIni > 0:
        d = posiciones[0] - posiciones[1]
        nd = np.linalg.norm(d)
        if nd > 1e-9:
            extIni = _muestras(posiciones[0], d / nd, largoIni)
            extIni.reverse()
    if largoFin > 0:
        d = posiciones[-1] - posiciones[-2]
        nd = np.linalg.norm(d)
        if nd > 1e-9:
            extFin = _muestras(posiciones[-1], d / nd, largoFin)

    todas = extIni + list(posiciones) + extFin
    esExt = np.zeros(len(todas), dtype=bool)
    if extIni:
        esExt[:len(extIni)] = True
    if extFin:
        esExt[len(todas) - len(extFin):] = True
    return todas, esExt


def _puntos_de_curva(curvaNode):
    p = curvaNode.GetCurvePointsWorld()
    if p is None or p.GetNumberOfPoints() == 0:
        return np.zeros((0, 3))
    return vtk_to_numpy(p.GetData()).astype(float)


def _largo_extension_puntas(curvaNode, otrasCurvas):
    """
    [C2] Cuanto hay que prolongar cada punta.

      - punta LIBRE (lejos de toda otra curva): EXTENSION_EXTREMOS_MM, para
        asegurar que llegue al borde del hueso.
      - punta que va a UNIRSE con otra curva (a menos de UNION_TOLERANCIA_MM):
        exactamente la distancia que falta para tocarla, mas SOLAPE_UNION_MM.
        Ni menos (quedaria un puente de hueso sin cortar en la union, y con eso
        no se separa nada) ni de mas (seria un espolon pasandose de la esquina).

    Devuelve (largoIni, largoFin, distIni, distFin) en mm; las distancias son
    para el reporte.
    """
    pts = _puntos_de_curva(curvaNode)
    if pts.shape[0] < 2:
        return 0.0, 0.0, np.inf, np.inf
    otros = [_puntos_de_curva(c) for c in otrasCurvas]
    otros = [o for o in otros if o.shape[0] > 0]
    if not otros:
        return (EXTENSION_EXTREMOS_MM, EXTENSION_EXTREMOS_MM, np.inf, np.inf)
    todos = np.vstack(otros)

    def _largo(p):
        d = float(np.linalg.norm(todos - p[None, :], axis=1).min())
        if d <= UNION_TOLERANCIA_MM:
            return d + SOLAPE_UNION_MM, d
        return EXTENSION_EXTREMOS_MM, d

    lIni, dIni = _largo(pts[0])
    lFin, dFin = _largo(pts[-1])
    return lIni, lFin, dIni, dFin


def _tangentes(pos, esCerrada):
    """Tangente unitaria de la polilinea en cada punto."""
    n = pos.shape[0]
    if n < 2:
        return np.tile(np.array([1.0, 0.0, 0.0]), (max(n, 1), 1))
    if esCerrada:
        t = np.roll(pos, -1, axis=0) - np.roll(pos, 1, axis=0)
    else:
        t = np.empty_like(pos)
        t[1:-1] = pos[2:] - pos[:-2]
        t[0] = pos[1] - pos[0]
        t[-1] = pos[-1] - pos[-2]
    nn = np.linalg.norm(t, axis=1)
    nn[nn < 1e-12] = 1.0
    return t / nn[:, None]


def _base_plano_normal(normales, tangentes):
    """
    Base ortonormal (n, b) del PLANO PERPENDICULAR A LA CURVA. La busqueda vive
    solo aca dentro: la unica libertad es cuanto se ladea el corte en la
    seccion transversal de la placa. Asi el corte no puede "descuadrarse" a lo
    largo de la curva.
    """
    proy = np.sum(normales * tangentes, axis=1)
    n = normales - tangentes * proy[:, None]
    nn = np.linalg.norm(n, axis=1)
    malos = nn < 1e-6
    if np.any(malos):
        n[malos] = normales[malos]
        nn = np.linalg.norm(n, axis=1)
    nn[nn < 1e-12] = 1.0
    n = n / nn[:, None]
    b = np.cross(tangentes, n)
    bn = np.linalg.norm(b, axis=1)
    bn[bn < 1e-12] = 1.0
    return n, b / bn[:, None]


# ============================================================
# PARTE 2 - DIRECCIONES (siempre perpendiculares)
# ============================================================
def _direcciones_de_curva(campo, posiciones, esExtension, esCerrada, distMuestreo):
    """
    Normal de superficie por punto + orientacion de signos con guarda + suavizado.
    Un solo punto malo NO puede invertir el signo del resto de la curva.
    """
    n = len(posiciones)
    crudas = [None] * n
    calidades = np.zeros(n)
    buenos = np.zeros(n, dtype=bool)

    for i, p in enumerate(posiciones):
        nr, cal, nu = campo.normal_en(p)
        crudas[i] = nr
        calidades[i] = cal
        buenos[i] = (nr is not None) and (nu >= MIN_VECINOS)

    idxValidos = [i for i in range(n) if crudas[i] is not None]
    if not idxValidos:
        return None, {"error": "ningun punto tiene superficie cerca"}
    for i in range(n):
        if crudas[i] is None:
            j = min(idxValidos, key=lambda k: abs(k - i))
            crudas[i] = np.array(crudas[j])

    candidatos = [i for i in range(n) if buenos[i] and not esExtension[i]]
    if not candidatos:
        candidatos = [i for i in range(n) if not esExtension[i]]
    if not candidatos:
        candidatos = idxValidos
    semilla = max(candidatos, key=lambda i: (calidades[i], -abs(i - n // 2)))

    normales = [np.array(v) for v in crudas]
    cosSalto = float(np.cos(np.radians(ANGULO_SALTO_MAX)))
    nSospechosos = 0

    def _propagar(rango):
        nonlocal nSospechosos
        ref = normales[semilla].copy()
        for i in rango:
            ni = normales[i]
            if not buenos[i]:
                normales[i] = ref.copy()
                continue
            if float(np.dot(ni, ref)) < 0.0:
                ni = -ni
            if float(np.dot(ni, ref)) < cosSalto:
                nSospechosos += 1
                normales[i] = ref.copy()
                continue
            normales[i] = ni
            ref = ni

    _propagar(range(semilla + 1, n))
    _propagar(range(semilla - 1, -1, -1))

    N = np.array(normales, dtype=float)

    circular = esCerrada
    if esCerrada and n > 2 and float(np.dot(N[-1], N[0])) < 0.0:
        circular = False
        print("    AVISO: el lazo no cierra con signo consistente; "
              "desactivo el suavizado circular.")

    w = max(1, int(round(VENTANA_SUAVIZADO_MM / max(distMuestreo, 1e-6))))
    tam = 2 * w + 1
    if n > tam:
        S = ndimage.uniform_filter1d(N, size=tam, axis=0,
                                     mode='wrap' if circular else 'nearest')
    else:
        S = N.copy()

    norm = np.linalg.norm(S, axis=1)
    flojas = norm < 1e-6
    if np.any(flojas):
        S[flojas] = N[flojas]
        norm = np.linalg.norm(S, axis=1)
    norm[norm < 1e-12] = 1.0
    S = S / norm[:, None]

    return S, {"n": n, "sospechosos": nSospechosos,
               "bajaCalidad": int(np.count_nonzero(~buenos)),
               "calidadMedia": float(np.mean(calidades)),
               "ventanaPuntos": w}


# ============================================================
# PARTE 3 - BARRIDO CON VERIFICACION [V1][V2]
# ============================================================
def _indice_de_salida(esHueso, esAbierto, g, mMargen):
    """
    [V1] Donde -y si- el rayo SALE del hueso.

    Se acepta la salida solo si se cumplen las DOS condiciones:
      (a) la muestra esta en ESPACIO ABIERTO (test 3D de burbuja: el hueco
          sobrevivio a la erosion, o sea no es un poro ni un espacio medular);
      (b) desde ahi hay mMargen muestras seguidas SIN hueso (el margen de
          seguridad a lo largo del rayo).

    Cualquiera de las dos sola se deja enganar en algun caso: el margen a lo
    largo del rayo no distingue un espacio medular de 3mm de la cavidad, y el
    test 3D podria dar un falso positivo justo en el borde de una cavidad. Las
    dos juntas son dificiles de enganar.

    Devuelve (corte, recorrido, atraviesa):
      corte     = indice de la salida (o T si nunca salio)
      recorrido = indice del ultimo hueso antes de la salida (espesor real)
      atraviesa = si de verdad salio del hueso
    """
    P, T = esHueso.shape

    # (b) ventana deslizante sin hueso
    if T > mMargen >= 1:
        c = np.concatenate([np.zeros((P, 1), np.int32),
                            np.cumsum(esHueso.astype(np.int32), axis=1)], axis=1)
        win = c[:, mMargen:] - c[:, :T - mMargen + 1]
        libre = np.zeros((P, T), dtype=bool)
        libre[:, :T - mMargen + 1] = (win == 0)
    else:
        libre = ~esHueso

    # (a) AND (b)
    valido = esAbierto & libre
    if g < T:
        valido[:, :g] = False
    else:
        valido[:] = False

    atraviesa = valido.any(axis=1)
    corte = np.where(atraviesa, np.argmax(valido, axis=1), T)

    idx = np.arange(T)[None, :]
    Hm = esHueso & (idx < corte[:, None])
    hay = Hm.any(axis=1)
    recorrido = np.where(hay, T - np.argmax(Hm[:, ::-1], axis=1), 0)
    return corte, recorrido, atraviesa


def _barrido(vol, pos, dirs, topeSinSalidaMM=None):
    """
    Estampa la cortina y verifica rayo por rayo.

    [V2] Si un rayo NO atraviesa, se lo corta en seco a topeSinSalidaMM en vez
    de dejarlo barrer hasta PROFUNDIDAD_MAX_MM. Eso es lo que impide que una
    direccion mala deje una lengua de corte de varios centimetros en diagonal
    (los "tajos random" de v14/v15).

    Devuelve (curtain, recorrido (N,2), atraviesa (N,2)).
    """
    paso = vol.paso
    N = pos.shape[0]
    ts = np.arange(0.0, PROFUNDIDAD_MAX_MM + paso, paso)
    T = int(ts.size)
    g = max(1, int(round(GRACIA_INICIO_MM / paso)))
    mMargen = max(1, int(round(MARGEN_SALIDA_MM / paso)))
    margenKerf = max(1, int(round(GROSOR_CORTE_MM / paso)))
    topeIdx = T if topeSinSalidaMM is None else max(
        1, min(int(round(topeSinSalidaMM / paso)), T))

    curtain = np.zeros(vol.dims, dtype=bool)
    recorrido = np.zeros((N, 2), dtype=float)
    atraviesa = np.zeros((N, 2), dtype=bool)

    for lado, signo in enumerate((1.0, -1.0)):
        m = (pos[:, None, :] +
             dirs[:, None, :] * (signo * ts)[None, :, None]).reshape(-1, 3)
        ijk, dentro, esHueso, esAbierto = vol.muestrear(m)
        H = esHueso.reshape(N, T)
        A = esAbierto.reshape(N, T)
        corte, rec, atr = _indice_de_salida(H, A, g, mMargen)

        recorrido[:, lado] = rec * paso
        atraviesa[:, lado] = atr

        # se estampa un poco mas alla de la salida para que el kerf quede bien
        # formado justo en la superficie de salida (ese sobrante cae en aire y
        # lo descarta el AND con el hueso)
        hasta = np.where(atr, np.minimum(corte + margenKerf, T), topeIdx)
        keep = (np.arange(T)[None, :] < hasta[:, None]).reshape(-1)
        sel = keep & dentro
        if np.any(sel):
            curtain[ijk[sel, 2], ijk[sel, 1], ijk[sel, 0]] = True

    return curtain, recorrido, atraviesa


def _medir_cruce(vol, pos, dirs, paso, g, mMargen):
    """
    Recorrido GEOMETRICO hasta salir del hueso (los dos lados sumados) y si
    logro salir de ambos.

    Es recorrido, NO cantidad de hueso: en diploe poroso, minimizar hueso
    premiaba a las direcciones oblicuas que enhebran los espacios medulares,
    que cruzan menos hueso pero recorren mucho mas camino. El recorrido es el
    espesor real de la placa en esa direccion y los poros no lo engañan.
    """
    P = pos.shape[0]
    ts = np.arange(0.0, PROFUNDIDAD_MAX_MM + paso, paso)
    T = int(ts.size)
    rec = np.zeros(P, dtype=float)
    sal = np.ones(P, dtype=bool)
    for signo in (1.0, -1.0):
        m = (pos[:, None, :] +
             dirs[:, None, :] * (signo * ts)[None, :, None]).reshape(-1, 3)
        _, _, H, A = vol.muestrear(m)
        _, r, s = _indice_de_salida(H.reshape(P, T), A.reshape(P, T), g, mMargen)
        rec += r * paso
        sal &= s
    return rec, sal


def _viterbi_angulos(costo, resDeg, esCerrada):
    """
    [C1] Elige la SECUENCIA de angulos de costo total minimo, con la
    restriccion dura de no girar mas de MAX_GIRO_POR_PASO_DEG entre puntos
    consecutivos. Programacion dinamica exacta, O(N x K^2).

    Esta es la pieza que hace imposible el tajo random. Elegir el minimo en
    cada punto POR SEPARADO (lo que hacian v14/v15) permite que dos vecinos
    caigan en minimos opuestos; aca la transicion prohibida tiene costo
    infinito, asi que esa secuencia directamente no existe en el espacio de
    soluciones.

    costo: (N, K). Devuelve el indice de angulo elegido por punto, (N,).
    """
    N, K = costo.shape
    if N == 0:
        return np.zeros(0, dtype=np.int32)
    if N == 1:
        return np.array([int(np.argmin(costo[0]))], dtype=np.int32)

    dk = max(1, int(round(MAX_GIRO_POR_PASO_DEG / max(resDeg, 1e-6))))
    salto = np.abs(np.arange(K)[:, None] - np.arange(K)[None, :])
    INF = 1e18
    trans = np.where(salto <= dk, PESO_SUAVIDAD * salto * resDeg, INF)

    # en curvas cerradas la costura se pone donde la eleccion es menos ambigua
    off = int(np.argmin(costo.min(axis=1))) if esCerrada and N > 2 else 0
    C = np.roll(costo, -off, axis=0) if off else costo

    acum = C[0].copy()
    prev = np.zeros((N, K), dtype=np.int32)
    cols = np.arange(K)
    for i in range(1, N):
        tot = acum[:, None] + trans          # (K anterior, K actual)
        pk = np.argmin(tot, axis=0)
        acum = tot[pk, cols] + C[i]
        prev[i] = pk

    idx = np.zeros(N, dtype=np.int32)
    idx[-1] = int(np.argmin(acum))
    for i in range(N - 1, 0, -1):
        idx[i - 1] = prev[i, idx[i]]
    return np.roll(idx, off) if off else idx


def _buscar_camino_mas_corto(vol, pos, normales, tangentes, esCerrada):
    """
    [C1] EL cambio de v17.

    Busca, para toda la curva a la vez, la SECUENCIA de direcciones que
    atraviesa el hueso por el camino mas corto, con la restriccion dura de que
    de un punto al siguiente el angulo no pueda girar mas de
    MAX_GIRO_POR_PASO_DEG.

    Por que asi y no punto por punto (que es lo que hacian v14/v15):
      Si en un punto el mejor camino es +40 grados y en el vecino es -40 (dos
      minimos casi empatados, comun en geometria complicada), eligiendo por
      separado quedan dos rayos VECINOS apuntando a lados opuestos. Eso es el
      tajo random. Y suavizar despues no arregla nada: el promedio de +40 y -40
      es 0, que no es ninguna de las dos y suele ser la direccion rasante.

      Con la restriccion de giro, esa situacion es IMPOSIBLE POR CONSTRUCCION.
      Y como se optimiza la suma de toda la curva, la solucion se "compromete":
      si al conjunto de la ceja le conviene resolverlo yendo hacia arriba, van
      todos hacia arriba; ninguno se escapa hacia la orbita porque ahi seria un
      poquito mas corto, porque cambiar de bando cuesta caro en el total.

    Se resuelve exacto con programacion dinamica (Viterbi) sobre la grilla de
    angulos: costo O(N x K^2), milisegundos.

    Devuelve (dirs, angulos, atraviesaElegido, diag).
    """
    N = pos.shape[0]
    n, b = _base_plano_normal(normales, tangentes)

    K = max(3, int(N_ANGULOS_BUSQUEDA))
    thetas = np.radians(np.linspace(-ANGULO_BUSQUEDA_MAX, ANGULO_BUSQUEDA_MAX, K))
    resDeg = (2.0 * ANGULO_BUSQUEDA_MAX) / (K - 1)

    paso = max(min(vol.spacing) * PASO_BUSQUEDA_VOXELES, 0.05)
    g = max(1, int(round(GRACIA_INICIO_MM / paso)))
    mMargen = max(1, int(round(MARGEN_SALIDA_MM / paso)))

    # --- costo de cada (punto, angulo) ---
    costo = np.empty((N, K), dtype=float)
    recor = np.empty((N, K), dtype=float)
    pasaKK = np.zeros((N, K), dtype=bool)
    for k, th in enumerate(thetas):
        d = np.cos(th) * n + np.sin(th) * b
        r = np.empty(N)
        s = np.empty(N, dtype=bool)
        for ini in range(0, N, BLOQUE_BUSQUEDA):
            fin = min(ini + BLOQUE_BUSQUEDA, N)
            rr, ss = _medir_cruce(vol, pos[ini:fin], d[ini:fin], paso, g, mMargen)
            r[ini:fin] = rr
            s[ini:fin] = ss
        ok = s & (r >= ESPESOR_MIN_VALIDO)
        # entre dos direcciones de recorrido parecido gana la mas perpendicular
        costo[:, k] = np.where(ok,
                               r * (1.0 + PESO_REGULARIZACION * (1.0 - np.cos(th))),
                               PENAL_SIN_SALIDA_MM)
        recor[:, k] = r
        pasaKK[:, k] = ok

    idx = _viterbi_angulos(costo, resDeg, esCerrada)
    ang = thetas[idx]
    fila = np.arange(N)
    pasa = pasaKK[fila, idx]
    rec = recor[fila, idx]

    dirs = np.cos(ang)[:, None] * n + np.sin(ang)[:, None] * b
    nn = np.linalg.norm(dirs, axis=1)
    nn[nn < 1e-12] = 1.0
    dirs = dirs / nn[:, None]

    grados = np.degrees(ang)
    giro = np.abs(np.diff(grados)) if N > 1 else np.zeros(1)
    diag = {
        "inclinacionMediana": float(np.median(np.abs(grados))),
        "inclinacionMax": float(np.max(np.abs(grados))),
        "giroMaxPorPaso": float(np.max(giro)) if giro.size else 0.0,
        "resolucionAngular": resDeg,
        "recorridoConBusqueda": float(np.median(rec[pasa])) if np.any(pasa) else 0.0,
    }
    return dirs, ang, pasa, diag


def _tirada_mas_larga(pasa, esCerrada):
    """
    [C3] Tirada CONTINUA mas larga de puntos que no atraviesan, en numero de
    puntos.

    Este es el numero que decide si el hueso se separa o no. El porcentaje
    engaña: la separacion no es un promedio, es una CADENA. Cortar el 95% de
    una cinta y dejar el 5% intacto no deja la cinta "95% cortada": la deja
    entera. 16 puntos seguidos a 0.37mm son ~6mm de hueso macizo que sostienen
    todo el colgajo.
    """
    mal = ~np.asarray(pasa, dtype=bool)
    if not mal.any():
        return 0
    x = np.concatenate([mal, mal]) if esCerrada else mal
    mejor = actual = 0
    for v in x:
        actual = actual + 1 if v else 0
        if actual > mejor:
            mejor = actual
    return int(min(mejor, mal.size))


def _resumen_penetracion(recorrido, atraviesa):
    """
    De los dos lados del barrido, uno sale al aire enseguida y el otro entra al
    hueso. Nos interesa el que entra: cuanto penetro y si logro atravesar.
    """
    lado = np.argmax(recorrido, axis=1)
    fila = np.arange(recorrido.shape[0])
    penetracion = recorrido[fila, lado]
    pasa = atraviesa[fila, lado]
    signo = np.where(lado == 0, 1.0, -1.0)
    return penetracion, pasa, signo


def _supersamplear(pos, dirs, esCerrada, prof, objetivoMM):
    """
    [B2] Los rayos estan a medio voxel en la SUPERFICIE, pero al entrar
    divergen (curvatura del craneo, y variacion de la direccion a lo largo de
    la curva). Si a profundidad se separan mas de un voxel, entre rayo y rayo
    queda hueso sin cortar y la cortina sale perforada: esa era la hilera de
    agujeritos.

    Aca se mide la separacion REAL a la profundidad que se va a usar y se
    insertan rayos interpolados hasta que no puedan separarse mas de
    'objetivoMM' a ninguna profundidad.
    """
    N = pos.shape[0]
    if N < 2:
        return pos, dirs, 1, 0.0
    if esCerrada:
        p1, d1 = pos, dirs
        p2, d2 = np.roll(pos, -1, axis=0), np.roll(dirs, -1, axis=0)
    else:
        p1, d1 = pos[:-1], dirs[:-1]
        p2, d2 = pos[1:], dirs[1:]

    sep = np.zeros(p1.shape[0])
    for t in (0.0, prof, -prof):
        sep = np.maximum(sep, np.linalg.norm((p2 + d2 * t) - (p1 + d1 * t), axis=1))
    sepMax = float(sep.max()) if sep.size else 0.0

    K = 1 if objetivoMM <= 0 else int(np.ceil(sepMax / objetivoMM))
    K = max(1, min(K, MAX_FACTOR_SUPERSAMPLEO))
    if K <= 1:
        return pos, dirs, 1, sepMax

    us = np.linspace(0.0, 1.0, K, endpoint=False)
    pl = (p1[:, None, :] + (p2 - p1)[:, None, :] * us[None, :, None]).reshape(-1, 3)
    dl = (d1[:, None, :] + (d2 - d1)[:, None, :] * us[None, :, None]).reshape(-1, 3)
    nn = np.linalg.norm(dl, axis=1)
    nn[nn < 1e-12] = 1.0
    dl = dl / nn[:, None]
    if not esCerrada:
        pl = np.vstack([pl, pos[-1:]])
        dl = np.vstack([dl, dirs[-1:]])
    return pl, dl, K, sepMax


def _dilatar_a_kerf(curtain, vol, grosorMM):
    """
    [B1] El radio correcto para una ranura de ancho W es W/2.
    v13/v14 usaban (W - spacing)/2: con spacing 0.45 y W=1.2 daba 0.375mm, y
    como el voxel vecino esta a 0.45mm no entraba ninguno. La ranura quedaba de
    UN SOLO VOXEL, mas fina que la de v12.
    """
    radio = max(0.0, grosorMM / 2.0)
    sp = min(vol.spacing)
    if radio < sp * 0.5:
        return curtain, sp
    dist = ndimage.distance_transform_edt(~curtain, sampling=vol.sampling)
    capas = int(np.floor(radio / sp))
    return (dist <= radio), (2 * capas + 1) * sp


# ============================================================
# PARTE 4 - ANILLO DE CORTE DE UNA CURVA
# ============================================================
def _preparar_curva(curvaNode, campo, esCerrada, otrasCurvas, distMuestreo,
                    vol=None):
    """
    Puntos resampleados + extensiones + direcciones. Comun a diagnosticar y
    cortar, para que las dos vean exactamente lo mismo.

    Si se pasa 'vol' y BUSCAR_CAMINO_MAS_CORTO esta activo, la direccion ya
    sale corregida por la busqueda coherente [C1].
    """
    pv = curvaNode.GetCurvePointsWorld()
    if pv is None or pv.GetNumberOfPoints() < 2:
        return None, None, {"error": "curva con menos de 2 puntos"}

    posiciones = _quitar_coincidentes(_resamplear_puntos(pv, distMuestreo), esCerrada)
    posiciones = list(posiciones)
    esExt = np.zeros(len(posiciones), dtype=bool)
    info = {}
    if not esCerrada:
        lIni, lFin, dIni, dFin = _largo_extension_puntas(curvaNode, otrasCurvas)
        posiciones, esExt = _extender_extremos(
            posiciones, lIni, lFin, distMuestreo, campo)
        info = {"largoExtIni": lIni, "largoExtFin": lFin,
                "distUnionIni": dIni, "distUnionFin": dFin}

    if len(posiciones) < 2:
        return None, None, {"error": "curva demasiado corta tras el resampleo"}

    dirs, diag = _direcciones_de_curva(campo, posiciones, esExt, esCerrada, distMuestreo)
    if dirs is None:
        return None, None, diag
    diag.update(info)

    pos = np.asarray(posiciones, dtype=float)
    if vol is not None and BUSCAR_CAMINO_MAS_CORTO:
        tang = _tangentes(pos, esCerrada)
        dirs, ang, pasaB, diagB = _buscar_camino_mas_corto(
            vol, pos, dirs, tang, esCerrada)
        diag.update(diagB)
        diag["busqueda"] = True
    else:
        diag["busqueda"] = False
    return pos, dirs, diag


def _anillo_de_corte(vol, curvaNode, campo, esCerrada, otrasCurvas,
                     distMuestreo, nIslasAntes=None):
    pos, dirs, diag = _preparar_curva(curvaNode, campo, esCerrada,
                                      otrasCurvas, distMuestreo, vol)
    if pos is None:
        return None, diag
    diag["puntosSinPasarXYZ"] = np.zeros((0, 3))

    # --- primera pasada: medir espesor y verificar quien atraviesa ---
    _, rec0, atr0 = _barrido(vol, pos, dirs, None)
    pen0, pasa0, _ = _resumen_penetracion(rec0, atr0)

    diag["nPuntos"] = int(pos.shape[0])
    diag["nNoAtraviesa"] = int(np.count_nonzero(~pasa0))
    diag["fracNoAtraviesa"] = float(np.count_nonzero(~pasa0) / max(pos.shape[0], 1))
    # [C3] lo que de verdad decide la separacion
    tirada = _tirada_mas_larga(pasa0, esCerrada)
    diag["tiradaPuntos"] = tirada
    diag["tiradaMM"] = tirada * distMuestreo
    if np.any(pasa0):
        espTipico = float(np.median(pen0[pasa0]))
        diag["espesorMediano"] = espTipico
        diag["espesorMax"] = float(np.max(pen0[pasa0]))
    else:
        espTipico = TOPE_SIN_SALIDA_MIN_MM
        diag["espesorMediano"] = 0.0
        diag["espesorMax"] = 0.0
    if np.any(~pasa0):
        diag["puntosSinPasarXYZ"] = pos[~pasa0]

    # [V2] tope para los rayos que no atraviesan: no pueden hacer dano lejos
    tope = max(TOPE_SIN_SALIDA_MIN_MM, TOPE_SIN_SALIDA_FACTOR * espTipico)
    diag["topeSinSalida"] = tope

    # --- [B2] continuidad de la cortina a la profundidad real de trabajo ---
    prof = max(float(np.percentile(pen0, 98)) if pen0.size else 5.0, 2.0)
    objetivo = OBJETIVO_CONTINUIDAD * min(vol.spacing)
    pos2, dirs2, K, sepMax = _supersamplear(pos, dirs, esCerrada, prof, objetivo)
    diag["factorSupersampleo"] = K
    diag["separacionMaxima"] = sepMax
    diag["separacionTope"] = bool(K >= MAX_FACTOR_SUPERSAMPLEO and sepMax / K > objetivo)

    curtain, _, _ = _barrido(vol, pos2, dirs2, tope)
    if not curtain.any():
        return None, dict(diag, error="la cortina quedo vacia (la curva no toca hueso?)")

    curtainDil, anchoReal = _dilatar_a_kerf(curtain, vol, GROSOR_CORTE_MM)
    diag["kerfReal"] = anchoReal

    anillo = np.logical_and(curtainDil, vol.mask)
    if not anillo.any():
        return None, dict(diag, error="el anillo quedo vacio tras el AND con el hueso")

    if VERIFICAR_SEPARACION and nIslasAntes is not None:
        _, nDesp = ndimage.label(vol.mask & ~anillo)
        diag["separo"] = bool(nDesp > nIslasAntes)
        diag["islasNuevas"] = int(nDesp - nIslasAntes)

    return anillo, diag


# ============================================================
# PARTE 5 - TAPONES
# ============================================================
def _tapon_de_curva(curvaNode, longitud):
    cap = vtk.vtkPolyData()
    try:
        slicer.modules.markups.logic().GetClosedCurveSurfaceArea(curvaNode, cap)
    except Exception as e:
        print(f"  no pude generar el tapon de '{curvaNode.GetName()}': {e}")
        return None
    if cap.GetNumberOfPoints() < 3:
        return None
    pts = np.array([cap.GetPoint(i) for i in range(cap.GetNumberOfPoints())])
    centro = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - centro)
    normal = vt[2] / np.linalg.norm(vt[2])
    tr = vtk.vtkTransform()
    tr.Translate(*(-normal * longitud / 2.0))
    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetInputData(cap)
    tf.SetTransform(tr)
    tf.Update()
    ext = vtk.vtkLinearExtrusionFilter()
    ext.SetInputData(tf.GetOutput())
    ext.SetExtrusionTypeToVectorExtrusion()
    ext.SetVector(*(normal * longitud))
    ext.CappingOn()
    ext.Update()
    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(ext.GetOutputPort())
    tri.Update()
    limp = vtk.vtkCleanPolyData()
    limp.SetInputConnection(tri.GetOutputPort())
    limp.Update()
    return limp.GetOutput()


def _mascaras_de_tapones(volumeNode, longitud, curvasCerradas):
    if not curvasCerradas:
        return {}
    previo = slicer.mrmlScene.GetFirstNodeByName("CranioPlan_Tapones")
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)
    segTap = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode',
                                                "CranioPlan_Tapones")
    segTap.CreateDefaultDisplayNodes()
    segTap.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    if segTap.GetDisplayNode():
        segTap.GetDisplayNode().SetVisibility(False)

    temporales, mapaNombre = [], {}
    for curva in curvasCerradas:
        poly = _tapon_de_curva(curva, longitud)
        if poly is None:
            continue
        mNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLModelNode", f"_Tapon_{curva.GetName()}")
        mNode.SetAndObservePolyData(poly)
        mNode.CreateDefaultDisplayNodes()
        mNode.GetDisplayNode().SetVisibility(False)
        temporales.append(mNode)
        antes = set(_ids_de_segmentos(segTap))
        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(mNode, segTap)
        for i in [x for x in _ids_de_segmentos(segTap) if x not in antes]:
            mapaNombre[i] = curva.GetName()

    mascaras = {}
    for segId, nombreCurva in mapaNombre.items():
        try:
            mascaras[nombreCurva] = slicer.util.arrayFromSegmentBinaryLabelmap(
                segTap, segId, volumeNode).astype(bool)
        except Exception as e:
            print(f"  no pude rasterizar el tapon de '{nombreCurva}': {e}")
    for m in temporales:
        slicer.mrmlScene.RemoveNode(m)
    slicer.mrmlScene.RemoveNode(segTap)
    return mascaras


# ============================================================
# PARTE 6 - UTILIDADES DE SEGMENTACION / MALLA
# ============================================================
def _volumen_de_referencia(segNode):
    rol = slicer.vtkMRMLSegmentationNode.GetReferenceImageGeometryReferenceRole()
    volNode = segNode.GetNodeReference(rol)
    if volNode is None:
        vols = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
        volNode = vols[-1] if vols else None
    return volNode


def _ids_de_segmentos(segNode):
    seg = segNode.GetSegmentation()
    return [seg.GetNthSegmentID(i) for i in range(seg.GetNumberOfSegments())]


def _rellenar_cavidades(arr, voxel_mm3, max_mm3):
    fondo = ~arr
    etiquetas, n = ndimage.label(fondo)
    if n == 0:
        return arr, 0
    ids_borde = set()
    for cara in (etiquetas[0, :, :], etiquetas[-1, :, :],
                 etiquetas[:, 0, :], etiquetas[:, -1, :],
                 etiquetas[:, :, 0], etiquetas[:, :, -1]):
        ids_borde.update(np.unique(cara).tolist())
    ids_borde.discard(0)
    conteos = np.bincount(etiquetas.ravel(), minlength=n + 1)
    a_rellenar = [i for i in range(1, n + 1)
                  if i not in ids_borde and conteos[i] * voxel_mm3 <= max_mm3]
    if a_rellenar:
        arr[np.isin(etiquetas, a_rellenar)] = True
    return arr, len(a_rellenar)


def _buscar_o_crear_carpeta(nombre):
    shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    sceneItemId = shNode.GetSceneItemID()
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(sceneItemId, hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if shNode.GetItemName(itemId) == nombre:
            return shNode, itemId
    return shNode, shNode.CreateFolderItem(sceneItemId, nombre)


def _subcarpeta(shNode, parentItemId, nombre):
    hijos = vtk.vtkIdList()
    shNode.GetItemChildren(parentItemId, hijos)
    for i in range(hijos.GetNumberOfIds()):
        itemId = hijos.GetId(i)
        if (shNode.GetItemName(itemId) == nombre
                and shNode.GetItemDataNode(itemId) is None):
            return itemId
    return shNode.CreateFolderItem(parentItemId, nombre)


def _limpiar_malla(modelNode, umbral=UMBRAL_LIMPIEZA_MALLA):
    poly = modelNode.GetPolyData()
    if poly is None or poly.GetNumberOfPoints() == 0:
        return
    conn = vtk.vtkPolyDataConnectivityFilter()
    conn.SetInputData(poly)
    conn.SetExtractionModeToAllRegions()
    conn.ColorRegionsOn()
    conn.Update()
    nRegiones = conn.GetNumberOfExtractedRegions()
    if nRegiones <= 1:
        return
    arrayRegiones = conn.GetOutput().GetPointData().GetArray("RegionId")
    conteo = {}
    for i in range(arrayRegiones.GetNumberOfTuples()):
        rid = int(arrayRegiones.GetTuple1(i))
        conteo[rid] = conteo.get(rid, 0) + 1
    mayor = max(conteo.values())
    conservar = [rid for rid, c in conteo.items() if c >= mayor * umbral]
    conn.SetExtractionModeToSpecifiedRegions()
    conn.InitializeSpecifiedRegionList()
    for rid in conservar:
        conn.AddSpecifiedRegion(rid)
    conn.Update()
    limpiador = vtk.vtkCleanPolyData()
    limpiador.SetInputConnection(conn.GetOutputPort())
    limpiador.Update()
    modelNode.SetAndObservePolyData(limpiador.GetOutput())
    print(f"    limpieza de malla: {nRegiones} regiones -> {len(conservar)}")


def _control_calidad(modelNode):
    poly = modelNode.GetPolyData()
    fe = vtk.vtkFeatureEdges()
    fe.SetInputData(poly)
    fe.BoundaryEdgesOn()
    fe.NonManifoldEdgesOn()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.Update()
    cerrada = (fe.GetOutput().GetNumberOfCells() == 0)
    volumen = 0.0
    if cerrada:
        tri = vtk.vtkTriangleFilter()
        tri.SetInputData(poly)
        tri.Update()
        masa = vtk.vtkMassProperties()
        masa.SetInputConnection(tri.GetOutputPort())
        masa.Update()
        volumen = masa.GetVolume()
    return cerrada, volumen


# ============================================================
# PARTE 7 - ENTRADAS, HIGIENE Y REPORTE
# ============================================================
def _curvas_de_la_escena():
    cerradas = list(slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode'))
    abiertas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')]
    if PREFIJO_CURVAS:
        cerradas = [c for c in cerradas if c.GetName().startswith(PREFIJO_CURVAS)]
        abiertas = [c for c in abiertas if c.GetName().startswith(PREFIJO_CURVAS)]
    return cerradas, abiertas


def _reportar_curvas(cerradas, abiertas):
    print("")
    print("--- CURVAS QUE SE VAN A USAR ---")
    if PREFIJO_CURVAS:
        print(f"  (filtrando por prefijo '{PREFIJO_CURVAS}')")
    for c in cerradas:
        print(f"  [lazo ] {c.GetName()}  ({c.GetNumberOfControlPoints()} pts de control)")
    for c in abiertas:
        print(f"  [linea] {c.GetName()}  ({c.GetNumberOfControlPoints()} pts de control)")
    print(f"  TOTAL: {len(cerradas) + len(abiertas)}")
    print("  >> Si hay MAS curvas de las que dibujaste, son viejas: borralas o")
    print("     usa PREFIJO_CURVAS.")
    print("")


def _preparar_entradas():
    try:
        segOrigen = slicer.util.getNode(NOMBRE_SEGMENTACION)
    except Exception:
        print(f"ERROR: no encuentro el nodo de segmentacion '{NOMBRE_SEGMENTACION}'.")
        return None
    try:
        craneoModel = slicer.util.getNode(NOMBRE_MODELO_CRANEO)
    except Exception:
        print(f"ERROR: no encuentro el modelo '{NOMBRE_MODELO_CRANEO}'. Corre enviar_a_planner().")
        return None
    if not craneoModel.IsA("vtkMRMLModelNode"):
        print(f"ERROR: '{NOMBRE_MODELO_CRANEO}' no es un Model node.")
        return None
    volumeNode = _volumen_de_referencia(segOrigen)
    if volumeNode is None:
        print("ERROR: no encuentro el volumen de referencia del CT.")
        return None
    segIdHueso = segOrigen.GetSegmentation().GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_HUESO)
    if not segIdHueso:
        print(f"ERROR: no hay segmento '{NOMBRE_SEGMENTO_HUESO}'. Corre confirmar_craneo().")
        return None
    cerradas, abiertas = _curvas_de_la_escena()
    if not cerradas and not abiertas:
        print("ERROR: no hay ninguna curva de corte en la escena.")
        return None
    return segOrigen, craneoModel, volumeNode, segIdHueso, cerradas, abiertas


def _muestreo_efectivo(volumeNode):
    if DIST_MUESTREO_MM > 0:
        return DIST_MUESTREO_MM
    return max(min(volumeNode.GetSpacing()) * FACTOR_MUESTREO, 0.05)


def _marcar_puntos(puntos, nombre, color, maximo=300):
    previo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)
    if puntos is None or len(puntos) == 0:
        return 0
    pts = np.asarray(puntos, dtype=float)
    if pts.shape[0] > maximo:
        pts = pts[np.linspace(0, pts.shape[0] - 1, maximo).astype(int)]
    fid = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode', nombre)
    fid.CreateDefaultDisplayNodes()
    for p in pts:
        try:
            fid.AddControlPoint(float(p[0]), float(p[1]), float(p[2]))
        except Exception:
            fid.AddFiducial(float(p[0]), float(p[1]), float(p[2]))
    d = fid.GetDisplayNode()
    if d:
        d.SetSelectedColor(*color)
        d.SetColor(*color)
        d.SetGlyphScale(1.5)
        d.SetTextScale(0.0)
    return int(pts.shape[0])


def _reportar_curva(nombre, tipo, diag):
    n = diag.get("nPuntos", 0)
    nNo = diag.get("nNoAtraviesa", 0)
    frac = diag.get("fracNoAtraviesa", 0.0)
    print(f"  '{nombre}' ({tipo}): {n} puntos")
    if diag.get("busqueda"):
        print(f"      camino mas corto: inclinacion mediana "
              f"{diag.get('inclinacionMediana', 0):.1f} deg, "
              f"max {diag.get('inclinacionMax', 0):.1f} deg")
        print(f"      giro maximo entre puntos vecinos: "
              f"{diag.get('giroMaxPorPaso', 0):.1f} deg "
              f"(limite {MAX_GIRO_POR_PASO_DEG}) -> sin abanicos")
    print(f"      ATRAVIESA EL HUESO: {n - nNo}/{n} puntos "
          f"({100.0*(1.0-frac):.1f}%)")
    # [C3] EL numero que decide
    tmm = diag.get("tiradaMM", 0.0)
    if tmm > 0:
        print(f"      >> PUENTE MAS LARGO SIN CORTAR: {tmm:.1f}mm "
              f"({diag.get('tiradaPuntos', 0)} puntos seguidos)")
        print("         ESTE es el numero que decide, no el porcentaje: la")
        print("         separacion es una cadena, no un promedio. Un solo puente")
        print("         continuo sostiene toda la pieza.")
    else:
        print("      >> sin puentes: el corte atraviesa en TODO el recorrido")
    print(f"      espesor atravesado: mediana {diag.get('espesorMediano', 0):.2f}mm, "
          f"max {diag.get('espesorMax', 0):.2f}mm")
    if nNo > 0:
        print(f"      {nNo} punto(s) no atraviesan; ahi el corte se detiene a "
              f"{diag.get('topeSinSalida', 0):.1f}mm y NO improvisa (marcados en rojo).")
    for lado, cl, dl in (("inicio", "largoExtIni", "distUnionIni"),
                         ("fin", "largoExtFin", "distUnionFin")):
        if cl in diag:
            d = diag[dl]
            if np.isfinite(d) and d <= UNION_TOLERANCIA_MM:
                print(f"      punta {lado}: UNION a {d:.1f}mm de otra curva -> "
                      f"extendida {diag[cl]:.1f}mm para CERRARLA")
            elif np.isfinite(d) and d < 20.0 and d > diag[cl]:
                print(f"      >> punta {lado}: hay otra curva a {d:.1f}mm pero se")
                print(f"         extendio solo {diag[cl]:.1f}mm -> QUEDA UN HUECO de "
                      f"{d - diag[cl]:.1f}mm.")
                print(f"         Si esas dos curvas tenian que unirse, acerca las puntas")
                print(f"         o subi UNION_TOLERANCIA_MM por encima de {d:.1f}.")
            else:
                print(f"      punta {lado}: libre -> extendida {diag[cl]:.1f}mm")
    if diag.get("espesorMax", 0) > ESPESOR_MAX_ESPERADO:
        print(f"      recorrido max > {ESPESOR_MAX_ESPERADO}mm: revisa si el hueso")
        print("         es realmente tan grueso ahi.")
    print(f"      continuidad: separacion max entre rayos "
          f"{diag.get('separacionMaxima', 0):.2f}mm -> "
          f"x{diag.get('factorSupersampleo', 1)} rayos")
    if diag.get("separacionTope"):
        print("      >> el supersampleo llego al tope: subi MAX_FACTOR_SUPERSAMPLEO.")
    if diag.get("sospechosos", 0) > 0:
        print(f"      {diag['sospechosos']} salto(s) de direccion frenados por la guarda")
    if "separo" in diag:
        if diag["separo"]:
            print(f"      SEPARACION: OK (+{diag.get('islasNuevas', 0)} pieza(s))")
        else:
            print("      >> SEPARACION: NO. Este corte quedo grabado pero NO partio")
            print("         el hueso. Mira el PUENTE MAS LARGO de arriba.")
            print("         Recorda: en un craneo cerrado una linea abierta sola NUNCA")
            print("         separa; solo separa un LAZO. Dos lineas separan unicamente")
            print("         si cierran el lazo entre las dos, en LOS DOS extremos.")


# ============================================================
# COMANDO DE DIAGNOSTICO
# ============================================================
def diagnosticar():
    """
    [V1] No corta nada. Dibuja un PICO por cada punto de la curva, con la
    LONGITUD REAL que penetra el corte en el hueso, y coloreado asi:

        AZUL  = ese rayo ATRAVIESA el hueso (el corte pasa de lado a lado)
        ROJO  = ese rayo NO atraviesa (el corte se queda adentro)

    O sea: mirando los picos ves literalmente hasta donde llega el corte y
    donde se queda corto, antes de cortar nada.
    """
    ent = _preparar_entradas()
    if ent is None:
        return
    segOrigen, craneoModel, volumeNode, segIdHueso, cerradas, abiertas = ent
    _reportar_curvas(cerradas, abiertas)

    curvas = list(cerradas) + list(abiertas)
    distMuestreo = _muestreo_efectivo(volumeNode)

    boneMask = slicer.util.arrayFromSegmentBinaryLabelmap(
        segOrigen, segIdHueso, volumeNode).astype(bool)
    print("Preparando el volumen (mascara de espacio abierto)...")
    vol = VolumenHueso(boneMask, volumeNode)

    print(f"Volumen: {volumeNode.GetName()}  spacing={vol.spacing}")
    print(f"Salida del hueso = espacio abierto > {UMBRAL_ESPACIO_ABIERTO_MM}mm "
          f"Y {MARGEN_SALIDA_MM}mm de margen libre")
    print("")

    campo = CampoNormales(craneoModel.GetPolyData())

    puntosVis = vtk.vtkPoints()
    lineas = vtk.vtkCellArray()
    escalar = vtk.vtkFloatArray()
    escalar.SetName("noAtraviesa")
    sinPasar = []

    print("--- DIAGNOSTICO POR CURVA ---")
    for curva in curvas:
        esCerrada = curva.IsA('vtkMRMLMarkupsClosedCurveNode')
        otras = [c for c in curvas if c is not curva]
        pos, dirs, diag = _preparar_curva(curva, campo, esCerrada, otras,
                                          distMuestreo, vol)
        if pos is None:
            print(f"  '{curva.GetName()}': {diag.get('error', 'sin resultado')}")
            continue

        _, rec, atr = _barrido(vol, pos, dirs, None)
        pen, pasa, signo = _resumen_penetracion(rec, atr)
        n = pos.shape[0]
        nNo = int(np.count_nonzero(~pasa))
        tirada = _tirada_mas_larga(pasa, esCerrada)
        tipo = "lazo" if esCerrada else "linea"
        print(f"  '{curva.GetName()}' ({tipo}): {n} puntos")
        if diag.get("busqueda"):
            print(f"      camino mas corto: inclinacion mediana "
                  f"{diag.get('inclinacionMediana', 0):.1f} deg, "
                  f"max {diag.get('inclinacionMax', 0):.1f} deg   |   "
                  f"giro max entre vecinos {diag.get('giroMaxPorPaso', 0):.1f} deg")
        print(f"      ATRAVIESA: {n - nNo}/{n} ({100.0*(n-nNo)/max(n,1):.1f}%)")
        if tirada:
            print(f"      >> PUENTE MAS LARGO SIN CORTAR: {tirada*distMuestreo:.1f}mm "
                  f"({tirada} puntos seguidos)  <-- ESTE decide la separacion")
        else:
            print("      >> sin puentes: atraviesa en TODO el recorrido")
        if np.any(pasa):
            print(f"      espesor atravesado: mediana {np.median(pen[pasa]):.2f}mm, "
                  f"max {np.max(pen[pasa]):.2f}mm")
        for lado, cl, dl in (("inicio", "largoExtIni", "distUnionIni"),
                             ("fin", "largoExtFin", "distUnionFin")):
            if cl in diag:
                dd = diag[dl]
                if np.isfinite(dd) and dd <= UNION_TOLERANCIA_MM:
                    print(f"      punta {lado}: UNION a {dd:.1f}mm de otra curva -> "
                          f"extendida {diag[cl]:.1f}mm para CERRARLA")
                elif np.isfinite(dd) and dd < 20.0 and dd > diag[cl]:
                    print(f"      >> punta {lado}: otra curva a {dd:.1f}mm, extendida "
                          f"solo {diag[cl]:.1f}mm -> HUECO de {dd - diag[cl]:.1f}mm")
                    print(f"         Acerca las puntas o subi UNION_TOLERANCIA_MM.")
                else:
                    print(f"      punta {lado}: libre -> extendida {diag[cl]:.1f}mm")
        if nNo:
            print(f"      {nNo} punto(s) no atraviesan (picos ROJOS)")
            sinPasar.append(pos[~pasa])

        # los picos se dibujan RECORTADOS al mismo tope que usa cortar(): un
        # rayo que no atraviesa recorreria los 25mm enteros y se veria como una
        # varilla larguisima disparada al aire, que asusta y no informa.
        espTip = float(np.median(pen[pasa])) if np.any(pasa) else TOPE_SIN_SALIDA_MIN_MM
        topeDib = max(TOPE_SIN_SALIDA_MIN_MM, TOPE_SIN_SALIDA_FACTOR * espTip)
        paso = max(1, n // 500)
        for i in range(0, n, paso):
            largo = max(min(float(pen[i]), topeDib), 0.6)
            a = pos[i]
            b = pos[i] + dirs[i] * signo[i] * largo
            ia = puntosVis.InsertNextPoint(*a)
            ib = puntosVis.InsertNextPoint(*b)
            ln = vtk.vtkLine()
            ln.GetPointIds().SetId(0, ia)
            ln.GetPointIds().SetId(1, ib)
            lineas.InsertNextCell(ln)
            v = 0.0 if pasa[i] else 1.0
            escalar.InsertNextValue(v)
            escalar.InsertNextValue(v)

    poly = vtk.vtkPolyData()
    poly.SetPoints(puntosVis)
    poly.SetLines(lineas)
    poly.GetPointData().SetScalars(escalar)

    previo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MODELO_PICOS)
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)
    vis = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", NOMBRE_MODELO_PICOS)
    vis.SetAndObservePolyData(poly)
    vis.CreateDefaultDisplayNodes()
    d = vis.GetDisplayNode()
    if d:
        d.SetVisibility(True)
        d.SetLineWidth(3)
        d.SetActiveScalarName("noAtraviesa")
        d.SetScalarVisibility(True)
        d.SetScalarRange(0.0, 1.0)
        d.SetAndObserveColorNodeID("vtkMRMLColorTableNodeFileColdToHotRainbow.txt")

    nNo = _marcar_puntos(np.vstack(sinPasar) if sinPasar else None,
                         NOMBRE_FID_SIN_PASAR, (1.0, 0.0, 0.0))

    print("")
    print(f"Modelo '{NOMBRE_MODELO_PICOS}': cada pico es el corte en ese punto,")
    print("dibujado con la PROFUNDIDAD REAL que alcanza dentro del hueso.")
    print("  AZUL = atraviesa de lado a lado.   ROJO = se queda adentro.")
    if nNo:
        print(f"{nNo} punto(s) rojos tambien marcados en '{NOMBRE_FID_SIN_PASAR}'")
        print("para que puedas hacer click y volar hasta ahi.")
    print("No se corto nada. Cuando estes conforme, corre cortar().")


# ============================================================
# COMANDO PRINCIPAL
# ============================================================
def cortar():
    ent = _preparar_entradas()
    if ent is None:
        return
    segOrigen, craneoModel, volumeNode, segIdHueso, curvasCerradas, curvasAbiertas = ent
    _reportar_curvas(curvasCerradas, curvasAbiertas)

    curvas = list(curvasCerradas) + list(curvasAbiertas)
    distMuestreo = _muestreo_efectivo(volumeNode)

    previo = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_SEG_TRABAJO)
    if previo is not None:
        slicer.mrmlScene.RemoveNode(previo)
    segTrabajo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode',
                                                    NOMBRE_SEG_TRABAJO)
    segTrabajo.CreateDefaultDisplayNodes()
    segTrabajo.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
    segTrabajo.GetSegmentation().CopySegmentFromSegmentation(
        segOrigen.GetSegmentation(), segIdHueso)
    idHueso = _ids_de_segmentos(segTrabajo)[0]

    boneMask = slicer.util.arrayFromSegmentBinaryLabelmap(
        segTrabajo, idHueso, volumeNode).astype(bool)

    print("Preparando el volumen (mascara de espacio abierto)...")
    vol = VolumenHueso(boneMask, volumeNode)
    voxel_mm3 = vol.voxel_mm3
    etiquetasPadre, nPadres = ndimage.label(boneMask)

    print(f"Volumen de referencia: {volumeNode.GetName()}  spacing={vol.spacing}")
    print(f"Config v17: camino mas corto "
          f"{'ON' if BUSCAR_CAMINO_MAS_CORTO else 'OFF'}  "
          f"cono +-{ANGULO_BUSQUEDA_MAX}deg  giro max {MAX_GIRO_POR_PASO_DEG}deg/paso")
    print(f"            kerf objetivo {GROSOR_CORTE_MM}mm   "
          f"salida: espacio abierto >{UMBRAL_ESPACIO_ABIERTO_MM}mm "
          f"+ margen {MARGEN_SALIDA_MM}mm")
    print(f"            sin salida -> corte en seco a "
          f"{TOPE_SIN_SALIDA_FACTOR}x el espesor tipico "
          f"(min {TOPE_SIN_SALIDA_MIN_MM}mm)")
    print(f"Islas de hueso antes del corte: {nPadres}")

    campo = CampoNormales(craneoModel.GetPolyData())
    taponesMasc = _mascaras_de_tapones(volumeNode, 2.0 * PROFUNDIDAD_MAX_MM,
                                       curvasCerradas)
    print(f"Tapones generados: {list(taponesMasc.keys()) if taponesMasc else 'ninguno'}")

    print("")
    print("--- CONSTRUCCION DE LOS CORTES ---")
    nombresAbiertas = set(c.GetName() for c in curvasAbiertas)
    anilloTotal = np.zeros_like(boneMask)
    sinPasar, noSepararon = [], []
    kerfReal = None
    for curva in curvas:
        esCerrada = curva.IsA('vtkMRMLMarkupsClosedCurveNode')
        otras = [c for c in curvas if c is not curva]
        anillo, diag = _anillo_de_corte(vol, curva, campo, esCerrada, otras,
                                        distMuestreo, nPadres)
        if anillo is None:
            print(f"  ADVERTENCIA: no se pudo construir el corte de "
                  f"'{curva.GetName()}': {diag.get('error', 'motivo desconocido')}")
            continue
        anilloTotal |= anillo
        kerfReal = diag.get("kerfReal", kerfReal)
        _reportar_curva(curva.GetName(), "lazo" if esCerrada else "linea", diag)
        pr = diag.get("puntosSinPasarXYZ")
        if pr is not None and len(pr):
            sinPasar.append(pr)
        if diag.get("separo") is False:
            noSepararon.append(curva.GetName())

    if not anilloTotal.any():
        print("ERROR: no se pudo construir ningun corte.")
        slicer.mrmlScene.RemoveNode(segTrabajo)
        return

    if kerfReal:
        print("")
        print(f"Ancho REAL de la ranura: {kerfReal:.2f}mm "
              f"({int(round(kerfReal/min(vol.spacing)))} voxeles)")

    n = _marcar_puntos(np.vstack(sinPasar) if sinPasar else None,
                       NOMBRE_FID_SIN_PASAR, (1.0, 0.0, 0.0))
    if n:
        print(f"{n} punto(s) que NO atraviesan, marcados en rojo en "
              f"'{NOMBRE_FID_SIN_PASAR}'.")
    if noSepararon:
        print(f"AVISO: {len(noSepararon)} corte(s) NO separaron el hueso: {noSepararon}")

    boneCut = np.logical_and(boneMask, np.logical_not(anilloTotal))
    print("")
    print(f"Corte total restado: {int(np.count_nonzero(anilloTotal))} voxeles "
          f"({np.count_nonzero(anilloTotal)*voxel_mm3/1000.0:.3f}cm3).")
    etiquetasHijas, nHijas = ndimage.label(boneCut)
    print(f"Islas tras el corte: {nHijas}")

    def curva_interior_de(mask, nvox):
        mejor, mejorFrac = None, 0.0
        for nombreCurva, tap in taponesMasc.items():
            if tap.shape != mask.shape:
                continue
            frac = np.count_nonzero(mask & tap) / float(nvox) if nvox else 0.0
            if frac > mejorFrac:
                mejor, mejorFrac = nombreCurva, frac
        return mejor if mejorFrac >= FRACCION_INTERIOR else None

    infoLab, hijasDeMadre = {}, {}
    for lab in range(1, nHijas + 1):
        m = (etiquetasHijas == lab)
        nvox = int(np.count_nonzero(m))
        if nvox < RUIDO_VOXELES:
            continue
        et = etiquetasPadre[m]
        et = et[et > 0]
        madre = int(np.bincount(et).argmax()) if et.size else 0
        infoLab[lab] = (nvox, madre, curva_interior_de(m, nvox))
        hijasDeMadre.setdefault(madre, []).append(lab)

    madreTieneTapa = {madre: any(infoLab[l][2] is not None for l in labs)
                      for madre, labs in hijasDeMadre.items()}

    piezas, contadorTapa = [], {}
    contadorResto = contadorFragmento = contadorHueso = 0
    for madre, labs in hijasDeMadre.items():
        creadaPorCorte = len(labs) > 1
        for lab in labs:
            nvox, _, curva = infoLab[lab]
            mask = (etiquetasHijas == lab)
            vol_cm3 = nvox * voxel_mm3 / 1000.0
            if curva is not None:
                contadorTapa[curva] = contadorTapa.get(curva, 0) + 1
                piezas.append((mask, f"Tapa_{curva}", True, curva, vol_cm3))
            elif creadaPorCorte and madreTieneTapa[madre]:
                contadorResto += 1
                piezas.append((mask, f"Resto_{contadorResto}", False, None, vol_cm3))
            elif creadaPorCorte:
                contadorFragmento += 1
                piezas.append((mask, f"Fragmento_{contadorFragmento}", False, None, vol_cm3))
            else:
                if vol_cm3 < VOLUMEN_MINIMO_CM3:
                    continue
                contadorHueso += 1
                piezas.append((mask, f"Hueso_{contadorHueso}", False, None, vol_cm3))

    if not piezas:
        print("ERROR: no sobrevivio ningun fragmento.")
        slicer.mrmlScene.RemoveNode(segTrabajo)
        return

    sufijos = {c: 0 for c, k in contadorTapa.items() if k > 1}
    piezasFinales = []
    for mask, nombre, esTapa, curva, v in piezas:
        if esTapa and curva in sufijos:
            nombre = f"{nombre}_{chr(ord('a') + sufijos[curva])}"
            sufijos[curva] += 1
        piezasFinales.append((mask, nombre, esTapa, curva, v))

    print("")
    print("--- PIEZAS ---")
    idsPorNombre, curvaPorNombre = {}, {}
    for mask, nombre, esTapa, curva, v in piezasFinales:
        arr, nRell = _rellenar_cavidades(mask.copy(), voxel_mm3, MAX_CAVIDAD_MM3)
        sid = segTrabajo.GetSegmentation().AddEmptySegment(nombre)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            arr.astype(np.uint8), segTrabajo, sid, volumeNode)
        idsPorNombre[nombre] = sid
        curvaPorNombre[nombre] = curva if esTapa else None
        volf = np.count_nonzero(arr) * voxel_mm3 / 1000.0
        print(f"  {nombre} [{'TAPA' if esTapa else 'resto/hueso'}]: "
              f"vol={volf:.4f}cm3 ({nRell} cavidad(es) rellenada(s))")

    segTrabajo.GetSegmentation().RemoveSegment(idHueso)
    nombresTapas = set(n for _, n, e, _, _ in piezasFinales if e)

    segTrabajo.GetSegmentation().SetConversionParameter("Smoothing factor", SUAVIZADO_EXPORT)
    segTrabajo.RemoveClosedSurfaceRepresentation()
    segTrabajo.CreateClosedSurfaceRepresentation()

    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)
    if not slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToModels(
            segTrabajo, list(idsPorNombre.values()), folderItemId):
        print("Fallo la exportacion a modelos. Revisa el Error log de Slicer.")
        return

    print("")
    print("--- CONTROL DE CALIDAD ---")
    nodoPorNombre = {}
    modelos = slicer.util.getNodesByClass("vtkMRMLModelNode")
    listos = total = 0
    for nombre in idsPorNombre.keys():
        cand = [m for m in modelos if m.GetName() == nombre
                and shNode.GetItemParent(shNode.GetItemByDataNode(m)) == folderItemId]
        if not cand:
            cand = [m for m in modelos if m.GetName() == nombre]
        if not cand:
            print(f"  ADVERTENCIA: no encontre el modelo '{nombre}'.")
            continue
        node = cand[0]
        nodoPorNombre[nombre] = node
        total += 1
        _limpiar_malla(node)
        cerrada, volumen = _control_calidad(node)
        print(f"  {nombre} [{'tapa de hueso' if nombre in nombresTapas else 'resto/hueso'}]: "
              f"{'WATERTIGHT (apto STL)' if cerrada else 'ABIERTA (NO apta para STL)'}  "
              f"vol={volumen/1000.0:.4f}cm3")
        if nombre in nombresTapas and node.GetDisplayNode():
            node.GetDisplayNode().SetColor(0.9, 0.25, 0.25)
        if cerrada:
            listos += 1

    itemsPorCurva = {}
    for nombre, node in nodoPorNombre.items():
        curva = curvaPorNombre.get(nombre)
        if curva is None:
            continue
        itemId = shNode.GetItemByDataNode(node)
        if itemId:
            itemsPorCurva.setdefault(curva, []).append(itemId)
    for curva, items in itemsPorCurva.items():
        sub = _subcarpeta(shNode, folderItemId, f"Tapas_{curva}")
        for itemId in items:
            shNode.SetItemParent(itemId, sub)

    for nodo in (craneoModel, segOrigen, segTrabajo):
        if nodo and nodo.GetDisplayNode():
            nodo.GetDisplayNode().SetVisibility(False)
    visN = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MODELO_PICOS)
    if visN is not None and visN.GetDisplayNode():
        visN.GetDisplayNode().SetVisibility(False)

    print("")
    print(f"Listo: {total} pieza(s), {listos} watertight, "
          f"{len(nombresTapas)} tapa(s) de hueso (en rojo).")
    print(f"Estan en la carpeta '{NOMBRE_CARPETA_SH}' del panel Data.")
    print("Para exportar STL:  exportar_stl('C:/ruta/de/salida')")


def exportar_stl(carpeta_destino):
    """Exporta a STL todos los fragmentos watertight."""
    import os
    if not os.path.isdir(carpeta_destino):
        os.makedirs(carpeta_destino)
    shNode, folderItemId = _buscar_o_crear_carpeta(NOMBRE_CARPETA_SH)

    def _recorrer(itemId, salida):
        hijos = vtk.vtkIdList()
        shNode.GetItemChildren(itemId, hijos)
        for i in range(hijos.GetNumberOfIds()):
            hid = hijos.GetId(i)
            node = shNode.GetItemDataNode(hid)
            if node is not None and node.IsA("vtkMRMLModelNode"):
                salida.append(node)
            else:
                _recorrer(hid, salida)

    nodos = []
    _recorrer(folderItemId, nodos)
    n = 0
    for node in nodos:
        cerrada, _ = _control_calidad(node)
        if not cerrada:
            print(f"  {node.GetName()}: NO watertight, no se exporta.")
            continue
        ruta = os.path.join(carpeta_destino, node.GetName() + ".stl")
        slicer.util.saveNode(node, ruta)
        print(f"  {node.GetName()} -> {ruta}")
        n += 1
    print(f"{n} STL exportado(s).")


def limpiar():
    """Borra el nodo de trabajo y los nodos de diagnostico."""
    for nombre in (NOMBRE_SEG_TRABAJO, NOMBRE_MODELO_PICOS, NOMBRE_FID_SIN_PASAR):
        nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
        while nodo is not None:
            slicer.mrmlScene.RemoveNode(nodo)
            nodo = slicer.mrmlScene.GetFirstNodeByName(nombre)
    print("Nodos de trabajo borrados.")


print("Bloque F v17 cargado.")
if BUSCAR_CAMINO_MAS_CORTO:
    print(f"  direccion: BUSCA EL CAMINO MAS CORTO (cono +-{ANGULO_BUSQUEDA_MAX} deg),")
    print(f"             pero la elige para la CURVA ENTERA con programacion dinamica")
    print(f"             y un giro maximo de {MAX_GIRO_POR_PASO_DEG} deg entre puntos vecinos:")
    print( "             por construccion no puede haber tajos random.")
else:
    print("  direccion: siempre perpendicular (busqueda desactivada)")
print(f"  verificacion por rayo: atraviesa el hueso? (espacio abierto "
      f">{UMBRAL_ESPACIO_ABIERTO_MM}mm + margen {MARGEN_SALIDA_MM}mm)")
print(f"  sin salida -> corte en seco a {TOPE_SIN_SALIDA_FACTOR}x el espesor "
      f"tipico: no puede hacer dano lejos del sitio")
print(f"  uniones entre curvas: se extienden hasta TOCAR (+{SOLAPE_UNION_MM}mm), "
      f"ya no quedan huecos")
print("  se reporta el PUENTE MAS LARGO sin cortar: ese numero decide, no el %")
print("Comandos:")
print("  diagnosticar()              -> PICOS: azul atraviesa, rojo NO. Sin cortar.")
print("  cortar()                    -> corte + verificacion por curva")
print("  exportar_stl('C:/carpeta')  -> exporta los fragmentos watertight")
print("  limpiar()                   -> borra los nodos de trabajo")
