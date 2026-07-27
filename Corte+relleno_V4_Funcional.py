# ============================================================
# CRANIOPLAN - BLOQUE F (CORTE) - SCRIPT DE CONSOLA
# ============================================================
# Versión del algoritmo: 2026-07-26-j (flap por partición geométrica del prisma del lazo: siempre separa)
#
# Pensado para probar el corte pegando el código en la consola de Python de
# 3D Slicer, sin tener instalado el módulo CranioPlan. Reemplaza al viejo
# "Bloque F v4": mismo formato (CONFIG arriba + comando cortar()), pero con el
# motor de corte nuevo:
#   - anillo de corte construido en vóxeles (barrido por la normal local
#     suavizada), inmune al 'twist' y a huecos por aliasing;
#   - curva CERRADA: el flap = hueso encerrado por el lazo, por PARTICIÓN
#     GEOMÉTRICA del prisma del lazo (no por conectividad) -> siempre separa,
#     aunque el lazo pase cerca de un agujero/sutura abierta;
#   - curva ABIERTA: identidad por piezas conexas.
#
# FLUJO (igual que veníamos):
#   1) Cargar la tomografía del paciente.
#   2) Correr el Bloque A+B por consola -> eliminar islas -> confirmar_craneo()
#      -> enviar_a_planner()  (deja el Model "Craneo_Final" en la escena).
#   3) Ir a Markups y dibujar una o varias curvas sobre "Craneo_Final"
#      (cerradas y/o abiertas; "Shortest distance on surface" + "Constrain to
#      Model" quedan más pegadas a la superficie).
#   4) Pegar este código en la consola y correr:  cortar()
#
# COMANDOS:
#   cortar()   -> corta el cráneo con TODAS las curvas de la escena, en orden.
#                 Crea "Craneo_restante" y "Fragmento_extraido_N" por cada flap.
#
# REQUISITO: scipy. Si no está, correr una vez en la consola de Slicer:
#   slicer.util.pip_install("scipy")
# ============================================================

import numpy as np
import vtk
import slicer


# ==================== CONFIGURACIÓN ====================
GROSOR_CORTE_MM        = 1.0    # ancho de hoja (kerf). Piso práctico ~2 vóxeles.
MARGEN_SEGURIDAD_MM    = 3.0    # margen extra de profundidad del corte
REDUCCION_MALLA        = 0.7    # decimación de las mallas de salida (0 = sin decimar)
VOL_MIN_FRAGMENTO_MM3  = 300.0  # solo curva ABIERTA: umbral de fragmento vs esquirla

NOMBRE_MODELO_CRANEO   = "Craneo_Final"       # Model node que deja enviar_a_planner()
NOMBRE_SEGMENTACION    = "Craneo_Automatico"  # fallback: nodo de segmentación del Bloque A
NOMBRE_SEGMENTO_HUESO  = "Craneo_Final"       # fallback: segmento dentro de ese nodo

CRANIOPLAN_VERSION = "2026-07-26-j (flap por partición geométrica del prisma del lazo: siempre separa)"


class LogicaCorteOsteotomia:
    """Pipeline de corte de osteotomía (Bloque F de CranioPlan)."""

    def contarPiezasConectadas(self, polyData):
        """Cantidad de componentes conectados de una malla."""
        if polyData is None or polyData.GetNumberOfPoints() == 0:
            return 0
        conectividad = vtk.vtkPolyDataConnectivityFilter()
        conectividad.SetInputData(polyData)
        conectividad.SetExtractionModeToAllRegions()
        conectividad.Update()
        return int(conectividad.GetNumberOfExtractedRegions())

    def _decimarMalla(self, polyData, reduccion):
        """
        Reduce la cantidad de triángulos de una malla, preservando la
        topología (no abre agujeros ni separa piezas).

        reduccion: fracción de triángulos a eliminar (0.0 a 1.0).
        0.0 desactiva la decimación.
        """
        if not reduccion or reduccion <= 0.0 or polyData is None:
            return polyData
        if polyData.GetNumberOfCells() == 0:
            return polyData

        decimador = vtk.vtkDecimatePro()
        decimador.SetInputData(polyData)
        decimador.SetTargetReduction(reduccion)
        decimador.PreserveTopologyOn()
        decimador.Update()

        resultado = vtk.vtkPolyData()
        resultado.DeepCopy(decimador.GetOutput())

        if resultado.GetNumberOfPoints() == 0:
            return polyData  # la decimación falló; devolvemos la original
        return resultado

    def crearModeloDesdePolyData(self, polyData, nombre, color):
        """
        Crea un vtkMRMLModelNode a partir de un vtkPolyData, con nombre
        único (GenerateUniqueName) y color/visualización estándar.
        """
        if polyData is None or polyData.GetNumberOfPoints() == 0:
            return None
        nombreUnico = slicer.mrmlScene.GenerateUniqueName(nombre)
        nodo = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', nombreUnico)
        nodo.SetAndObservePolyData(polyData)
        nodo.CreateDefaultDisplayNodes()
        nodo.GetDisplayNode().SetColor(*color)
        nodo.GetDisplayNode().SetScalarVisibility(False)
        return nodo

    def _resamplearPuntos(self, puntosOriginales, distanciaMuestreoMM):
        """
        Recorre la polilínea de una curva y devuelve puntos (numpy arrays)
        espaciados uniformemente cada distanciaMuestreoMM. Implementación
        propia, sin depender de utilidades internas de Slicer.
        """
        n = puntosOriginales.GetNumberOfPoints()
        if n < 2:
            return []

        original = [np.array(puntosOriginales.GetPoint(i)) for i in range(n)]

        resampleados = [original[0]]
        distanciaAcumulada = 0.0
        puntoAnterior = original[0]

        for i in range(1, n):
            puntoActual = original[i]
            segmento = puntoActual - puntoAnterior
            largoSegmento = np.linalg.norm(segmento)

            while distanciaAcumulada + largoSegmento >= distanciaMuestreoMM:
                falta = distanciaMuestreoMM - distanciaAcumulada
                direccion = segmento / largoSegmento if largoSegmento > 1e-9 else segmento
                nuevoPunto = puntoAnterior + direccion * falta
                resampleados.append(nuevoPunto)
                puntoAnterior = nuevoPunto
                segmento = puntoActual - puntoAnterior
                largoSegmento = np.linalg.norm(segmento)
                distanciaAcumulada = 0.0

            distanciaAcumulada += largoSegmento
            puntoAnterior = puntoActual

        if np.linalg.norm(resampleados[-1] - original[-1]) > 1e-6:
            resampleados.append(original[-1])

        return resampleados

    def _normalPromedioEnPunto(self, punto, pointLocator, normalesArray, radioMM=3.0):
        """
        Normal de superficie promediada sobre los vértices dentro de un
        radio, en vez de tomar la de una sola celda (más estable sobre
        mallas de marching cubes).
        """
        listaIds = vtk.vtkIdList()
        pointLocator.FindPointsWithinRadius(radioMM, punto.tolist(), listaIds)

        if listaIds.GetNumberOfIds() == 0:
            idCercano = pointLocator.FindClosestPoint(punto.tolist())
            if idCercano < 0:
                return None
            listaIds.InsertNextId(idCercano)

        acumulado = np.zeros(3)
        for j in range(listaIds.GetNumberOfIds()):
            acumulado += np.array(normalesArray.GetTuple3(listaIds.GetId(j)))

        norma = np.linalg.norm(acumulado)
        if norma < 1e-6:
            idCercano = pointLocator.FindClosestPoint(punto.tolist())
            if idCercano < 0:
                return None
            normal = np.array(normalesArray.GetTuple3(idCercano))
            norma = np.linalg.norm(normal)
            return normal / norma if norma > 1e-9 else None

        return acumulado / norma

    def _quitarPuntosCoincidentes(self, posiciones, esCerrada, toleranciaMM=0.2):
        """
        Elimina puntos consecutivos que están (casi) en la misma posición.
        Un tramo de longitud cero rompe el cálculo de la tangente (NaN) y
        arruina la pared de corte sin disparar ningún error.
        """
        limpias = []
        for p in posiciones:
            if not limpias or np.linalg.norm(p - limpias[-1]) > toleranciaMM:
                limpias.append(p)

        if esCerrada and len(limpias) > 2:
            if np.linalg.norm(limpias[-1] - limpias[0]) <= toleranciaMM:
                limpias.pop()

        return limpias

    def _anilloDeCorteEnVoxeles(self, boneMask, curvaNode, mallaHueso, volumeNode,
                                 grosorMM, margenSeguridadMM, esCerrada,
                                 distanciaMuestreoMM=0.5, radioNormalMM=3.0):
        """
        Construye el anillo (curva cerrada) o canal (curva abierta) de corte
        DIRECTAMENTE en el espacio de vóxeles, en vez de mallar un tubo fino
        y esperar que la vóxelización lo llene. Por cada punto de la curva se
        barre a lo largo de la normal del hueso (para atravesar todo el
        espesor) y luego se dilata el resultado el grosor de hoja. Así el
        anillo es continuo por construcción: inmune al 'twist' del lazo y a
        los huecos por aliasing, que eran la causa de que el disco no se
        separara.

        Devuelve una máscara booleana (misma forma que boneMask) con los
        vóxeles de HUESO a quitar, o None si no se pudo construir.
        """
        from scipy import ndimage

        puntosCurva = curvaNode.GetCurvePointsWorld()
        if puntosCurva is None or puntosCurva.GetNumberOfPoints() < 2:
            return None

        posiciones = self._resamplearPuntos(puntosCurva, distanciaMuestreoMM)
        posiciones = self._quitarPuntosCoincidentes(posiciones, esCerrada)
        posiciones = list(posiciones)
        if esCerrada and len(posiciones) >= 3:
            posiciones.append(posiciones[0])  # cerrar el lazo explícitamente
        n = len(posiciones)
        if n < 2:
            return None
        print(
            f"CranioPlan DIAGNÓSTICO: curva {'cerrada' if esCerrada else 'abierta'} "
            f"resampleada a {n} puntos para el barrido en vóxeles."
        )

        # Normales del hueso (para barrer a lo largo del espesor).
        normalesFilter = vtk.vtkPolyDataNormals()
        normalesFilter.SetInputData(mallaHueso)
        normalesFilter.ComputePointNormalsOn()
        normalesFilter.ComputeCellNormalsOff()
        normalesFilter.SplittingOff()
        normalesFilter.ConsistencyOn()
        normalesFilter.AutoOrientNormalsOn()
        normalesFilter.Update()
        mallaN = normalesFilter.GetOutput()
        normalesArray = mallaN.GetPointData().GetNormals()
        if normalesArray is None:
            print("CranioPlan DIAGNÓSTICO: no se pudieron calcular normales del hueso.")
            return None
        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(mallaN)
        pointLocator.BuildLocator()

        rasToIjk = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(rasToIjk)

        dims = boneMask.shape  # (z, y, x)
        espaciado = volumeNode.GetSpacing()
        pasoMM = max(min(espaciado) * 0.5, 0.1)
        PROFUNDIDAD_MM = 8.0  # perpendicular al hueso: no ensancha el corte,
                              #  solo garantiza atravesar el espesor completo

        # --- Dirección de extrusión por punto: normal local SUAVIZADA ---
        # Extruir a lo largo de la normal local (perpendicular al hueso)
        # mantiene el corte fino y pegado a la línea trazada. Pero las normales
        # crudas cerca de huecos/zonas finas tienen picos que producían gubias
        # anchas (sobre todo en curvas abiertas). Se suavizan con un promedio
        # móvil a lo largo de la curva: se conserva la curvatura general y se
        # matan los picos. Reemplaza al 'eje global' anterior, que cortaba de
        # más en las zonas inclinadas.
        normalesPunto = []
        for p in posiciones:
            normalesPunto.append(
                self._normalPromedioEnPunto(p, pointLocator, normalesArray, radioNormalMM)
            )

        indicesValidos = [i for i, nrm in enumerate(normalesPunto) if nrm is not None]
        if not indicesValidos:
            print("CranioPlan DIAGNÓSTICO: no se pudo calcular ninguna normal sobre la curva.")
            return None
        # Rellenar faltantes con la normal válida más cercana.
        for i in range(n):
            if normalesPunto[i] is None:
                j = min(indicesValidos, key=lambda k: abs(k - i))
                normalesPunto[i] = np.array(normalesPunto[j])

        # Alinear todas al sentido medio (evita promediar normales opuestas).
        media = np.zeros(3)
        for i in indicesValidos:
            media = media + normalesPunto[i]
        if np.linalg.norm(media) > 1e-9:
            media = media / np.linalg.norm(media)
            for i in range(n):
                if float(np.dot(normalesPunto[i], media)) < 0:
                    normalesPunto[i] = -normalesPunto[i]

        # Promedio móvil a lo largo de la curva (circular si es cerrada).
        ventana = 3  # a cada lado => 7 puntos ~ 3 mm
        normalesSuaves = []
        for i in range(n):
            acum = np.zeros(3)
            for d in range(-ventana, ventana + 1):
                if esCerrada:
                    k = (i + d) % n
                else:
                    k = min(max(i + d, 0), n - 1)
                acum = acum + normalesPunto[k]
            norma = np.linalg.norm(acum)
            normalesSuaves.append(acum / norma if norma > 1e-9 else normalesPunto[i])

        curtain = np.zeros(dims, dtype=bool)

        def marcar(rasXYZ):
            ijk = [0.0, 0.0, 0.0, 0.0]
            rasToIjk.MultiplyPoint([float(rasXYZ[0]), float(rasXYZ[1]), float(rasXYZ[2]), 1.0], ijk)
            ii = int(round(ijk[0]))
            jj = int(round(ijk[1]))
            kk = int(round(ijk[2]))
            if 0 <= kk < dims[0] and 0 <= jj < dims[1] and 0 <= ii < dims[2]:
                curtain[kk, jj, ii] = True

        for i, p in enumerate(posiciones):
            direccion = normalesSuaves[i]
            t = -PROFUNDIDAD_MM
            while t <= PROFUNDIDAD_MM:
                marcar(np.asarray(p, dtype=float) + direccion * t)
                t += pasoMM

        if not curtain.any():
            return None

        # Dilatar el grosor de hoja (kerf). Estructura 3x3x3 (26-conexa): es
        # la que garantiza un anillo CONTINUO sin huecos diagonales sobre la
        # trayectoria curva. La 6-conexa dejaba fugas de 1 vóxel y el disco no
        # se separaba. El ancho lateral queda ~1.5 mm; el corte ya es preciso
        # porque la extrusión es perpendicular al hueso (sin inclinación).
        radioVox = max(1, int(round((grosorMM / 2.0) / min(espaciado))))
        estructura = np.ones((3, 3, 3), dtype=bool)
        curtainDil = ndimage.binary_dilation(curtain, structure=estructura, iterations=radioVox)

        anillo = np.logical_and(curtainDil, boneMask)
        if not anillo.any():
            return None
        return anillo

    def _puntosDentroPoligono2D(self, P, poligono):
        """Vectorizado: para P (Nx2) devuelve un bool array indicando si cada
        punto cae dentro del polígono (Mx2), por ray casting."""
        x = P[:, 0]
        y = P[:, 1]
        dentro = np.zeros(len(P), dtype=bool)
        n = len(poligono)
        j = n - 1
        for i in range(n):
            xi, yi = poligono[i]
            xj, yj = poligono[j]
            cond = ((yi > y) != (yj > y)) & \
                   (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-12) + xi)
            dentro = np.logical_xor(dentro, cond)
            j = i
        return dentro

    def _particionPorPrismaDelLazo(self, boneMask, loopRAS, volumeNode, slabMM=25.0):
        """
        Parte el hueso (boneMask) en (flapMask, restanteMask) según el PRISMA
        del lazo cerrado: un vóxel pertenece al flap si su proyección cae
        DENTRO del polígono del lazo (proyectado a su plano de mejor ajuste) Y
        está a menos de ±slabMM del plano. Todo lo demás va al restante.

        Es una partición puramente GEOMÉTRICA: no usa conectividad, así que
        SIEMPRE separa correctamente el hueso encerrado por el lazo, aunque el
        anillo de corte no lo haya desconectado (que es lo que fallaba cerca de
        los agujeros/suturas abiertas, donde quedaban puentes de hueso).

        La franja ±slabMM evita que el prisma (infinito) capture hueso de la
        pared opuesta del cráneo que quede sobre el mismo eje.
        """
        pts = np.asarray(loopRAS, dtype=float)
        if pts.shape[0] < 3:
            return None, None
        c0 = pts.mean(axis=0)
        _, _, vh = np.linalg.svd(pts - c0)
        u = vh[0]
        v = vh[1]
        eje = vh[2]
        poligono = np.column_stack([(pts - c0) @ u, (pts - c0) @ v])

        indices = np.argwhere(boneMask)  # N x 3 (z, y, x)
        if indices.shape[0] == 0:
            return None, None

        M = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(M)
        Mnp = np.array([[M.GetElement(r, c) for c in range(4)] for r in range(4)])
        # IJK homogéneo por vóxel: (i=x, j=y, k=z, 1)
        ijk = np.column_stack([
            indices[:, 2], indices[:, 1], indices[:, 0], np.ones(len(indices))
        ]).astype(float)
        ras = (ijk @ Mnp.T)[:, :3]

        rel = ras - c0
        s = rel @ eje
        a = rel @ u
        b = rel @ v
        enSlab = np.abs(s) <= slabMM
        dentroPoli = self._puntosDentroPoligono2D(np.column_stack([a, b]), poligono)
        esFlap = enSlab & dentroPoli

        flapMask = np.zeros_like(boneMask)
        restMask = np.zeros_like(boneMask)
        zc = indices[:, 0]
        yc = indices[:, 1]
        xc = indices[:, 2]
        flapMask[zc[esFlap], yc[esFlap], xc[esFlap]] = True
        noFlap = ~esFlap
        restMask[zc[noFlap], yc[noFlap], xc[noFlap]] = True
        return flapMask, restMask

    def generarOsteotomia(self, mallaRemanente, curvaNode, volumeNode,
                            indiceInicialFragmento=1,
                            grosorMM=1.0, margenSeguridadMM=3.0,
                            volumenMinimoFragmentoMM3=300.0,
                            reduccionMallaFragmentos=0.7):
        """
        Ejecuta un corte de osteotomía sobre el cráneo remanente.

        mallaRemanente: vtkPolyData con TODO el hueso que todavía no se
        extrajo (el cráneo restante actual). Puede tener varias piezas
        conexas (placas separadas por suturas abiertas). El corte se aplica
        solo sobre este hueso.

        SEGUIMIENTO DE IDENTIDAD — corrección de fondo (26/07/2026):
        --------------------------------------------------------------
        El bug de raíz de todas las versiones anteriores era decidir "qué es
        cráneo restante" y "qué es fragmento extraído" por TAMAÑO: la pieza
        más grande = restante, todas las demás = extraídas. Eso está mal
        porque un cráneo pediátrico ya viene partido en varias placas
        grandes ANTES de cortar (suturas abiertas). El ranking por tamaño
        tomaba esas placas naturales y las llamaba "fragmentos extraídos",
        aunque el corte nunca las tocó.

        La solución correcta es rastrear la IDENTIDAD de cada pieza:
          1. Se voxeliza el hueso de entrada y se etiquetan sus componentes
             conexos ANTES de cortar (scipy.ndimage.label -> "piezas madre").
          2. Se resta la pared de corte (solo QUITA vóxeles) y se re-separan
             las islas -> "piezas hija".
          3. Como restar solo quita vóxeles, cada hija es subconjunto de
             exactamente una madre. Se mapea cada hija a su madre por la
             etiqueta mayoritaria bajo su máscara.
          4. Por cada madre:
               - si produjo <=1 hija real -> el corte NO la separó: va
                 ENTERA al cráneo restante (aunque no sea la más grande).
               - si produjo >=2 hijas reales -> el corte SÍ la separó: la
                 mayor queda en el restante, las otras son flaps extraídos,
                 y las esquirlas (< volumenMinimoFragmentoMM3) vuelven al
                 restante (conserva el volumen óseo).
          5. Craneo_restante = unión (vtkAppendPolyData) de todas las piezas
             del bucket restante, en UN solo modelo (puede tener varias
             placas). Un modelo Fragmento_extraido por flap.

        Así, con un único corte de lazo cerrado sobre una placa, el
        resultado es exactamente: 1 fragmento extraído + el resto del cráneo
        (todas las demás placas + el remanente de la placa cortada), sin
        importar cuántas placas naturales haya ni cuál sea la más grande.

        indiceInicialFragmento: número con el que arranca la numeración de
        los fragmentos de ESTE corte (para que en cortes sucesivos los
        nombres sigan la cuenta: Fragmento_extraido_1, _2, _3, ...).

        volumenMinimoFragmentoMM3: umbral absoluto (300 mm³ = 0.3 cm³ por
        defecto) para distinguir un flap real de una esquirla de la
        vóxelización. Provisional: falta validar con más casos del Garrahan
        que no descarte un fragmento pediátrico legítimamente chico.

        MOTOR DE CORTE: resta volumétrica (por vóxeles) con el Segment
        Editor, no booleano de mallas (vtkBooleanOperationPolyDataFilter
        falla sistemáticamente en mallas craneales reales). Precisión
        limitada al vóxel (0.5 mm), más fino que el grosor de hoja (~1 mm).

        Devuelve un dict:
          {"restante": modelNode, "fragmentos": [modelNode, ...],
           "piezasCreadas": int, "piezasAntes": int}
        o None si el corte no se pudo calcular.
        """
        print(f"CranioPlan {CRANIOPLAN_VERSION}: Bloque F (generarOsteotomia).")

        if mallaRemanente is None or curvaNode is None or volumeNode is None:
            print("CranioPlan DIAGNÓSTICO: falta el hueso, la curva o el volumen.")
            return None
        if mallaRemanente.GetNumberOfPoints() == 0:
            print("CranioPlan DIAGNÓSTICO: la malla de hueso remanente está vacía.")
            return None

        try:
            from scipy import ndimage
        except ImportError:
            print(
                "CranioPlan DIAGNÓSTICO: scipy no está disponible. El seguimiento de "
                "identidad de piezas lo necesita; sin él no se puede distinguir un flap "
                "real de una placa natural. Se aborta el corte."
            )
            return None

        esCerrada = bool(curvaNode.IsA("vtkMRMLMarkupsClosedCurveNode"))
        print(
            f"CranioPlan DIAGNÓSTICO: curva {curvaNode.GetClassName()} -> "
            f"{'CERRADA' if esCerrada else 'ABIERTA'}."
        )

        piezasAntesMalla = self.contarPiezasConectadas(mallaRemanente)
        print(
            f"CranioPlan: hueso remanente de entrada — {piezasAntesMalla} pieza(s) "
            "conexa(s) según la malla."
        )

        # --- Segmentación temporal donde se hace el corte por vóxeles ---
        segCorte = slicer.mrmlScene.AddNewNodeByClass(
            'vtkMRMLSegmentationNode', 'CranioPlan_Corte_Temporal'
        )
        segCorte.CreateDefaultDisplayNodes()
        segCorte.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        idHueso = segCorte.AddSegmentFromClosedSurfaceRepresentation(
            mallaRemanente, "Hueso", [0.9, 0.8, 0.6]
        )
        if not idHueso:
            print("CranioPlan DIAGNÓSTICO: no se pudo importar el hueso a la segmentación.")
            slicer.mrmlScene.RemoveNode(segCorte)
            return None

        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segEditorNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentEditorNode')
        segmentEditorWidget.setMRMLSegmentEditorNode(segEditorNode)
        segmentEditorWidget.setSegmentationNode(segCorte)
        segmentEditorWidget.setSourceVolumeNode(volumeNode)

        # === PASO 1 de identidad: etiquetar las piezas madre ANTES de cortar ===
        boneArrAntes = slicer.util.arrayFromSegmentBinaryLabelmap(segCorte, idHueso, volumeNode)
        boneMaskAntes = boneArrAntes.astype(bool)  # astype copia; no es vista viva
        labeledAntes, nMadres = ndimage.label(boneMaskAntes)
        print(f"CranioPlan: {nMadres} pieza(s) madre etiquetada(s) por vóxeles ANTES del corte.")

        # === PASO 2: construir el anillo de corte EN VÓXELES y restarlo ===
        # No se malla un tubo delgado (eso producía el 'twist' del lazo y
        # huecos por aliasing que dejaban el corte incompleto). El anillo se
        # arma barriendo la curva a lo largo de la normal del hueso y
        # dilatándola el grosor de hoja: continuo por construcción.
        anilloCorte = self._anilloDeCorteEnVoxeles(
            boneMaskAntes, curvaNode, mallaRemanente, volumeNode,
            grosorMM, margenSeguridadMM, esCerrada
        )
        if anilloCorte is None:
            print("CranioPlan DIAGNÓSTICO: no se pudo construir el anillo de corte en vóxeles.")
            segmentEditorWidget = None
            slicer.mrmlScene.RemoveNode(segEditorNode)
            slicer.mrmlScene.RemoveNode(segCorte)
            return None

        nVoxAnillo = int(np.count_nonzero(anilloCorte))
        print(
            f"CranioPlan DIAGNÓSTICO: anillo de corte = {nVoxAnillo} vóxeles de hueso a quitar "
            f"(grosor {grosorMM:.1f} mm)."
        )

        boneCut = np.logical_and(boneMaskAntes, np.logical_not(anilloCorte)).astype(np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            boneCut, segCorte, idHueso, volumeNode
        )

        # === PASO 3: descomponer el resultado del corte ===
        # Ya no se usa el Segment Editor para separar por islas: la
        # clasificación se hace sobre las máscaras de vóxeles (numpy), más
        # robusto. Para curva CERRADA se parte por geometría (prisma del lazo);
        # para ABIERTA, por identidad de piezas conexas.
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segEditorNode)

        espaciado = volumeNode.GetSpacing()
        volumenVoxelMM3 = espaciado[0] * espaciado[1] * espaciado[2]
        boneCutMask = boneCut.astype(bool)

        flapMasks = []        # una máscara booleana por flap extraído
        restanteMask = None

        if esCerrada:
            # ---- CLASIFICACIÓN GEOMÉTRICA POR EL PRISMA DEL LAZO ----
            # El flap es el hueso ENCERRADO por el lazo (dentro del polígono
            # proyectado y dentro de ±slab del plano). Partición pura por
            # geometría: NO depende de que el anillo haya desconectado el disco
            # por conectividad (que fallaba cerca de agujeros dejando puentes).
            # Por eso SIEMPRE separa el flap del resto.
            puntosLazo = curvaNode.GetCurvePointsWorld()
            loopRAS = None
            if puntosLazo is not None and puntosLazo.GetNumberOfPoints() >= 3:
                loopRAS = np.array([puntosLazo.GetPoint(i)
                                    for i in range(puntosLazo.GetNumberOfPoints())])

            if loopRAS is None:
                print("CranioPlan DIAGNÓSTICO: no pude leer los puntos del lazo; todo al restante.")
                restanteMask = boneCutMask
            else:
                flapMask, restMask = self._particionPorPrismaDelLazo(
                    boneCutMask, loopRAS, volumeNode
                )
                nFlap = int(np.count_nonzero(flapMask)) if flapMask is not None else 0
                nRest = int(np.count_nonzero(restMask)) if restMask is not None else 0
                print(
                    "CranioPlan DIAGNÓSTICO: partición geométrica por el prisma del lazo -> "
                    f"flap {nFlap * volumenVoxelMM3 / 1000.0:.2f} cm3 ({nFlap} vóx), "
                    f"restante {nRest * volumenVoxelMM3 / 1000.0:.2f} cm3."
                )
                UMBRAL_FLAP_VOXELES = 30
                if flapMask is not None and nFlap >= UMBRAL_FLAP_VOXELES:
                    flapMasks.append(flapMask)
                    restanteMask = restMask
                else:
                    print("CranioPlan DIAGNÓSTICO: casi no hay hueso dentro del lazo; 0 flaps.")
                    restanteMask = boneCutMask
        else:
            # ---- CURVA ABIERTA: identidad por conectividad ----
            # Una línea abierta solo separa hueso si cruza una placa de lado a
            # lado. Se etiquetan las piezas conexas tras el corte y se comparan
            # con las madre: si una madre se partió en >=2 piezas reales, la
            # mayor queda en el restante y las otras son flaps.
            labeledCut, nHijas = ndimage.label(boneCutMask)
            restanteMask = np.zeros_like(boneCutMask)
            UMBRAL_RUIDO_VOXELES = 30
            porMadre = {}
            for lab in range(1, nHijas + 1):
                m = (labeledCut == lab)
                vox = int(np.count_nonzero(m))
                if vox < UMBRAL_RUIDO_VOXELES:
                    restanteMask = np.logical_or(restanteMask, m)  # ruido -> restante
                    continue
                etiquetas = labeledAntes[m]
                etiquetas = etiquetas[etiquetas > 0]
                madre = int(np.bincount(etiquetas).argmax()) if etiquetas.size else 0
                porMadre.setdefault(madre, []).append((vox, m))

            print("CranioPlan DIAGNÓSTICO: resultado del corte por pieza madre:")
            for madre in sorted(porMadre.keys()):
                piezas = sorted(porMadre[madre], key=lambda x: x[0], reverse=True)
                reales = [(v, m) for (v, m) in piezas
                          if v * volumenVoxelMM3 >= volumenMinimoFragmentoMM3]
                esquirlas = [(v, m) for (v, m) in piezas
                             if v * volumenVoxelMM3 < volumenMinimoFragmentoMM3]
                if len(reales) <= 1:
                    for (v, m) in piezas:
                        restanteMask = np.logical_or(restanteMask, m)
                    detalle = ", ".join(f"{v * volumenVoxelMM3 / 1000.0:.2f}" for (v, m) in piezas)
                    print(f"    madre {madre}: NO separada ({len(reales)} real; {detalle} cm3) -> restante.")
                else:
                    restanteMask = np.logical_or(restanteMask, reales[0][1])
                    for (v, m) in reales[1:]:
                        flapMasks.append(m)
                    for (v, m) in esquirlas:
                        restanteMask = np.logical_or(restanteMask, m)
                    vols = ", ".join(f"{v * volumenVoxelMM3 / 1000.0:.2f}" for (v, m) in reales)
                    print(
                        f"    madre {madre}: SEPARADA en {len(reales)} reales ({vols} cm3). "
                        f"Mayor al restante; {len(reales) - 1} extraída(s)."
                    )

        # === PASO 4: exportar cada máscara (restante y flaps) a un modelo ===
        def _mascaraAModelo(mascara, nombre, color):
            if mascara is None or not mascara.any():
                return None
            segTmp = slicer.mrmlScene.AddNewNodeByClass(
                'vtkMRMLSegmentationNode', 'CranioPlan_tmp_export'
            )
            segTmp.CreateDefaultDisplayNodes()
            segTmp.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
            sid = segTmp.GetSegmentation().AddEmptySegment("m")
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                mascara.astype(np.uint8), segTmp, sid, volumeNode
            )
            segTmp.CreateClosedSurfaceRepresentation()
            pd = vtk.vtkPolyData()
            segTmp.GetClosedSurfaceRepresentation(sid, pd)
            ok = pd.GetNumberOfPoints() > 0
            pdcopy = vtk.vtkPolyData()
            if ok:
                pdcopy.DeepCopy(pd)
            slicer.mrmlScene.RemoveNode(segTmp)
            if not ok:
                return None
            pdcopy = self._decimarMalla(pdcopy, reduccionMallaFragmentos)
            return self.crearModeloDesdePolyData(pdcopy, nombre, color)

        slicer.mrmlScene.RemoveNode(segCorte)

        restanteModel = _mascaraAModelo(restanteMask, "Craneo_restante", (0.9, 0.8, 0.6))
        if restanteModel is None:
            print("CranioPlan DIAGNÓSTICO: el cráneo restante quedó vacío tras el corte.")
            return None

        COLORES_EXTRAIDOS = [
            (0.55, 0.75, 0.85),  # celeste
            (0.75, 0.85, 0.55),  # verde claro
            (0.85, 0.55, 0.75),  # rosado
            (0.95, 0.75, 0.45),  # naranja suave
        ]
        fragmentosExtraidos = []
        for k, fmask in enumerate(flapMasks):
            indice = indiceInicialFragmento + k
            nombre = f"Fragmento_extraido_{indice}"
            color = COLORES_EXTRAIDOS[(indice - 1) % len(COLORES_EXTRAIDOS)]
            modelo = _mascaraAModelo(fmask, nombre, color)
            if modelo is not None:
                fragmentosExtraidos.append(modelo)
                print(f"CranioPlan: {modelo.GetName()} creado.")

        print(
            f"CranioPlan: corte terminado. {len(fragmentosExtraidos)} fragmento(s) "
            "extraído(s) por este corte."
        )
        print("=" * 60)

        return {
            "restante": restanteModel,
            "fragmentos": fragmentosExtraidos,
            "piezasCreadas": len(fragmentosExtraidos),
            "piezasAntes": nMadres,
        }


# ============================================================
# ORQUESTADOR DE CONSOLA
# ============================================================

def _volumen_de_referencia():
    """Volumen CT de referencia. Usa el del nodo de segmentación si existe;
    si no, el último ScalarVolume cargado."""
    seg = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_SEGMENTACION)
    if seg is not None and seg.IsA("vtkMRMLSegmentationNode"):
        rol = slicer.vtkMRMLSegmentationNode.GetReferenceImageGeometryReferenceRole()
        vol = seg.GetNodeReference(rol)
        if vol is not None:
            return vol
    vols = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
    return vols[-1] if vols else None


def _malla_craneo_inicial():
    """Devuelve el polydata del cráneo a cortar. Prioriza el Model
    'Craneo_Final'; si no está, lo exporta desde el segmento homónimo."""
    m = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MODELO_CRANEO)
    if (m is not None and m.IsA("vtkMRMLModelNode") and m.GetPolyData()
            and m.GetPolyData().GetNumberOfPoints() > 0):
        pd = vtk.vtkPolyData()
        pd.DeepCopy(m.GetPolyData())
        return pd

    seg = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_SEGMENTACION)
    if seg is not None and seg.IsA("vtkMRMLSegmentationNode"):
        segId = seg.GetSegmentation().GetSegmentIdBySegmentName(NOMBRE_SEGMENTO_HUESO)
        if segId:
            seg.CreateClosedSurfaceRepresentation()
            pd = vtk.vtkPolyData()
            seg.GetClosedSurfaceRepresentation(segId, pd)
            if pd.GetNumberOfPoints() > 0:
                out = vtk.vtkPolyData()
                out.DeepCopy(pd)
                return out
    return None


def _curvas_de_la_escena():
    """Todas las curvas de markups: primero cerradas, después abiertas."""
    cerradas = list(slicer.util.getNodesByClass('vtkMRMLMarkupsClosedCurveNode'))
    abiertas = [n for n in slicer.util.getNodesByClass('vtkMRMLMarkupsCurveNode')
                if not n.IsA('vtkMRMLMarkupsClosedCurveNode')]
    return cerradas + abiertas


def _limpiar_salidas_previas():
    for nombre in ["Craneo_restante"]:
        n = slicer.mrmlScene.GetFirstNodeByName(nombre)
        while n is not None:
            slicer.mrmlScene.RemoveNode(n)
            n = slicer.mrmlScene.GetFirstNodeByName(nombre)
    for n in list(slicer.util.getNodesByClass('vtkMRMLModelNode')):
        if n.GetName().startswith("Fragmento_extraido"):
            slicer.mrmlScene.RemoveNode(n)


def cortar():
    """Corta el cráneo con TODAS las curvas de la escena, en orden."""
    volumeNode = _volumen_de_referencia()
    if volumeNode is None:
        print("ERROR: no encuentro el volumen CT de referencia. Carga la tomografía primero.")
        return
    malla = _malla_craneo_inicial()
    if malla is None:
        print(f"ERROR: no encuentro el cráneo ('{NOMBRE_MODELO_CRANEO}' como Model ni como "
              f"segmento en '{NOMBRE_SEGMENTACION}'). Corre confirmar_craneo() + enviar_a_planner().")
        return
    curvas = _curvas_de_la_escena()
    if not curvas:
        print("ERROR: no hay ninguna curva de corte en la escena. Dibuja al menos una en Markups.")
        return

    print(f"Cortando con {len(curvas)} curva(s): {[c.GetName() for c in curvas]}")
    _limpiar_salidas_previas()

    logica = LogicaCorteOsteotomia()
    nFrag = 0
    for curva in curvas:
        respaldo = vtk.vtkPolyData()
        respaldo.DeepCopy(malla)

        # borrar el Craneo_restante previo para que el nombre quede limpio
        viejo = slicer.mrmlScene.GetFirstNodeByName("Craneo_restante")
        if viejo is not None and viejo.IsA("vtkMRMLModelNode"):
            slicer.mrmlScene.RemoveNode(viejo)

        res = logica.generarOsteotomia(
            respaldo, curva, volumeNode,
            indiceInicialFragmento=nFrag + 1,
            grosorMM=GROSOR_CORTE_MM,
            margenSeguridadMM=MARGEN_SEGURIDAD_MM,
            volumenMinimoFragmentoMM3=VOL_MIN_FRAGMENTO_MM3,
            reduccionMallaFragmentos=REDUCCION_MALLA,
        )
        if res is None:
            print(f"  '{curva.GetName()}': el corte no se pudo calcular; se continúa con las demás.")
            continue
        nFrag += res["piezasCreadas"]
        if res["restante"] is not None and res["restante"].GetPolyData() is not None:
            malla = res["restante"].GetPolyData()

    # ocultar el cráneo original para ver los fragmentos
    m = slicer.mrmlScene.GetFirstNodeByName(NOMBRE_MODELO_CRANEO)
    if m is not None and m.IsA("vtkMRMLModelNode") and m.GetDisplayNode():
        m.GetDisplayNode().SetVisibility(False)

    print(f"Listo. {nFrag} fragmento(s) extraído(s) en total. "
          "Están como 'Fragmento_extraido_N' + 'Craneo_restante' en el panel Data.")


print("Bloque F (corte) - script de consola cargado.")
print(f"  grosor de corte: {GROSOR_CORTE_MM} mm   margen: {MARGEN_SEGURIDAD_MM} mm")
print("Comando:  cortar()   (usa todas las curvas dibujadas sobre Craneo_Final)")
