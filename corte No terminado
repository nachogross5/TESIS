import numpy as np
import vtk

# ==== CONFIGURACION - AJUSTAR ANTES DE CORRER ====
nombreCurva = "CC"
nombreCraneo = "Craneo_Final"
anchoCorte = 2.0
profundidadCorte = 15.0
radioBusquedaNormal = 3.0
mSamplesPorSegmento = 5   # muestras a lo largo de CADA segmento (entre tus puntos) para promediar la normal de esa plancha

curvaNode = slicer.util.getNode(nombreCurva)
craneoNode = slicer.util.getNode(nombreCraneo)
craneoPoly = craneoNode.GetPolyData()

normalsFilter = vtk.vtkPolyDataNormals()
normalsFilter.SetInputData(craneoPoly)
normalsFilter.ComputePointNormalsOn()
normalsFilter.ComputeCellNormalsOff()
normalsFilter.SplittingOff()
normalsFilter.ConsistencyOn()
normalsFilter.AutoOrientNormalsOn()
normalsFilter.Update()
craneoConNormales = normalsFilter.GetOutput()
normalesArray = craneoConNormales.GetPointData().GetNormals()

locator = vtk.vtkPointLocator()
locator.SetDataSet(craneoConNormales)
locator.BuildLocator()

# SOLO los puntos de control (los que vos marcaste a mano)
nPuntos = curvaNode.GetNumberOfControlPoints()
puntosCurva = []
for i in range(nPuntos):
    p = [0.0, 0.0, 0.0]
    curvaNode.GetNthControlPointPositionWorld(i, p)
    puntosCurva.append(np.array(p))
print("Puntos de control usados:", nPuntos)

def normalPromedioEnPunto(p, radio):
    idList = vtk.vtkIdList()
    locator.FindPointsWithinRadius(radio, p, idList)
    if idList.GetNumberOfIds() == 0:
        idList.InsertNextId(locator.FindClosestPoint(p))
    acumulado = np.zeros(3)
    for j in range(idList.GetNumberOfIds()):
        acumulado += np.array(normalesArray.GetTuple3(idList.GetId(j)))
    norma = np.linalg.norm(acumulado)
    if norma < 1e-6:
        return np.array(normalesArray.GetTuple3(locator.FindClosestPoint(p)))
    return acumulado / norma

def normalPromedioSegmento(p1, p2, muestras, radio):
    acumulado = np.zeros(3)
    for k in range(muestras):
        t = k / float(muestras - 1) if muestras > 1 else 0.5
        acumulado += normalPromedioEnPunto(p1 + t * (p2 - p1), radio)
    norma = np.linalg.norm(acumulado)
    return acumulado / norma if norma > 1e-6 else acumulado

# Una normal por SEGMENTO (entre punto i y punto i+1) - loop cerrado, wraparound
normalesSegmento = []
for i in range(nPuntos):
    p1 = puntosCurva[i]
    p2 = puntosCurva[(i+1) % nPuntos]
    normalesSegmento.append(normalPromedioSegmento(p1, p2, mSamplesPorSegmento, radioBusquedaNormal))

# Normal de union en cada punto = promedio de las 2 planchas que se tocan ahi (NO es suavizado de la curva, es solo para que las 2 planchas compartan el borde)
normalesJunta = []
for i in range(nPuntos):
    suma = normalesSegmento[i-1] + normalesSegmento[i]  # segmento anterior + segmento siguiente
    norma = np.linalg.norm(suma)
    normalesJunta.append(suma / norma if norma > 1e-6 else normalesSegmento[i])

# Tangente en cada punto (wraparound, loop cerrado)
tangentes = []
for i in range(nPuntos):
    anterior = puntosCurva[i-1]
    siguiente = puntosCurva[(i+1) % nPuntos]
    t1 = puntosCurva[i] - anterior
    t2 = siguiente - puntosCurva[i]
    n1 = np.linalg.norm(t1)
    n2 = np.linalg.norm(t2)
    t = (t1/n1 if n1 > 1e-6 else t1) + (t2/n2 if n2 > 1e-6 else t2)
    normaT = np.linalg.norm(t)
    tangentes.append(t/normaT if normaT > 1e-6 else np.array([1.0,0.0,0.0]))

# w = direccion de ancho, con correccion de continuidad de signo (evita el twist)
wCrudos = []
for i in range(nPuntos):
    w = np.cross(normalesJunta[i], tangentes[i])
    normaW = np.linalg.norm(w)
    w = w / normaW if normaW > 1e-6 else np.array([1.0, 0.0, 0.0])
    wCrudos.append(w)

wCorregidos = [wCrudos[0]]
for i in range(1, nPuntos):
    wActual = wCrudos[i]
    if np.dot(wActual, wCorregidos[i-1]) < 0:
        wActual = -wActual
    wCorregidos.append(wActual)
if np.dot(wCorregidos[-1], wCorregidos[0]) < 0:
    print("ADVERTENCIA: el loop tiene un twist impar - revisar geometria de la curva")

# Anillos SOLO en tus puntos de control (nPuntos anillos, nPuntos planchas)
puntosHerramienta = vtk.vtkPoints()
celdas = vtk.vtkCellArray()
anillos = []
for i in range(nPuntos):
    n = normalesJunta[i]
    w = wCorregidos[i]
    c = puntosCurva[i]
    p1 = c + (anchoCorte/2.0)*w + (profundidadCorte/2.0)*n
    p2 = c - (anchoCorte/2.0)*w + (profundidadCorte/2.0)*n
    p3 = c - (anchoCorte/2.0)*w - (profundidadCorte/2.0)*n
    p4 = c + (anchoCorte/2.0)*w - (profundidadCorte/2.0)*n
    idx1 = puntosHerramienta.InsertNextPoint(p1.tolist())
    idx2 = puntosHerramienta.InsertNextPoint(p2.tolist())
    idx3 = puntosHerramienta.InsertNextPoint(p3.tolist())
    idx4 = puntosHerramienta.InsertNextPoint(p4.tolist())
    anillos.append((idx1, idx2, idx3, idx4))

for i in range(nPuntos):
    a = anillos[i]
    b = anillos[(i+1) % nPuntos]
    for k in range(4):
        k2 = (k+1) % 4
        quad = vtk.vtkPolygon()
        quad.GetPointIds().SetNumberOfIds(4)
        quad.GetPointIds().SetId(0, a[k])
        quad.GetPointIds().SetId(1, a[k2])
        quad.GetPointIds().SetId(2, b[k2])
        quad.GetPointIds().SetId(3, b[k])
        celdas.InsertNextCell(quad)

herramientaPoly = vtk.vtkPolyData()
herramientaPoly.SetPoints(puntosHerramienta)
herramientaPoly.SetPolys(celdas)

triangulador = vtk.vtkTriangleFilter()
triangulador.SetInputData(herramientaPoly)
triangulador.Update()

limpiador = vtk.vtkCleanPolyData()
limpiador.SetInputData(triangulador.GetOutput())
limpiador.Update()
herramientaFinal = limpiador.GetOutput()

herramientaNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "HerramientaCorte")
herramientaNode.SetAndObservePolyData(herramientaFinal)
herramientaNode.CreateDefaultDisplayNodes()
herramientaNode.GetDisplayNode().SetColor(1.0, 0.6, 0.0)
herramientaNode.GetDisplayNode().SetOpacity(0.5)

featureEdges = vtk.vtkFeatureEdges()
featureEdges.SetInputData(herramientaFinal)
featureEdges.BoundaryEdgesOn()
featureEdges.NonManifoldEdgesOn()
featureEdges.FeatureEdgesOff()
featureEdges.ManifoldEdgesOff()
featureEdges.Update()
numBordesLibres = featureEdges.GetOutput().GetNumberOfCells()
print("Bordes libres/no-manifold en HerramientaCorte:", numBordesLibres)
if numBordesLibres > 0:
    print("ATENCION: la herramienta no quedo cerrada. Revisar antes de confiar en los fragmentos.")

booleanFiltro = vtk.vtkBooleanOperationPolyDataFilter()
booleanFiltro.SetOperationToDifference()
booleanFiltro.SetInputData(0, craneoPoly)
booleanFiltro.SetInputData(1, herramientaFinal)
booleanFiltro.Update()
resultadoCorte = booleanFiltro.GetOutput()

conectividad = vtk.vtkPolyDataConnectivityFilter()
conectividad.SetInputData(resultadoCorte)
conectividad.SetExtractionModeToAllRegions()
conectividad.ColorRegionsOn()
conectividad.Update()
numRegiones = conectividad.GetNumberOfExtractedRegions()
print("Fragmentos generados:", numRegiones)

for r in range(numRegiones):
    conectividad.SetExtractionModeToSpecifiedRegions()
    conectividad.InitializeSpecifiedRegionList()
    conectividad.AddSpecifiedRegion(r)
    conectividad.Update()
    limpiadorFrag = vtk.vtkCleanPolyData()
    limpiadorFrag.SetInputConnection(conectividad.GetOutputPort())
    limpiadorFrag.Update()
    fragmentoNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Fragmento_%d" % (r+1))
    fragmentoNode.SetAndObservePolyData(limpiadorFrag.GetOutput())
    fragmentoNode.CreateDefaultDisplayNodes()

print("Listo.")
