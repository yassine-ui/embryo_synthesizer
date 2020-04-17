import vtk

ColorBackground = [0.0, 0.0, 0.0]

FirstobjPath = r"Skull.obj"

reader = vtk.vtkOBJReader()

reader.SetFileName(FirstobjPath)
reader.Update()
polyData = reader.GetOutput()

colors = vtk.vtkUnsignedCharArray()
colors.SetNumberOfComponents(3)
colors.SetNumberOfTuples(polyData.GetNumberOfCells())
for c in range(polyData.GetNumberOfCells()):
    colors.SetTuple(c, [255, 0, 0])
polyData.GetCellData().SetScalars(colors)


transform = vtk.vtkTransform()
transform.Translate(20, 0, 0)
#transform.RotateWXYZ(45,0,1,0)
transformFilter=vtk.vtkTransformPolyDataFilter()
transformFilter.SetTransform(transform)
transformFilter.SetInputConnection(reader.GetOutputPort())
transformFilter.Update()
polyData2 = transformFilter.GetOutput()


colors = vtk.vtkUnsignedCharArray()
colors.SetNumberOfComponents(3)
colors.SetNumberOfTuples(polyData2.GetNumberOfCells())
for c in range(polyData2.GetNumberOfCells()):
    colors.SetTuple(c, [0, 0, 255])
polyData2.GetCellData().SetScalars(colors)


appendFilter = vtk.vtkAppendPolyData()
if vtk.VTK_MAJOR_VERSION <= 5:
    appendFilter.AddInputConnection(polyData.GetProducerPort())
    appendFilter.AddInputConnection(polyData2.GetProducerPort())
else:
    appendFilter.AddInputData(polyData)
    appendFilter.AddInputData(polyData2)

appendFilter.Update()

#  Remove any duplicate points.
cleanFilter = vtk.vtkCleanPolyData()
cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
cleanFilter.Update()

# Create a mapper and actor
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(cleanFilter.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

plane = vtk.vtkPlane()
plane.SetOrigin(0, -14, 0)
plane.SetNormal(0, 1, 0)

#create cutter
cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(appendFilter.GetOutputPort())
cutter.Update()

cutStrips = vtk.vtkStripper() ; #Forms loops (closed polylines) from cutter
cutStrips.SetInputConnection(cutter.GetOutputPort())
cutStrips.Update()
cutPoly = vtk.vtkPolyData() ; #This trick defines polygons as polyline loop
cutPoly.SetPoints((cutStrips.GetOutput()).GetPoints())
cutPoly.SetPolys((cutStrips.GetOutput()).GetLines())

cutMapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    cutMapper.SetInput(cutPoly)
else:
    cutMapper.SetInputData(cutPoly)

cutActor = vtk.vtkActor()
#cutActor.GetProperty().SetColor(1, 1, 1)
#cutActor.GetProperty().SetEdgeColor(0, 1, 0)

cutActor.GetProperty().SetLineWidth(2)
cutActor.GetProperty().EdgeVisibilityOn()
#cutActor.GetProperty().SetOpacity(0.7)
cutActor.SetMapper(cutMapper)




ren = vtk.vtkRenderer()
ren.SetBackground(ColorBackground)

renWin = vtk.vtkRenderWindow()

renWin.AddRenderer(ren)


iren = vtk.vtkRenderWindowInteractor()

iren.SetRenderWindow(renWin)

# Assign actor to the renderer

#actor.GetProperty().SetColor(1, 0, 1)

#ren.AddActor(actor)
#ren.AddActor(actor2)
ren.AddActor(cutActor)
# Enable user interface interactor

iren.Initialize()

renWin.Render()

iren.Start()
