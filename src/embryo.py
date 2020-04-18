import vtk
from cell import Cell
import numpy as np
from vtk.vtkFiltersKitPython import vtkAppendPolyData


class Embryo:
    def __init__(self, semi_axe, num_cells, cells_params):
        self.semi_axe = semi_axe
        self.num_cells = num_cells
        # self.cells_params = cells_params
        self.cells = [None] * num_cells

        for i in range(num_cells):
            self.cells[i] = Cell(cells_params[i][0], cells_params[i][1], cells_params[i][2], cells_params[i][3],
                                 cells_params[i][4],
                                 cells_params[i][5], cells_params[i][6], cells_params[i][7], cells_params[i][8])

    def get_actor(self):
        # create source
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(100)
        sphere.SetPhiResolution(50)

        transform = vtk.vtkTransform()
        transform.Translate(0, 0, 0)
        transform.Scale(self.semi_axe*2, self.semi_axe*2, self.semi_axe*2)

        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInputConnection(sphere.GetOutputPort())
        transformer.SetTransform(transform)
        transformer.Update()

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(transformer.GetOutputPort())

        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.3)


        plane = vtk.vtkPlane()
        plane.SetOrigin(0, 0, 0)
        plane.SetNormal(0, 0, 1)

        # create cutter
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputConnection(transformer.GetOutputPort())
        cutter.Update()

        cutStrips = vtk.vtkStripper();  # Forms loops (closed polylines) from cutter
        cutStrips.SetInputConnection(cutter.GetOutputPort())
        cutStrips.Update()

        cutPoly = vtk.vtkPolyData()  # This trick defines polygons as polyline loop
        cutPoly.SetPoints((cutStrips.GetOutput()).GetPoints())
        cutPoly.SetPolys((cutStrips.GetOutput()).GetLines())
        cutMapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            cutMapper.SetInput(cutPoly)
        else:
            cutMapper.SetInputData(cutPoly)

        cutActor = vtk.vtkActor()
        cutActor.GetProperty().SetColor(1, 1, 1)
        cutActor.GetProperty().SetEdgeColor(0, 1, 0)

        cutActor.GetProperty().SetLineWidth(1)
        cutActor.GetProperty().EdgeVisibilityOn()
        # cutActor.GetProperty().SetOpacity(0.7)
        cutActor.SetMapper(cutMapper)

        return actor, cutActor

    def set_actors(self, renderer):
        actor, cutActor = self.get_actor()
        renderer.AddActor(actor)
        renderer.AddActor(cutActor)
        for i in range(self.num_cells):
            actor, cutActor = self.cells[i].get_actor()
            renderer.AddActor(actor)
            renderer.AddActor(cutActor)