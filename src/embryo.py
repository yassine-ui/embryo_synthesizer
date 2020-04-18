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
        print(semi_axe)

    def get_poly_data(self, sphere, transform):
        transform.Push()
        transform.Scale(self.semi_axe, self.semi_axe, self.semi_axe)
        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInputConnection(sphere.GetOutputPort())
        transformer.SetTransform(transform)
        transformer.Update()
        transform.Pop()
        poly_data = transformer.GetOutput()
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetNumberOfTuples(poly_data.GetNumberOfCells())
        r = np.random.randint(256)
        g = np.random.randint(256)
        b = np.random.randint(256)
        for c in range(poly_data.GetNumberOfCells()):
            colors.SetTuple(c, [r, g, b])
        poly_data.GetCellData().SetScalars(colors)
        return poly_data

    def append_poly_data(self, sphere, append_filter):
        transform = vtk.vtkTransform()
        append_filter.AddInputData(self.get_poly_data(sphere, transform))
        for i in range(self.num_cells):
            append_filter.AddInputData(self.cells[i].get_poly_data(sphere, transform))
