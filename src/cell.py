import vtk
import numpy as np

class Cell:
    def __init__(self, a, b, c, x, y, z, t, u, v, new_x=None, new_y=None, new_z=None, new_t=None, new_u=None,
                 new_v=None, is_alive=True, velocity=None, acceleration=None):
        self.a = a
        self.b = b
        self.c = c
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.u = u
        self.v = v
        self.new_x = new_x
        self.new_y = new_y
        self.new_z = new_z
        self.new_t = new_t
        self.new_u = new_u
        self.new_v = new_v
        self.is_alive = is_alive
        self.velocity = velocity
        self.acceleration = acceleration

    def set_new_position(self, new_x, new_y, new_z):
        self.new_x = new_x
        self.new_y = new_y
        self.new_z = new_z

    def set_new_rotation_angle(self, new_t, new_u, new_v):
        self.new_t = new_t
        self.new_u = new_u
        self.new_v = new_v

    def set_acceleration(self, acceleration):
        self.acceleration = acceleration

    def set_velocity(self, velocity):
        self.velocity = velocity

    def set_is_alive(self, is_alive):
        self.is_alive = is_alive

    def get_poly_data(self, sphere, transform):
        transform.Push()
        transform.Scale(self.a, self.b, self.c)
        transform.Translate(self.x, self.y, self.z)
        transform.RotateX(self.t)
        transform.RotateY(self.u)
        transform.RotateZ(self.v)
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