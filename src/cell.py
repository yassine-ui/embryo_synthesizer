import vtk
import numpy as np

class Cell:
    def __init__(self, a, b, c, x, y, z, t, u, v, new_x=None, new_y=None, new_z=None, new_t=None, new_u=None,
                 new_v=None, is_alive=True, velocity=None, acceleration=None):
        self.red = np.random.randn()
        self.green = np.random.randn()
        self.blue = np.random.randn()
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

    def get_actor(self):
        # create source
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(100)
        sphere.SetPhiResolution(100)

        #sphere.SetCenter(self.x - (0 - (self.a / 2)), self.y - (0 - (self.b / 2)), self.z - (0 - (self.c / 2)))
        transform = vtk.vtkTransform()
        transform.Translate(self.x, self.y, self.z)
        transform.RotateX(self.t)
        transform.RotateY(self.u)
        transform.RotateZ(self.v)
        transform.Scale(self.a*2, self.b*2, self.c*2)

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
        actor.GetProperty().SetColor(self.red, self.green, self.blue)

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
        cutActor.GetProperty().SetColor(self.red, self.green, self.blue)
        cutActor.GetProperty().SetEdgeColor(self.red, self.green, self.blue)

        #cutActor.GetProperty().SetLineWidth(1)
        #cutActor.GetProperty().EdgeVisibilityOn()
        # cutActor.GetProperty().SetOpacity(0.7)
        cutActor.SetMapper(cutMapper)

        return actor, cutActor
