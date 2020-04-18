import sys
import vtk
import numpy as np
import math
import os
import pathlib
import subprocess

from Quadric import CloseWin
from embryo import Embryo

class Simulator:
    def __init__(self):
        self.colors = vtk.vtkNamedColors()
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(self.colors.GetColor3d('Silver'))
        self.renwin = vtk.vtkRenderWindow()
        self.renwin.SetWindowName("CallBack")
        self.renwin.AddRenderer(self.renderer)
        # An interactor
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.renwin)
        self.use_function_callback = False

        (RSphere, NbEllipsoides, Ellipsoides) = self.ReadResPacking('solution.txt')
        self.emb = Embryo(RSphere, NbEllipsoides, Ellipsoides)

        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(100)
        sphere.SetPhiResolution(50)
        append_filter = vtk.vtkAppendPolyData()
        self.emb.append_poly_data(sphere, append_filter)

        clean_filter = vtk.vtkCleanPolyData()
        clean_filter.SetInputConnection(append_filter.GetOutputPort())
        clean_filter.Update()

        # Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(clean_filter.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.3)
        self.renderer.AddActor(actor)
        self.interactor.Initialize()
        self.renwin.Render()
        self.interactor.Start()

    def ReadResPacking(self, filename):
        with open(filename, 'r') as infile:
            data = infile.read()
        my_list = data.splitlines()
        l = 0
        print(my_list[l])
        RSphere = np.float(my_list[0])
        NbEllipsoides = np.int(my_list[1])
        Ellipsoides = [float(x) for x in ' '.join(my_list[2:]).split()]
        Ellipsoides = np.reshape(Ellipsoides, (NbEllipsoides, 9))
        infile.close()
        return (RSphere, NbEllipsoides, Ellipsoides)

if __name__ == "__main__":
    sim = Simulator()