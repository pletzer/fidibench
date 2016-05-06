#!/usr/bin/env python

import argparse
import vtk

parser = argparse.ArgumentParser(description='Generate picture.')
parser.add_argument('filename', default='',
                   help='VTK file')

args = parser.parse_args()


reader = vtk.vtkRectilinearGridReader()
contour = vtk.vtkContourFilter()
mapper = vtk.vtkPolyDataMapper()
actor = vtk.vtkActor()

# Connect
contour.SetInputConnection(reader.GetOutputPort())
mapper.SetInputConnection(contour.GetOutputPort())

reader.SetFileName(args.filename)
reader.Update()

# Create the renderer, the render window, and the interactor. The
# renderer draws into the render window, the interactor enables mouse-
# and keyboard-based interaction with the scene.
aRenderer = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(aRenderer)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()