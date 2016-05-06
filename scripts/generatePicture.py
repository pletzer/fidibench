#!/usr/bin/env python

import argparse
import vtk

parser = argparse.ArgumentParser(description='Generate picture.')
parser.add_argument('--filename', dest='filename', default='',
                   help='VTK file')

args = parser.parse_args()


reader = vtk.vtkRectilinearGridReader()
contour = vtk.vtkContourFilter()
mapper = vtk.vtkPolyDataMapper()
actor = vtk.vtkActor()

# Connect
contour.SetInputConnection(reader.GetOutputPort())
mapper.SetInputConnection(contour.GetOutputPort())
actor.SetMapper(mapper)

reader.SetFileName(args.filename)
reader.Update()

contour.SetNumberOfContours(1)
contour.SetValue(0, 0.5)

# Create the renderer, the render window, and the interactor. The
# renderer draws into the render window, the interactor enables mouse-
# and keyboard-based interaction with the scene.
ren = vtk.vtkRenderer()
ren.AddActor(actor)
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()