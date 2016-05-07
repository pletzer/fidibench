#!/usr/bin/env python

import argparse
import vtk

parser = argparse.ArgumentParser(description='Generate picture.')
parser.add_argument('--filename', dest='filename', default='',
                   help='VTK file')

args = parser.parse_args()


reader = vtk.vtkRectilinearGridReader()
cell2pt = vtk.vtkCellDataToPointData()
contour = vtk.vtkContourFilter()
mapper = vtk.vtkPolyDataMapper()
actor = vtk.vtkActor()

reader.SetFileName(args.filename)
reader.ReadAllScalarsOn()

numCont = 11
fmin, fmax = 0., 1.
contour.SetNumberOfContours(numCont)
for i in range(numCont):
    contour.SetValue(i, fmin + (fmax - fmin)*i/float(numCont - 1))
contour.GenerateTrianglesOn()
contour.ComputeGradientsOn()

mapper.ScalarVisibilityOff()

actor.GetProperty().SetColor(1.0, 0., 0.)

# Connect
pdata = reader.GetOutput()
cell2pt.SetInputConnection(reader.GetOutputPort())
contour.SetInputConnection(cell2pt.GetOutputPort())
mapper.SetInputConnection(contour.GetOutputPort())
actor.SetMapper(mapper)

mapper.Update()

#mapper.SetScalar('f')
#print 'min/max f range: ', mapper.GetRange()
#pdata = contour.GetOutput()
#print pdata

numScalars = reader.GetNumberOfScalarsInFile()
for i in range(numScalars):
    print 'scalar field: ', reader.GetScalarsNameInFile(i)

print contour.GetOutput()

# Create the renderer, the render window, and the interactor. The
# renderer draws into the render window, the interactor enables mouse-
# and keyboard-based interaction with the scene.
ren = vtk.vtkRenderer()
ren.AddActor(actor)
ren.SetBackground(0.5, 0.5, 0.5)
renWin = vtk.vtkRenderWindow()
renWin.SetSize(600, 500)
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()