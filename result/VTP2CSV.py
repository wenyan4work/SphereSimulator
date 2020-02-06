import numpy as np
import scipy as sp
import vtk
import glob as glob
import re
import os


def getFrameNumber_lambda(filename): return int(
    re.search('_(.+?).pvtp', filename).group(1))


def vtp2csv(fileIn, fileOut):
    reader = vtk.vtkGenericDataObjectReader()
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName(fileIn)
    reader.Update()

    point_obj = reader.GetOutput()
    points = point_obj.GetPoints()

    table = vtk.vtkDataObjectToTable()
    table.SetInputData(point_obj)
    table.Update()
    table.GetOutput().AddColumn(points.GetData())
    table.Update()

    writer = vtk.vtkDelimitedTextWriter()
    writer.SetInputConnection(table.GetOutputPort())
    writer.SetFileName(fileOut)
    writer.Update()
    writer.Write()


FileList = glob.glob('./result*/Sphere_*.pvtp')
# sort as numerical order
FileList.sort(key=getFrameNumber_lambda)
print(FileList)


for fileVTP in FileList:
    fileIn = fileVTP
    pre, ext = os.path.splitext(os.path.basename(fileIn))
    fileOut = './SphereCSV/'+pre+'.csv'
    vtp2csv(fileIn, fileOut)
