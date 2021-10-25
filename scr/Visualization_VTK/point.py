import os
import json
import numpy as np
import vtk


def make_point(point_list):
    # A polyvertex is a cell represents a set of 0D vertices
    ug = vtk.vtkUnstructuredGrid()   
    points = vtk.vtkPoints()
    polyVertex = vtk.vtkPolyVertex()
    
    numberOfVertices = len(point_list)
    polyVertex.GetPointIds().SetNumberOfIds(numberOfVertices)
    
    for i in range(numberOfVertices):
        points.InsertNextPoint(point_list[i])
        polyVertex.GetPointIds().SetId(i, i)

    ug.SetPoints(points)
    ug.InsertNextCell(polyVertex.GetCellType(), polyVertex.GetPointIds())

    return ug



def save_point(point_dic):
    i=0
    writer = vtk.vtkUnstructuredGridWriter()    
    if os.path.exists(os.path.join(os.getcwd(), 'VTKpoints')):
        print('\nVTKpoints exists in the current folder\nPlease remove it and try again.')
        exit()
    else:
        os.mkdir(os.path.join(os.getcwd(), 'VTKpoints'))
        
    for k in point_dic.keys():
        uGrids=make_point(point_dic[k])
        writer.SetFileName(os.path.join('VTKpoints','pointVTK_'+str(i)+'.vtk'))
        writer.SetInputData(uGrids)
        writer.Write()
        i+=1


if __name__ == '__main__':
    with open(os.path.join('..','Mesh_cage','Nodes_position.json'), 'r') as f:
       point_position = json.load(f)
    f.close()
        
    save_point(point_position)
    
