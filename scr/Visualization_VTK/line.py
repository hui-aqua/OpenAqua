import os
import json
import numpy as np
import vtk


def make_line(point_list,line_list):
    # ug is the container to store the grid object
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    line = vtk.vtkLine()   # 2-point elements
    
    
   # add points
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements
    for j in range(len(line_list)):
        for i in range(2):
            line.GetPointIds().SetId(i, line_list[j][i])
        ug.InsertNextCell(line.GetCellType(), line.GetPointIds())    
    return ug


def save_lines(point_dic:dict,
               elements:list,
               component:str):
    i=0
    writer = vtk.vtkUnstructuredGridWriter()    
    if os.path.exists(os.path.join(os.getcwd(), 'VTKlines_'+component)):
        print('\nVTKlines_'+component+' exists in the current folder\nPlease remove it and try again.')
        exit()
    else:
        os.mkdir(os.path.join(os.getcwd(), 'VTKlines_'+component))
        
    for k in point_dic.keys():
        uGrids=make_line(point_dic[k],elements[component])
        writer.SetFileName(os.path.join('VTKlines_'+component,'lineVTK_'+component+str(i)+'.vtk'))
        writer.SetInputData(uGrids)
        writer.Write()
        i+=1

    
    

if __name__ == '__main__':
    with open(os.path.join('..', 'Mesh_cage', 'Nodes_position.json'), 'r') as f:
        point_position = json.load(f)
    f.close()
    with open(os.path.join('..','Mesh_cage','element.json'), 'r') as f:
       element = json.load(f)
    f.close()

    save_lines(point_position,element,'rope')
    save_lines(point_position, element, 'pipes_top')

