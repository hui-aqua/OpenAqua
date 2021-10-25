import os
import json
import numpy as np
import vtk


def make_trigangle(point_list,triangle_list):
    # ug is the container to store the grid object
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    triangles = vtk.vtkTriangle()  # 3-point elements
    pixel=vtk.vtkPixel() # 4-point elements
    
    # add points   
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements
    for j in range(len(triangle_list)):
        if len(triangle_list[j])==3: # add 3-point elements
            for i in range(3): 
                triangles.GetPointIds().SetId(i, triangle_list[j][i])
            ug.InsertNextCell(triangles.GetCellType(), triangles.GetPointIds())    
        elif len(triangle_list[j])==4:  # add 4-point elements
            for i in range(4): 
                pixel.GetPointIds().SetId(i,triangle_list[j][i])
            ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())    
    return ug


def save_surf(point_dic,elements):
    i=0
    writer = vtk.vtkUnstructuredGridWriter()  
    if os.path.exists(os.path.join(os.getcwd(), 'VTKsurfs')):
        print('\nVTKsurfs exists in the current folder\nPlease remove it and try again.')
        exit()
    else:
        os.mkdir(os.path.join(os.getcwd(), 'VTKsurfs'))
    
    
    for k in point_dic.keys():
        uGrids=make_trigangle(point_dic[k],elements)
        writer.SetFileName(os.path.join('VTKsurfs','surfTK_'+str(i)+'.vtk'))
        writer.SetInputData(uGrids)
        writer.Write()
        i+=1  
    

if __name__ == '__main__':
    # read nodes' positions
    with open(os.path.join('..','Mesh_cage','Nodes_position.json'), 'r') as f:
       point_position = json.load(f)
    f.close()
    # read element 
    with open(os.path.join('..','Mesh_cage','element.json'), 'r') as f:
       element = json.load(f)
    f.close()

    save_surf(point_position,element['netting'])
