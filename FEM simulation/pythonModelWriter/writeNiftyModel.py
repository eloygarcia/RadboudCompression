import sys
sys.path.append('/aux/')

from aux.auxiliarDefinitions import NHParameters, TIParameters, elasticParameters
from aux.auxiliarDefinitions import computeNHparameters, computeTIparameters

import numpy as np
from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader
)
import argparse

import xml.etree.cElementTree as ET

material = 'NH'
global materialParameters
materialParameters = elasticParameters
anisotropyDirection = [1,0,0]
def writeNiftyModel(meshpath, outputpath,thickness, gravity, offset, material, anisotropyRatio):
    global materialParameters
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(meshpath)
    reader.Update()

    ## writing xml
    root = ET.Element("Model")

    ## System Params
    doc = ET.SubElement(root, "SystemParams")
    ET.SubElement(doc, "TimeStep").text = "0.00001"
    ET.SubElement(doc, "TotalTime").text = "1"
    ET.SubElement(doc, "DampingCoeff").text = "0.065"
    ET.SubElement(doc, "Density").text = "0.000001" # "1"
    ET.SubElement(doc, "DoDeformableCollision").text = "0"

    ## Nodes
    nodes = ET.SubElement(root, "Nodes")
    nodes.set('DOF',str(3))
    nodes.set('NumNodes', str( reader.GetOutput().GetNumberOfPoints()))
    nodes.text = " ".join([str(round(item,3)) for item in np.array(reader.GetOutput().GetPoints().GetData()).flatten().tolist()])

    ## Elements
    elements = ET.SubElement(root, "Elements")
    elements.set('Type', 'T4')
    elements.set('NumEls', str( reader.GetOutput().GetNumberOfCells()))
    temp = np.array(reader.GetOutput().GetCells().GetData()).reshape((reader.GetOutput().GetNumberOfCells(), 5))
    elements.text = " ".join([str(item) for item in temp[:,1:].flatten().tolist()])

    ## Elements Sets
    listOfElementSets = []
    listOfMaterials = []
    listOfParameters = []

    temp_elementSets = np.array(reader.GetOutput().GetCellData().GetArray('materials'))

    for i in range(len(np.unique(temp_elementSets))):
        listOfElementSets.append(ET.SubElement(root, "ElementSet"))

        tempi_ele = np.where(temp_elementSets == np.unique(temp_elementSets)[i])[0]
        listOfElementSets[i].set("Size", str(len(tempi_ele)))
        listOfElementSets[i].text = " ".join([str(item) for item in tempi_ele])

        listOfMaterials.append( ET.SubElement( listOfElementSets[i], "Material"))
        listOfMaterials[i].set("Type",material)

        if material == 'NH':
            materialParameters = NHParameters
        elif material == 'TI':
            materialParameters = computeTIparameters(elasticParameters, anisotropyRatio, anisotropyDirection)

        listOfParameters.append( ET.SubElement(listOfMaterials[i], "ElasticParams"))
        listOfParameters[i].set( "NumParams", str( len(materialParameters[0])))
        listOfParameters[i].text = " ".join([str(item) for item in materialParameters[i] ])

    ## Boundary conditions
    bc = ET.SubElement(root,"Constraint")
    temp_BC = np.where(np.array(reader.GetOutput().GetPointData().GetArray('boundaryConditions')))[0]
    bc.set('DOF',str(0))
    # bc.set('DOF',"all")
    bc.set('Type', 'Fix')
    bc.set('NumNodes', str( len(temp_BC)))

    bc_nodes = ET.SubElement(bc,"Nodes")
    bc_nodes.text = " ".join([str(item) for item in temp_BC])

    ## Gravity
    grav = ET.SubElement(root,"Constraint")
    #grav.set("LoadShape", "POLY345")
    grav.set("LoadShape", "RAMP") # linear
    grav.set('NumNodes', str( reader.GetOutput().GetNumberOfPoints()))
    grav.set("Type","Gravity")

    grav_nodes = ET.SubElement(grav, "Nodes")
    grav_nodes.text = " ".join([str(item) for item in range( reader.GetOutput().GetNumberOfPoints())])

    acc_mag = ET.SubElement(grav, "AccelerationMagnitude")
    acc_mag.text = str(gravity)

    acc_dir = ET.SubElement(grav, "AccelerationDirection")
    acc_dir.text = " ".join([str(item) for item in np.array([0,1,0])])

    ## Contact Plates

    listOfContactPlates = []
    listOfVertex_a = []
    listOfVertex_b = []
    listOfVertex_c = []
    listOfDisplacements = []
    listOfSlavesNodes = []

    bbox = reader.GetOutput().GetPoints().GetBounds()
    # offset = 20
    # thickness = 70

    va = []
    va.append([bbox[0]-500, bbox[3], bbox[4]-500])
    va.append([bbox[0]-500, bbox[2], bbox[4]-500])

    vb = []
    vb.append([bbox[1]+500, bbox[3], bbox[4]-500])
    vb.append([bbox[0]-500, bbox[2], bbox[5]+500])

    vc = []
    vc.append([bbox[0]-500, bbox[3], bbox[5]+500])
    vc.append([bbox[1]+500, bbox[2], bbox[4]-500])

    disp = []
    disp.append([0, -offset, 0])
    disp.append([0, (bbox[3]-(bbox[2]+offset))-thickness ,0])
    #disp.append([0, offset,0])
    for i in range(2):
        listOfContactPlates.append( ET.SubElement(root, "ContactPlate"))

        listOfVertex_a.append( ET.SubElement(listOfContactPlates[i], "a"))
        listOfVertex_a[i].text = " ".join([str(round(item,3)) for item in va[i]])

        listOfVertex_b.append(ET.SubElement(listOfContactPlates[i], "b"))
        listOfVertex_b[i].text = " ".join([str(round(item,3)) for item in vb[i]])

        listOfVertex_c.append(ET.SubElement(listOfContactPlates[i], "c"))
        listOfVertex_c[i].text = " ".join([str(round(item,3)) for item in vc[i]])

        listOfDisplacements.append(ET.SubElement(listOfContactPlates[i], "Disp"))
        listOfDisplacements[i].text = " ".join([str(round(item,3)) for item in disp[i]])

        listOfSlavesNodes.append( ET.SubElement(listOfContactPlates[i], "SlvNodes"))
        listOfSlavesNodes[i].set("NumNodes", str(reader.GetOutput().GetNumberOfPoints()))
        listOfSlavesNodes[i].text =  " ".join([str(item) for item in range( reader.GetOutput().GetNumberOfPoints())])


    ## Write xml
    tree = ET.ElementTree(root)
    tree.write(outputpath)

if __name__ == "__main__":
    import sys