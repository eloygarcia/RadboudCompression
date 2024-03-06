# https://www.linkedin.com/pulse/converting-class-object-xmldictjson-string-muhammad-shoaib/
import logging
import json
from typing import Any, List, Dict
from xml.etree import ElementTree as et


LOG = logging.getLogger(__name__)


from dataclasses import dataclass
import numpy as np
from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader
)

class writerNiftySim:
    def __init__(self):
        self.outputputfilename=None
        self.f = 0


    def setFilename(self, filename):
        self.outputputfilename = filename

    def setVtkFEMmesh(self, filename):
        reader = vtkUnstructuredGridReader()
        reader.SetFileName( filename )
        reader.Update()


## https://docs.python.org/3/tutorial/classes.html
# https://www.linkedin.com/pulse/converting-class-object-xmldictjson-string-muhammad-shoaib/
@dataclass
class systemParameters:
    timeStep: float = 0.0005
    totalTime: float = 10
    DampgingCoefficient: float = 0.65
    density: int = 1
    DoDeformabeCollision: int = 0

@dataclass
class nodesInfo:
    numberOfNodes: int
    nodes : np.array([], dtype='float')
    degreesOfFreedom: int = 3

@dataclass
class elementsInfo:
    numberOfElements: int
    elements: np.array([], dtype='int')
    typeOfelements: str ='T4'

@dataclass
class material:
    elasticParams: np.array([], dtype='float')
    type: str = 'NH'
    numberOfElasticParams: int = 0

@dataclass
class elementSet:
    size: int
    elements: np.array([], dtype='int')
    material: material # Type material

@dataclass
class boundaryConditions:
    degreeOfFreedom: int
    numberOfNodes: int
    type: str
    listofnodes: list

@dataclass
class contactPlate:
    a: list
    b: list
    c: list
    displacement: list
    numberOfNodes: int
    slaveNodes: list

#@dataclass
class gravity:
    numberOfNodes: int
    nodes: list
    accelerationDirection: list = list([0,1,0])
    type: str ="Gravity"
    loadShape: str = "POLY345"
    accelerationMagnitud: int = 150


@dataclass
class FEMmodel:
    nodesinformation: nodesInfo
    elementsinformation: elementsInfo
    listOfElementsSets: list
    listOfBoundaryConditions: list
    listOfPaddels:list
    gravityInformation: gravity
    parameters: systemParameters


""" All classes must be drived from the following base class.
"""


# https://www.linkedin.com/pulse/converting-class-object-xmldictjson-string-muhammad-shoaib/
class XmlElement:
    """Base Class"""
    _NAME = None  # Name of the Class that you want to appear in the XML.
    _TEXT = None  # Text description of the class if available.


def get_members(obj: object) -> List[str]:
    """ Helper function
        Get class members other than functions and private members.
    """
    members = []
    for attr in dir(obj):
        if not callable(getattr(obj, attr)) and \
                not attr.startswith("_") and \
                not attr.startswith("__"):
            members.append(attr)
    return members


def to_xml(obj: Any) -> Any:
    """ Convert a class object and its members to XML.

    Each class member is treated as a tag to the current XML-element.
    Each member object is treated as a new sub-element.
    Each 'list' member is treated as a new list tag.
    """

    if isinstance(obj, dict):
        raise Exception("Dictionary type is not supported.")

    root = None
    tags = {}

    subelements = {}  # type: Dict[Any, Any]
    for member in get_members(obj):
        item = getattr(obj, member)
        member = member.replace('_', '-')

        # if object is None, add empty tag
        if item is None:
            subelements[member] = et.Element(member)
        else:
            if not isinstance(item, (str, XmlElement, list, set, tuple)):
                raise Exception("Attributes must be an expected type, but was: {}".format(type(item)))

            # Add list sub-elements
            if isinstance(item, (list, set, tuple)):
                subelements[member] = []
                for list_object in item:
                    subelements[member].append(to_xml(list_object))
            # Add sub-element
            elif isinstance(item, XmlElement):
                subelements[member] = to_xml(item)
            # Add element's tag name
            else:
                tags[member] = item

    try:
        if obj._NAME:
            root = et.Element(obj._NAME, tags)
        else:
            raise Exception("Name attribute can't be empty.")
    except (AttributeError, TypeError) as ex:
        print("Attribute value or type is wrong. %s: %s", obj, ex)
        raise

    # Add sub elements if any
    if subelements:
        for name, values in subelements.items():
            if isinstance(values, list):  # if list of elements. Add all sub-elements
                sub = et.SubElement(root, name)
                for value in values:
                    sub.append(value)
            else:  # single sub-child or None
                if values is None:  # if None, add empty tag with name
                    sub = et.SubElement(root, name)
                else:  # else add object
                    root.append(values)

    try:
        if obj._TEXT:
            root.text = obj._TEXT
    except AttributeError as ex:
        print("Attribute does not exists. %s: %s", obj, ex)
        raise

    return root


def object_to_xml(obj: Any) -> Any:
    """ Convert the given class object to xml document str.
    """
    return et.tostring(element=to_xml(obj), encoding="UTF-8")