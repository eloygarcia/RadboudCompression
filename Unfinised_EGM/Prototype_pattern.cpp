#include <iostream>
#include <stdio.h>
using std::string;
//#include "internalMesh.h"

class internalMesh
{

protected:
    int * elements_;
    float * points_;
    string name_;


public:
    void setElements( int * elements){ elements_ = elements;};
    int * getElements() {return elements_;};

    void setPoints( float * points){ points_ = points;};
    float * getPoints() {return points_;};

    internalMesh(){};
    //internalMesh( string name):name_(name){};
    internalMesh( string name):name_(name){};
    //internalMesh(const internalMesh&); //:name_(name){};

    virtual ~internalMesh(){};

    virtual internalMesh* Clone() const {return 0;};  // atención a esta línea.
    virtual void Method(){};
};


// Patrón Prototype
class concreteMesh : public internalMesh
{
protected:
    int id_;

public:
    concreteMesh(string name, float id):internalMesh(name),id_(id){};
    //concreteMesh(float id)::internalMesh(),id_(id){};
    internalMesh * Clone() const override{
        return new concreteMesh(*this);
    };
};


// main
int main()
{
    float points[3] = {1.0,1.0,1.0};
    int elements[3] = {1,2,3};
    string m ="mesh";

    //std::cout << "Points: [" << points[0] << ", " << points[1] << ", " << points[2] << "] " << std::endl;
    //std::cout << "Elements: [" << elements[0] << ", " << elements[1] << ", " << elements[2] << "] " << std::endl;

    concreteMesh * mesh = new concreteMesh(m,2);
    mesh->setPoints(points);
    mesh->setElements(elements);
    std::cout << mesh->getPoints() << std::endl;
    std::cout << mesh->getElements() << std::endl;

    internalMesh * mesh2 = mesh->Clone();
    std::cout << mesh2->getPoints() << std::endl;
    std::cout << mesh2->getElements() << std::endl;

    // Print mesh1
    float * pt1 = mesh->getPoints();
    int * el1 = mesh->getElements();

    std::cout << "Points: [" << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "] " << std::endl;
    std::cout << "Elements: [" << el1[0] << ", " << el1[1] << ", " << el1[2] << "] " << std::endl;

    // Print mesh2
    float * pt = mesh2->getPoints();
    int * el = mesh2->getElements();

    std::cout << "Points: [" << pt[0] << ", " << pt[1] << ", " << pt[2] << "] " << std::endl;
    std::cout << "Elements: [" << el[0] << ", " << el[1] << ", " << el[2] << "] " << std::endl;

    return 0;
}