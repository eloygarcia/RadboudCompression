#include <iostream>
#include <stdio.h>
using std::string;
//#include "internalMesh.h"

// Classe mesh
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

// Patron command
class CommandTranslacion
{
    protected:
        CommandTranslacion(){};
    public:
        virtual ~CommandTranslacion(){};
        virtual void Ejecutar() const {};
};

// Ahora mismo se produce un Aliasing de memoria !!!
// Podemos hacer una copia y devolverla
class Translacion : public CommandTranslacion
{
    private:
        float * translation_;
        concreteMesh * mesh_;
        concreteMesh * mesh1;
    public:
        explicit Translacion( float * translation, internalMesh * mesh):translation_(translation),mesh_(mesh){};
        void Ejecutar() const override {
            internalMesh * mesh1 = mesh_->Clone();
            float * temp = mesh1->getPoints();
            for(int i=0; i<3; i++) temp[i]+=translation_[i];
            mesh1->setPoints(temp);
        };
};

// main
int main()
{
    float points[3] = {3.0, 2.0, 1.0};
    int elements[3] = {5,5,5};
    string m ="mesh";

    concreteMesh * mesh = new concreteMesh(m,1);
    mesh->setPoints(points);
    mesh->setElements(elements);

    // Print mesh1
    float * pt1 = mesh->getPoints();
    int * el1 = mesh->getElements();

    std::cout << "Points: [" << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "] " << std::endl;
    std::cout << "Elements: [" << el1[0] << ", " << el1[1] << ", " << el1[2] << "] " << std::endl;

    float trs[3] = {-5.0,14.0,5.0};
    CommandTranslacion * translacion = new Translacion(trs,mesh);
    //Translacion * translacion = new Translacion(trs, mesh);
    //std::cout << translacion << std::endl;
    translacion->Ejecutar();

    // Print mesh1
    float * pt2 = mesh->getPoints();
    int * el2 = mesh->getElements();

    std::cout << "Points: [" << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << "] " << std::endl;
    std::cout << "Elements: [" << el2[0] << ", " << el2[1] << ", " << el2[2] << "] " << std::endl;

    return 0;
}