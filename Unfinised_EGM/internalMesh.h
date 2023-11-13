//#ifndef __internalMesh__
//#define __internalMesh__

class internalMesh
{
public:
    internalMesh(void){};
    virtual ~internalMesh(void){};

private:
    int * elements_;
    float * points_;

public:
    void setElements( int * elements){ elements_ = elements;};
    int * getElements() {return elements_;};

    void setPoints( float * points){ points = points_;};
    float * getPoints() {return points_;};

    virtual internalMesh *Clone() const = 0;
};

// Patrón Prototype
class concreteMesh : public internalMesh
{
private:
    int id_;

public:
    internalMesh *Clone() const override{
        return new concreteMesh(*this);
    }
};

// #endif