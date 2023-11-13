// Patrones de diseño, E. GAmma et al. Pearson Education
// C.1 Lista, página 340
#include <iostream>
#include <stdio.h>

#define CAPACIDAD_LISTA_DETERMINADA 512

template <class T>
class Lista
{
    public:
        Lista(long length = CAPACIDAD_LISTA_DETERMINADA){};
        Lista(Lista&){}; // Redefine el constructor de copia predeterminado para que los datos miembros se inicialicen adecuadamente.
        ~Lista(){};
        Lista& operator=(const Lista&); // Implementa la operación de asignacioón para asignar correctamente los datos miembros.

        long Contar() const;
        T& Obtener(long indice) const;
        T& Primero() const;
        T& Ultimo() const;
        bool Incluye(const T&) const;

        void Insertar(const T&) const{};
        void InsertarAlPrincipio(const T&){};

        void Eliminar(const T&) {};
        void EliminarUltimo() {};
        void EliminarPrimero() {};
        void EliminarTodos() {};

        //T& Cima() const {return T&;};
        //void Meter(const T&) ;
        //T& Sacar();
};

template<class T>
class Iterador
{
    public:
        virtual void Primero() = 0;
        virtual void Siguiente() = 0;
        virtual bool HaTerminado() const=0;
        virtual T ElementoActual() const=0;
    protected:
        Iterador();
};

template<class T>
class IteradorLista : public Iterador<T>
{
    public:
        IteradorLista(const Lista<T>* lista);
        virtual void Primero();
        virtual void Siguiente();
        virtual bool HaTerminado() const;
        virtual T ElementoActual() const;
};

/*
class Punto
{
    public:
        static const Punto Cero;
        Punto( float* {0.0, 0.0, 0.0} );

        friend Punto operator+(const Punto&, const Punto&);
};
*/

int main()
{
    Lista<int> lista;
    lista.Insertar(2);
    std::cout << lista[0] << std::endl;

    //Punto* pt = new Punto();
    //Punto* pt2 = new Punto({1.0,1.0,1.0});

    //std::cout << pt[0] << std::endl;
    return 0;
}