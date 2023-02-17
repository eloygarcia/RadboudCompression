#include <iostream>
#include <fstream>

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkStaticCleanUnstructuredGrid.h"

void Usage(char* argv[])
{
    std::cout << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << argv[0] << " <mesh_file.vtk> " << std::endl;
    std::cout << std::endl;
    std::cout << "CUIDADO!!"<< std::endl;
    std::cout << "Esta función sobreescribe la malla original!" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    if(argc<2 || strcmp(argv[1],"--h")==0 || strcmp(argv[1],"--help")==0 ) {
        Usage(argv);
        return EXIT_SUCCESS;
    }

    std::string mesh = argv[1];

    vtkSmartPointer< vtkUnstructuredGridReader> meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( mesh.c_str() );
        meshReader->Update();

    std::cout << std::endl;
    std::cout << "Initial number of points: "<< meshReader->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "Initial number of cells: "<< meshReader->GetOutput()->GetNumberOfCells() << std::endl;

    vtkSmartPointer< vtkStaticCleanUnstructuredGrid > cleaningMesh = vtkSmartPointer< vtkStaticCleanUnstructuredGrid>::New();
        cleaningMesh->SetInputData(meshReader->GetOutput());
        cleaningMesh->ToleranceIsAbsoluteOn();
        cleaningMesh->SetTolerance(0.0);
        cleaningMesh->RemoveUnusedPointsOn();
        cleaningMesh->Update();

    std::cout << std::endl;
    std::cout << "Final number of points: "<< cleaningMesh->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "Final number of cells: "<< cleaningMesh->GetOutput()->GetNumberOfCells() << std::endl;

    vtkSmartPointer< vtkUnstructuredGridWriter > writer = vtkSmartPointer< vtkUnstructuredGridWriter>::New();
        writer->SetFileName( mesh.c_str() );
        writer->SetInputData( cleaningMesh->GetOutput());
        writer->Update();

    return EXIT_SUCCESS;
}