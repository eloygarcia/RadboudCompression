#include "auxiliarDefinitions.h"
#include "MechanicalProperties.h"
#include "xmlModelWriter.h"
//#include "NiftySimEjecutable.h"

#include <unistd.h>

#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"

int main( int argc, char* argv[])
{
    if(argc<3){
        std::cout << "usage: " << argv[0] << " mesh_filename ouput_filename" << std::endl;
        return EXIT_FAILURE;
    }

    std::string mesh = argv[1];
    std::string output = argv[2];

    vtkSmartPointer< vtkUnstructuredGridReader> meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( mesh.c_str() );
        meshReader->Update();

    // =============================== EXTRACCION DE LA MALLA ===========================
	// Points
	std::vector<double> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = meshReader->GetOutput()->GetPoints();
	int NoP = meshReader->GetOutput()->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};
	for (int i=0; i<NoP; i++)
	{
		i_points->GetPoint(i,temp_point);

		initial_points.push_back( temp_point[0] );
		initial_points.push_back( temp_point[1] );
		initial_points.push_back( temp_point[2] );
	}
	
	std::vector<double> final_points;
	for(int xx = 0; xx< initial_points.size(); xx++){	final_points.push_back( 0 ); }
	
	// Elements
	int accumm = 0;
	int minimum = 999999999;
	std::vector<int> elements; // = input_mri->getElements_model();
	vtkCellArray * i_cells = meshReader->GetOutput()->GetCells(); // LA PUTA QUE LOS PARIO !!
	vtkIdList * pts = vtkIdList::New();

    for(int i=0; i<i_cells->GetNumberOfCells(); i++)
	{
		i_cells->GetCell(accumm, pts);
		for (int j=0; j<pts->GetNumberOfIds(); j++)
		{
			elements.push_back( pts->GetId(j) );
			//if( pts->GetId(j) < minimum ) minimum = pts->GetId(j);
		}
		accumm = accumm + 1 + pts->GetNumberOfIds();
	}
    

    // Tissue!
	//vtkDataArray* tisue = u_grid_reader->GetOutput()->GetCellData()->GetArray(0);
	std::vector<int> tissues; // = input_mri->getTissueElements_model();
	for(int i=0; i<i_cells->GetNumberOfCells(); i++)
	{
		tissues.push_back( 0 );
	}
std::cout << "LLega aquí 0" << std::endl;
	// BoundaryConditions!
	vtkDataArray* BoundCond = meshReader->GetOutput()->GetPointData()->GetArray(0);
    std::cout << BoundCond->GetSize() << std::endl;
    std::cout << i_points->GetNumberOfPoints() << std::endl;
    //std::cout << i_cells->GetNumberOfCells() << std::endl;
    //usleep(10000000);
	double tuple = 0;
	std::vector<int> boundaryConditions; // = input_mri->getBoundaryConditions_model();
	std::cout << "Boundary Condition! "  << std::endl;
	for(int i=0; i<i_points->GetNumberOfPoints(); i++)
	{
		tuple = BoundCond->GetTuple1(i);
        //std::cout << BoundCond->GetTuple1(i) << std::endl;
        //usleep(1000);
		boundaryConditions.push_back( (int) tuple );
	}
    //std::cout << "HAsta aqui" << std::endl;
    //usleep(5000);
    /**/
std::cout << "LLega aquí 1" << std::endl;
    xmlModelWriter * model = new xmlModelWriter;

    model->addNodes(initial_points,"3");
    model->addElements(elements,"T4");

    std::vector<double> material_parameters;
    float nuPoisson = 0.495;
	float E_fat = 4.46*1000 ;
	float E_gland = 15.1*1000 ;
    //double a0 = shearModulus( E_fat, nuPoisson);
		material_parameters.push_back(shearModulus((E_fat+E_gland)/2 , nuPoisson ));
	//double a1 = bulkModulus( E_fat, nuPoisson);
		material_parameters.push_back(bulkModulus((E_fat+E_gland)/2 , nuPoisson ));

    std::vector<int> unique_element_vector;
    for( int i=0; i<(elements.size()/4); i++) unique_element_vector.push_back(i); 
    model->addElementSet( unique_element_vector, "NH", material_parameters ); 

std::cout << "LLega aquí" << std::endl;
    // Boundary Conditions
    std::vector<int> bC_list;
	for(int j=0; j<boundaryConditions.size(); j++)
	{
		if(boundaryConditions[j]!=0) bC_list.push_back(j);	
	}
	model->addConstraint(bC_list,"Fix","2");
    /*for(int i=1;i<4; i++)
    {
        bC_list.clear();
        for(int j=0; j<boundaryConditions.size(); j++)
        {
            if(boundaryConditions[j]==i) bC_list.push_back(j);
        }
        if(i<2) model->addConstraint(bC_list, "Fix", "2");
        else model->addConstraint(bC_list, "Fix","2");
    }
	*/
std::cout << "Pasa bien"<< std::endl;
    

    std::vector<int> unique_node;
    for (int i=0; i<i_points->GetNumberOfPoints(); i++) unique_node.push_back(i);
    std::vector<double> accelerationDirection;
        accelerationDirection.push_back(0);
        accelerationDirection.push_back(1);
        accelerationDirection.push_back(0);
    model->addGravity(unique_node, "250", accelerationDirection );

    // =============================== System Params ========================
	model->addSystemParams("<TimeStep>", "0.0001");
	model->addSystemParams("<TotalTime>", "1");
	model->addSystemParams("<DampingCoeff>", "0.069");
	model->addSystemParams("<Density>", "1");
	model->addSystemParams("<DoDeformableCollision>", "0");

	// ======================== Model Writer =========================== 
    std::string niftysim_filename = output ;
    model->writeModel( (char*) niftysim_filename.c_str());


    return EXIT_SUCCESS;
}