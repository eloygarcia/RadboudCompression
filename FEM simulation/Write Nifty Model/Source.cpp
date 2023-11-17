#include "auxiliarDefinitions.h"
#include "MechanicalProperties.h"
#include "xmlModelWriter.h"
#include "NiftySimEjecutable.h"
#include "auxiliarFunctions_Transformations.h"

#include <unistd.h>

#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkCellData.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"



int main( int argc, char* argv[])
{
    if(argc<3){
        std::cout << "usage: " << argv[0] << " mesh_filename ouput_filename thickness (in mm.)" << std::endl;
        return EXIT_FAILURE;
    }

    std::string mesh = argv[1];
    std::string output = argv[2];
	float thick = atof(argv[3]);

	// Vtk
    vtkSmartPointer< vtkUnstructuredGridReader> meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( mesh.c_str() );
        meshReader->Update();

    std::cout << "malla leida!" << std::endl;

    // =============================== EXTRACCION DE LA MALLA ===========================

	// Extracting Points
	std::vector<double> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = meshReader->GetOutput()->GetPoints();
	int NoP = meshReader->GetOutput()->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};
	/* TO DO */
	/* Tal vez podamos cambiar estas líneas para poder hacer de un modo mśa elegante la transformación
	/*   vtkPoints a std::vector
	/* */
	for (int i=0; i<NoP; i++)
	{
		i_points->GetPoint(i,temp_point);

		initial_points.push_back( temp_point[0] );
		initial_points.push_back( temp_point[1] );
		initial_points.push_back( temp_point[2] );
	}

	std::cout << "Number of points: " << NoP << std::endl;

	// Extracting Elements
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
    
    std::cout << "Number of Cells" << i_cells->GetNumberOfCells()  << std::endl;

    // Extracting tissue Tissue!
	vtkDataArray* materialtissues = meshReader->GetOutput()->GetCellData()->GetArray(0); // creo que es exta línea.
	std::vector<int> tissues; // = input_mri->getTissueElements_model();

	std::cout << "CHECK THIS BEFORE CONTINUING" << std::endl;
	std::cout << "Number of cells: " << i_cells->GetNumberOfCells() << std::endl;
	std::cout << "Name of cell data: " << materialtissues->GetName() << std::endl;
	std::cout << "Number of cell data tuples: " << materialtissues->GetNumberOfTuples() << std::endl;

	for(int i=0; i<i_cells->GetNumberOfCells(); i++)
	{
		tissues.push_back( materialtissues->GetTuple1(i) );
	}
	//LabelElements(reader->GetOutput(), initial_points, elements, tissues);
		auto biggest = std::max_element(std::begin( tissues ), std::end( tissues ));
	std::cout << "El valor biggest es de : " << *biggest << std::endl;



	// Extracting nodes BoundaryConditions!
	vtkDataArray* BoundCond = meshReader->GetOutput()->GetPointData()->GetArray(0);
	double tuple = 0;
	std::vector<int> boundaryConditions; // = input_mri->getBoundaryConditions_model();

	//std::cout << "Boundary Condition! "  << std::endl;
	//std::cout << "Number of points: " << i_points->GetNumberOfPoints() << std::endl;

    std::cout << "if the software fails here, please, check the boundary conditions" << std::endl;
	for(int i=0; i<i_points->GetNumberOfPoints() ; i++)
	{
		tuple = BoundCond->GetTuple1(i);  // Falta un try-catch si no hay boundary conditions
		if(tuple==1)
		{
			boundaryConditions.push_back( (int)i );
		}
	}


// =========================== Model Parameters. ============================
    std::cout << "Entra Model parameters" << std::endl;
	modelParameters myParameters;
		myParameters.name = 'phantom';

		myParameters.number_of_tissues =*biggest;
		//myParameters.number_of_tissues = 1;

		std::vector<float> cen;
			cen.push_back(0);
			cen.push_back(0);
		myParameters.cen = cen; // centroid 
		myParameters.side = 'L'; //temp_side; 
		myParameters.angle = 0; //angle; 

		myParameters.nodes =  initial_points;
		myParameters.boundingBox = get_boundingBox( initial_points ); // intensity-based registration
		// myParameters.boundingBox = get_boundingBox( final_points ); // intensity-based registration
		myParameters.elements = elements;

			 	myParameters.subsets_on = 1;
				myParameters.number_of_tissues = *biggest;
				myParameters.elementTissue = tissues;

		myParameters.constraint_bC_list = boundaryConditions;
		myParameters.constraint_type = "Fix";
		myParameters.constraint_axis = "0";

		myParameters.position = "CC"; // position;
		myParameters.thickness = thick;
		myParameters.pectoralAngle = 0; // input_mammogram->getPectoralAngle();

		myParameters.paddleAngle = 0;

	float nuPoisson = 0.499;
	float E_fat = 4.46*1000 ;
	float E_gland = 15.1*1000 ;

        // 1. Fatty Tissue
		double a0 = shearModulus( E_fat, nuPoisson);
			myParameters.material_parameters[0] = a0;
		double a1 = bulkModulus( E_fat, nuPoisson);
			myParameters.material_parameters[1] = a1;

		// 2. Glandular Tissue
		double a2 = shearModulus( E_gland, nuPoisson);
			myParameters.material_parameters[2] = a2;
		double a3 = bulkModulus( E_gland, nuPoisson );
			myParameters.material_parameters[3] = a3;

        // 3. Skin Tissue
		double a4 = shearModulus( 4*E_gland, nuPoisson);
			myParameters.material_parameters[4] = a4;
		double a5 = bulkModulus( 4*E_gland, nuPoisson);
			myParameters.material_parameters[5] = a5;

        // 4. Pectoral Muscle tissue
		double a6 = shearModulus( 8*E_gland, nuPoisson);
			myParameters.material_parameters[6] = a6;
		double a7 = bulkModulus( 8*E_gland, nuPoisson);
			myParameters.material_parameters[7] = a7;


	std::cout << std::endl;
	std::cout << "Tejido Glandular: " << std::endl;
	std::cout << "\t modulo de Young: E="<< E_gland <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	std::cout << "\t shear modulus: s="<< a2 <<", bulk modulus: m=" << a3 << std::endl;
	std::cout << std::endl;
	std::cout << "Tejido Fatty: " << std::endl;
	std::cout << "\t modulo de Young: E="<< E_fat <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	std::cout << "\t shear modulus: s="<< a0 <<", bulk modulus: m=" << a1 << std::endl;
	std::cout << std::endl;

	NiftySimEjecutable * nifty = new NiftySimEjecutable();
		nifty->SetModelParameters( myParameters	);
		nifty->SetOutputDir( output );
	int sys_i =	nifty->Update();
	if(sys_i==1) return EXIT_FAILURE;

    std::cout << "Exit success" << std::endl;
    return EXIT_SUCCESS;
}
