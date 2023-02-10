#include "auxiliarDefinitions.h"
#include "MechanicalProperties.h"
#include "xmlModelWriter.h"
#include "NiftySimEjecutable.h"
#include "auxiliarFunctions_Transformations.h"

#include <unistd.h>

#include "itkImage.h"
#include "itkImageFileReader.h"

#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"

typedef itk::Image<unsigned char,3> ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::Point< double, ImageType::ImageDimension > PointType;


void LabelElements(ImageType::Pointer image, std::vector<double> initial_points, std::vector<int> elements, std::vector<int> & tissues)
{
	double vertex_1[3] = {0.0,0.0,0.0};
	double vertex_2[3] = {0.0,0.0,0.0};
	double vertex_3[3] = {0.0,0.0,0.0};
	double vertex_4[3] = {0.0,0.0,0.0};

	// double barycenter[3] = {0.0,0.0,0.0};

	ImageType::IndexType pixelIndex;
	PointType barycenter;
	bool isInside = true;

	ImageType::PixelType pixelValue;
	int count_good =0;
	int count_bad = 0;

	for(int i=0; i<(elements.size()/4); i++)
	{
		vertex_1[0] = initial_points[ 3*elements[4*i] ];
		vertex_1[1] = initial_points[ 3*elements[4*i]+1 ];
		vertex_1[2] = initial_points[ 3*elements[4*i]+2 ];

		vertex_2[0] = initial_points[ 3*elements[4*i+1] ];
		vertex_2[1] = initial_points[ 3*elements[4*i+1]+1 ];
		vertex_2[2] = initial_points[ 3*elements[4*i+1]+2 ];

		vertex_3[0] = initial_points[ 3*elements[4*i+2] ];
		vertex_3[1] = initial_points[ 3*elements[4*i+2]+1 ];
		vertex_3[2] = initial_points[ 3*elements[4*i+2]+2 ];

		vertex_4[0] = initial_points[ 3*elements[4*i+3] ];
		vertex_4[1] = initial_points[ 3*elements[4*i+3]+1 ];
		vertex_4[2] = initial_points[ 3*elements[4*i+3]+2 ];
	
		barycenter[0] = (vertex_1[0] + vertex_2[0] + vertex_3[0] + vertex_4[0]) /4;
		barycenter[1] = (vertex_1[1] + vertex_2[1] + vertex_3[1] + vertex_4[1]) /4;
		barycenter[2] = (vertex_1[2] + vertex_2[2] + vertex_3[2] + vertex_4[2]) /4;

		isInside = image->TransformPhysicalPointToIndex( barycenter, pixelIndex );

		if ( isInside )
		{
			pixelValue = image->GetPixel( pixelIndex);
			if(pixelValue == 0)	tissues[i] = 3;
			else if(pixelValue >3 ) tissues[i] = 2;
			else tissues[i] = (int) pixelValue;

			count_good +=1;
			
		} else {
			count_bad +=1;
		}
	}
	std::cout << "Good : " << count_good <<  std::endl;
	std::cout << "Bad : " << count_bad << std::endl;
}

int main( int argc, char* argv[])
{
    if(argc<4){
        std::cout << "usage: " << argv[0] << " mesh_filename image_filename ouput_filename thickness" << std::endl;
        return EXIT_FAILURE;
    }

    std::string mesh = argv[1];
	std::string image = argv[2];
    std::string output = argv[3];
	float thick = atof(argv[4]);

	// Itk
	ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( image );
		try{
			reader->Update();
		}catch(itk::ExceptionObject & excp ){
			std::cout << "Reader Exception Object!" << std::endl;
			std::cout << excp << std::endl;
			return EXIT_FAILURE;
		}

    std::cout << "imagen leida!"<< std::endl;

	// Vtk
    vtkSmartPointer< vtkUnstructuredGridReader> meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( mesh.c_str() );
        meshReader->Update();

    std::cout << "malla leida!" << std::endl;

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

	std::cout << "Number of points: " << NoP << std::endl;

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
    
    std::cout << "Number of Cells" << i_cells->GetNumberOfCells()  << std::endl;
    // Tissue!
	//vtkDataArray* tisue = u_grid_reader->GetOutput()->GetCellData()->GetArray(0);
	std::vector<int> tissues; // = input_mri->getTissueElements_model();
	for(int i=0; i<i_cells->GetNumberOfCells(); i++)
	{
		tissues.push_back( 1 );
	}
	LabelElements(reader->GetOutput(), initial_points, elements, tissues);
	
	auto biggest = std::max_element(std::begin( tissues ), std::end( tissues ));
	std::cout << "El valor biggest es de : " << *biggest << std::endl;

	// BoundaryConditions!
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

		double a0 = shearModulus( E_fat, nuPoisson);
			myParameters.material_parameters[0] = a0;
		double a1 = bulkModulus( E_fat, nuPoisson);
			myParameters.material_parameters[1] = a1;
		double a2 = shearModulus( E_gland, nuPoisson);
			myParameters.material_parameters[2] = a2;
		double a3 = bulkModulus( E_gland, nuPoisson );
			myParameters.material_parameters[3] = a3;

		double a4 = shearModulus( 4*E_gland, nuPoisson);
			myParameters.material_parameters[4] = a4;
		double a5 = bulkModulus( 4*E_gland, nuPoisson);
			myParameters.material_parameters[5] = a5;

	std::cout << std::endl;
	std::cout << "Tejido Glandular: " << std::endl;
	std::cout << "\t modulo de Young: E="<< E_gland <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	std::cout << "\t shear modulus: s="<< a2 <<", bulk modulus: m=" << a3 << std::endl;
	std::cout << std::endl;
	std::cout << "Tejido Fatty: " << std::endl;
	std::cout << "\t modulo de Young: E="<< E_fat <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	std::cout << "\t shear modulus: s="<< a0 <<", bulk modulus: m=" << a1 << std::endl;
	std::cout << std::endl;

		// Los parametros se han de actualizar al hacer cada nueva transformaci'on!!!

	NiftySimEjecutable * nifty = new NiftySimEjecutable();
		nifty->SetModelParameters( myParameters	);
		nifty->SetOutputDir( output );
	int sys_i =	nifty->Update();
	if(sys_i==1) return EXIT_FAILURE;


    return EXIT_SUCCESS;
}
