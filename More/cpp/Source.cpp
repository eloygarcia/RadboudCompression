#include "auxiliarHeaders.h"

// usage
void aval_metrics(){
	std::cout << "Available metrics:" << std::endl;
	std::cout << "\t 0-MI;\n\t 1-NCC;\n\t 2-MSQ;\n\t 3-Dice;\n\t 4-GC;\n\t 5-PI;\n\t 6-GD;\n\t 7-LC" << std::endl;
	std::cout << std::endl;
}

void Usage(char* argv[])
{
	std::cout << " " << std::endl;
	std::cout << "MRI-mammo Registration Project " << std::endl;
	std::cout << argv[0] <<" mammogram_fn densityMap_fn MRIsegmentation_fn MRImask_fn newSegmentation_fn outputDir metric" << std::endl;
	std::cout << std::endl;
	
	aval_metrics();
}

// ========================== MAIN ==============================
int main( int argc, char* argv[])
{
	std::cout << std::endl;
// ===================== Usage =========================
//	if(argc!=5 || argc!=6)
//	{
//		Usage(argv);
//		return EXIT_FAILURE;
//	}

//======================================================
	std::cout << std::endl;
	std::cout << "Mammogram file name:  \t" << argv[1] << std::endl;
	std::cout << std::endl;
	std::cout << "Density Map file name:\t" << argv[2] << std::endl;
	std::cout << std::endl;
	std::cout << "MRI file name:  \t" << argv[3] << std::endl;
	std::cout << std::endl;

// ================== Read Data ========================
	// mammogram
	std::string inputMammogramfilename = argv[1];
	std::string inputDensityMapfilename = argv[2];
	// mri
	std::string inputSegmentationfilename = argv[3];

	std::string inputGridfilename = argv[4];
	// New Segmentation.
	// std::string newSegmentation_filename = argv[5];

	std::string output_dir = argv[5];

	// int metric = atoi( argv[6] );
	//int metric = 0; // i.e. Mutual Information
	 int metric = 1; // i.e. Normalized Cross-Correlation
	 // 0-MI; 1-NCC; 2-MSQ; 3-Dice; 4-GC; 5-PI; 6-GD; 7-LC
	if(metric>7)
	{
		aval_metrics();
		return EXIT_FAILURE;
	}

// ================= Mammogram actions. =================
	mammogram * input_mammogram = new mammogram;
		input_mammogram->read_mammogram_with_metadata( inputMammogramfilename );
		input_mammogram->otsu_segmentation();

	std::cout << std::endl;
	std::cout << "Mammogram Information: " << std::endl;
		// ==================== Mammmogram information =====================
		std::string temp_side = input_mammogram->getImageLaterality(true);
		if ((temp_side!="R ") && (temp_side!="L ")) {
			std::cout << "Check side breast !" << std::endl; 
			return EXIT_FAILURE; 
		}
		input_mammogram->getAnodeTargetMaterial(true);
		input_mammogram->getFilterMaterial(true);
		int KVP = input_mammogram->getKVP(true);
		int source2detector = input_mammogram->getDistanceSourceToDetector(true);
		float focalSpot = input_mammogram->getFocalSpot(true);
		std::string position = input_mammogram->getViewPossition(true);
			if((position.c_str()[0] =='C') & (position.c_str()[1]=='C')) {position.clear(); position ="CC";}
			else if((position.c_str()[0] =='M') & (position.c_str()[1]=='L') & (position.c_str()[2]=='O')) {position.clear(); position="MLO";}
		int thick = input_mammogram->getBodyPartThickness(true);
		input_mammogram->getCompressionForce(true);
		int angle = input_mammogram->getPositionerPrimaryAngle(true);
		std::vector<float> cen = input_mammogram->getCentroid(true);
		if((position.c_str()[0] =='C') & (position.c_str()[1]=='C'))
		{
			if (temp_side=="L ")	cen[0] = cen[0]-20;
			else if (temp_side=="R ")	cen[0] = cen[0]+20;
		} /*else if((position.c_str()[0] =='M') & (position.c_str()[1]=='L') & (position.c_str()[2]=='O')) {
			if (temp_side=="L ")	cen[0] = cen[0]-10;
			else if (temp_side=="R ")	cen[0] = cen[0]+10;
			cen[1] = cen[1]+10;
		}*/

		input_mammogram->read_densityMap( inputDensityMapfilename );
		if((position.c_str()[0] =='M') & (position.c_str()[1]=='L') & (position.c_str()[2]=='O')) {
			input_mammogram->read_pectoralmuscleMask( argv[6] );  // Recordar cambiar si se cambia argc!!
		}

		std::string outputmammogram = output_dir + "/" + position + "_mammogram.mhd";
		input_mammogram->write_mammogram(outputmammogram);
		std::string outputdensitymap = output_dir + "/"+ position + "_densitymap.mhd";
		input_mammogram->write_densitymap(outputdensitymap);

// =================== MRI actions ======================
	// Reader3floatType::Pointer segmentation_reader = Reader3floatType::New();
	ReaderGradientType::Pointer segmentation_reader = ReaderGradientType::New();
		segmentation_reader->SetFileName( inputSegmentationfilename );
		try{
			segmentation_reader->Update();
		} catch (itk::ExceptionObject & excp ) {
			std::cout << "MRI segmentation file Reader Exception!! " << std::endl;
			std::cout << excp << std::endl;
			std::cout << std::endl;
			return EXIT_FAILURE;
		}

	vtkSmartPointer<vtkUnstructuredGridReader> u_grid_reader_A = vtkSmartPointer<vtkUnstructuredGridReader>::New();
		u_grid_reader_A->SetFileName( inputGridfilename.c_str() );
		u_grid_reader_A->Update();

	std::string initial_grid_name = output_dir + "\\InitialUnstructuredGrid.vtk";
	vtkSmartPointer<vtkUnstructuredGridWriter> u_grid_writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		u_grid_writer->SetInputData( u_grid_reader_A->GetOutput() );
		u_grid_writer->SetFileName( initial_grid_name.c_str() );
		u_grid_writer->Update();

	//=========================== INICIALIZACION DE LA TRANSFORMATION ============================
	// EXTAMOS AQUI !!
	std::vector<double> translation;
		translation.push_back(0);
		translation.push_back(0);
		translation.push_back(0);
	
	//float RotX = M_PI_2 ; // M_PI_2 -- Rotaci'on para ponerlo en orden. ;
	float RotX = 0.;
	float RotY = 0.; // 
	float RotZ = 0.;

//	float angle_rad = grad2rad( angle );

/*	if (temp_side=="R ") {
		RotZ = M_PI_2; // Right M_PI_2;
	}else if (temp_side=="L "){
		RotZ = -M_PI_2; // Left -M_PI_2;
	} else {
		return EXIT_FAILURE;
	}
*/
// =============================== EXTRACCION DE LA MALLA ===========================
	/* Feature-based registration */  
	// Points
	std::vector<double> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = u_grid_reader_A->GetOutput()->GetPoints();
	int NoP = u_grid_reader_A->GetOutput()->GetNumberOfPoints();
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
	vtkCellArray * i_cells = u_grid_reader_A->GetOutput()->GetCells(); // LA PUTA QUE LOS PARIO !!
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
	vtkDataArray* tisue = u_grid_reader_A->GetOutput()->GetCellData()->GetArray(0);
	std::vector<int> tissues; // = input_mri->getTissueElements_model();
	double tis=0;
	for( int i=0; i<tisue->GetSize(); i++)
	{
		//tis = tisue->GetTuple1(i); // varios materiales
		//tissues.push_back( tis ); // varios materiales
		tissues.push_back( 1 ); // material unico
	}

	// BoundaryConditions!
	vtkDataArray* BoundCond = u_grid_reader_A->GetOutput()->GetPointData()->GetArray(0);
	double tuple = 0;
	std::vector<int> boundaryConditions; // = input_mri->getBoundaryConditions_model();
	std::cout << "Boundary Condition! "  << std::endl;
	for(int i=0; i<i_points->GetNumberOfPoints(); i++)
	{
		tuple = BoundCond->GetTuple1(i);
		if(tuple==1)
		{
			boundaryConditions.push_back( (int)i );
			//std::cout << i << " " ;
		}
	}

// =========================== Model Parameters. ============================

	auto biggest = std::max_element(std::begin( tissues ), std::end( tissues ));
	std::cout << "El valor biggest es de : " << *biggest << std::endl;
	
	modelParameters myParameters;
		myParameters.number_of_tissues =*biggest;

		myParameters.cen = cen;
		myParameters.side = temp_side;
		myParameters.angle = angle;

		myParameters.nodes = initial_points;
		myParameters.boundingBox = get_boundingBox( final_points ); // intensity-based registration
		myParameters.elements = elements;
			 if(*biggest >=2){
			 	myParameters.subsets_on = 1;
				myParameters.number_of_tissues = *biggest;
				myParameters.elementTissue = tissues;
			 } else {
				myParameters.subsets_on = 0;
				myParameters.number_of_tissues = *biggest;
				myParameters.elementTissue = tissues; 
			 }

		myParameters.constraint_bC_list = boundaryConditions;
		myParameters.position = position;
		myParameters.thickness = (float)thick;
		myParameters.pectoralAngle = input_mammogram->getPectoralAngle();

	float nuPoisson = 0.495;
	float E_fat = 4.46*1000 ;
	float E_gland = 15.1*1000 ;

		double a0 = shearModulus( E_fat, nuPoisson);
			myParameters.material_parameters[0] = a0;
		double a1 = bulkModulus( E_fat, nuPoisson);
			myParameters.material_parameters[1] = a1;
//		double a2 = shearModulus( E_gland, nuPoisson);
//			myParameters.material_parameters[2] = a2;
//		double a3 = bulkModulus( E_gland, nuPoisson );
//			myParameters.material_parameters[3] = a3;

		double a2 = anisotropyTerm( E_fat, 1);
			myParameters.material_parameters[2] = a2;

	//std::cout << std::endl;
	//std::cout << "Tejido Glandular: " << std::endl;
	//std::cout << "\t modulo de Young: E="<< E_gland <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	//std::cout << "\t shear modulus: s="<< a2 <<", bulk modulus: m=" << a3 << std::endl;
	//std::cout << std::endl;
	//std::cout << "Tejido Fatty: " << std::endl;
	//std::cout << "\t modulo de Young: E="<< E_fat <<", ratio de Poisson: nu=" << nuPoisson << std::endl;
	//std::cout << "\t shear modulus: s="<< a0 <<", bulk modulus: m=" << a1 << std::endl;
	//std::cout << std::endl;

std::cout << std::endl;
std::cout << "Anisotropic hyperelastic neo-hookean material model: " << std::endl;
std::cout << "\t modulo de Young: E="<< E_fat << std::endl;
std::cout <<"\t ratio de Poisson: nu=" << nuPoisson << std::endl; 
std::cout << "\t Anisotropic behaviour: eta=" << a2 << std::endl;
std::cout << std::endl;


		// Los parametros se han de actualizar al hacer cada nueva transformaci'on!!!
	
// ===================== TARGET IMAGE ! ===================
/*
	Image2floatType::Pointer target_image =  input_mammogram->get_dense_tissue(true)  ; // Intensity-based registration
	std::string target_name = output_dir + "\\" + position + "_target_image.nrrd";
	Writer2floatType::Pointer writer = Writer2floatType::New();
		writer->SetInput( target_image );
		writer->SetFileName( target_name );
		writer->SetUseCompression( true );
		try{
			writer->Update();
		} catch(itk::ExceptionObject & err) {
			std::cout << std::endl;
			std::cout << "Target Image Writer exception! " << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}
*/
	std::string target_name = output_dir + "\\" + position + "_target_image.nrrd";
	input_mammogram->write_densitymap(target_name);

// ====================== Next Level !  ==============================
	Image2floatType::Pointer target_image =  input_mammogram->get_dense_tissue(true)  ; // Intensity-based registration

	//IntensityBased * ib = new IntensityBased;
	//	ib->SetTargetImage( target_image );
	//	//ib->SetSourceImage( input_mri->getImage_NewSegmentation() );
	//	ib->SetSourceImage( segmentation_reader->GetOutput() );
	//	ib->SetModelParameters( myParameters);
	//	ib->SetOutputDirectory( output_dir );
	//	
	//	ib->SetMetric( metric ); // 0-MI; 1-NCC; 2-MSQ; 3-Dice; 4-GC; 5-PI; 6-GD; 7-LC
	//			// NCC - valores negativos, un mayor valor negativo corresponde con valores más negativos.
	//			// NCC -> minimize !


	GradientBased * ib = new GradientBased;
		ib->SetTargetImage( target_image );
		//ib->SetSourceImage( input_mri->getImage_NewSegmentation() );
		ib->SetSourceImage( segmentation_reader->GetOutput() );
		ib->SetModelParameters( myParameters);
		ib->SetOutputDirectory( output_dir );
		
		ib->SetMetric( metric ); // 0-MI; 1-NCC; 2-MSQ; 3-Dice; 4-GC; 5-PI; 6-GD; 7-LC
				// NCC - valores negativos, un mayor valor negativo corresponde con valores más negativos.
				// NCC -> minimize !



	/* Inicialización de variables !*/
/**/	
	int numberofparameters = 9;									// Anisotropy ratio = 4 //
	float currentPoint[9] = {translation[0], translation[1], RotZ, RotX, RotY, 1.25*E_fat, 0, nuPoisson, (float) thick};
		// E_gland/E_fat, nuPoisson, (float)thick}; // --> {  translationX, translationY, RotationZ, RotationX, RotationY, ShearFatty, BulkFatty, ShearGlandular, BulkGlandular, thickness} 
	// float stepSize[5] = {5.,.5, 0.1, 0.1, 0.5}; // intensity-based 
	// float acceleration = 1.1;
	float epsilon = 0.0;
		if((metric==0) || (metric==1) || (metric==4) || (metric==7))		epsilon = .000001;
		else if ( (metric==2) || (metric==5) || (metric==6) )		epsilon = 1.;
		else if ( (metric==3) )		epsilon = 0.000001;
	
	int n_iter = 50;

		ib->SetTransformation( currentPoint );

/* SIMULATED ANNEALING - uncompleted*/
	float searchSpace[18] = {-75, 75,
								-50, 50, 
								- M_PI*.2, M_PI*.2,
								- M_PI*.15, M_PI*.15,
								- M_PI*.15, M_PI*.15,
								E_fat, E_gland, //1000,7500,
								0,0.000001, //1,10, antes (0,2)
								nuPoisson,1.000001*nuPoisson, // antes(0.495,0.499)
								.75*(float)thick, 1.25*(float)thick}; // antes(0.5*thick, 1.5*thick)
				// Recuerda {RotZ, RotX, RotY}
	// Creo que sería más cómo utilizar el modulo de young y el ratio de poisson con el fin de que no se realizaran simulaciones no realistas!!

	SimulatedAnnealing * simulatedAnnealing = new SimulatedAnnealing();
		simulatedAnnealing->setNumberOfParameters( numberofparameters );
		simulatedAnnealing->setInitialPoint( currentPoint );
		simulatedAnnealing->defineSearchSpace( searchSpace );
		
		// SimulatedAnnealing->setTemperature( 2 );
		simulatedAnnealing->setTemperature( 2000000 );

		simulatedAnnealing->setEpsilon( epsilon );
		simulatedAnnealing->setMaxIter( n_iter );

		// 0-MI; 1-NCC; 2-MSQ; 3-Dice; 4-GC; 5-PI; 6-GD; 7-LC
			if( (metric==0) || (metric==1) || (metric==2) || (metric==4) || (metric==5) || (metric==6) )		simulatedAnnealing->setMinimizeOn();
			else if ( (metric==3) || (metric==7) )		simulatedAnnealing->setMaximizeOn();

		simulatedAnnealing->setOptimFunction( ib );
		simulatedAnnealing->setOutputDir( output_dir );

		simulatedAnnealing->setPosition( position );

int sys_a = simulatedAnnealing->Update();
	if(sys_a!=0) {
		std::cout << std::endl;
		std::cout << "El proceso de convergencia ha fallado!!  " << std::endl;
		return EXIT_FAILURE; 
	} ;



// Extraer la imagen !
/*
	std::cout << std::endl;
	std::cout << " ESTA EXTRAYENDO LA IMAGEN, CON LO CUAL, HA LLEGADO A UNA CONVERGENCIA !!! " << std::endl;
	std::cout << std::endl;

	std::string image_name = output_dir + "\\" + position + "_final_image.nrrd";
	Writer2floatType::Pointer writer_image_final = Writer2floatType::New();
		writer_image_final->SetInput( ib->GetOutputImage()  );
		writer_image_final->SetFileName( image_name );
		writer_image_final->SetUseCompression( true );
		try{
			writer_image_final->Update();
		} catch(itk::ExceptionObject & err) {
			std::cout << std::endl;
			std::cout << "Ouptut Image Writer exception! " << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}

	std::cout << std::endl;
	std::cout << " EXIT SUCCESS ! " << std::endl;
	std::cout << std::endl;
*/
	return EXIT_SUCCESS;
}