#ifndef __niftysimEjecutable_cpp
#define __niftysimEjecutable_cpp

#include "NiftySimEjecutable.h"

NiftySimEjecutable::NiftySimEjecutable(void)
{
	sys_a = 0;
	de = 0.0;

	// ==== De: Compression Function =======
	niftysim_fn = "/NiftySim-model.xml";
	//princ = "\"C:\\niftysim\\source\\Debug\\niftysim.exe\" -x ";

	// ========= De:   ====================
	niftysimOutput_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	UGridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	comprrGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	pt = new double[3];
	d = new double[3];
	temp = new double[3];

	auxPoint = new double[3];

	//new_model = new xmlModelWriter;
}

NiftySimEjecutable::~NiftySimEjecutable(void)
{
	origen.~vector();
	spacing.~vector();
	size.~vector();

	niftysim_fn.erase( niftysim_fn.begin(), niftysim_fn.end() );
	princ.erase( princ.begin(), princ.end() );
	niftysim_filename.erase( niftysim_filename.begin(), niftysim_filename.end() );
	action.erase( action.begin(), action.end() );

	temp_mat.~vector();

	unique_element_vector.~vector();
	temp_vect.~vector();

	bb.~vector();

	m_i_points.~vector();
	compressedgridname.erase( compressedgridname.begin(), compressedgridname.end() );
	vtk_niftysim_filename.erase( vtk_niftysim_filename.begin(), vtk_niftysim_filename.end() );

	// delete[] m_parameters;
	m_outputDir.erase( m_outputDir.begin(), m_outputDir.end() );
	m_finalpoints.~vector();

	delete [] pt;
	delete [] d;
	delete [] temp;
	delete [] auxPoint;

	a.~vector();
	b.~vector();
	c.~vector();
	disp.~vector();

}

// private
void NiftySimEjecutable::CompressionFunction(modelParameters myParameters, std::string output_dir)
{	
	// ================================ NiftySim ==============================
	std::cout << std::endl;
	std::cout << " ======================== NiftySim ======================== " << std::endl;
	std::cout << std::endl;

	//std::cout << "ANTES DE ESCRIBIR EL XML" << std::endl;
	//std::string princ = "\"C:\\niftysim\\source\\Debug\\niftysim.exe\" -x ";
	// =========================================================
	// std::string niftysim_fn = "\\NiftySim-model.xml";
	 
	// niftysim_filename = output_dir + niftysim_fn;
	niftysim_filename = output_dir;
	// =========================================================
	
	action = princ + output_dir + niftysim_fn + " -sport -t -v -export-mesh " +
							output_dir + "\\" + myParameters.position + "_NiftySim-meshgrid.vtk";

	// ======================= NiftySim Write Model ===================	
	xmlModelWriter * new_model = new xmlModelWriter;
	// Nodes.
	new_model->addNodes( myParameters.nodes , myParameters.nodes_Dof);
	// Elements.
	new_model->addElements(myParameters.elements, myParameters.type);
	// Elements Set = 1.
	if(myParameters.number_of_tissues == 1) {
		//std::vector<int> unique_element_vector;
		//std::vector<double> temp_mat;
		unique_element_vector.clear() ; temp_mat.clear();
		for(unsigned int i_1 =0; i_1<(myParameters.elements.size()/4); i_1++) { unique_element_vector.push_back(i_1); };
		if(myParameters.material_parameters.size() > 2) { 
			std::cout << "Found " << myParameters.material_parameters.size() << " material parameteres. Just taking 2 first." << std::endl;
				temp_mat.push_back( myParameters.material_parameters[0] );
				temp_mat.push_back( myParameters.material_parameters[1] );
			new_model->addElementSet( unique_element_vector, myParameters.material, temp_mat ); 
		} else {
			new_model->addElementSet( unique_element_vector, myParameters.material, myParameters.material_parameters ); 
		}
	}else {
		//std::vector<int> temp_vect;
		//std::vector<double> temp_mat;
		temp_vect.clear() ; temp_mat.clear();
		for(int i=0; i<myParameters.number_of_tissues; i++) {
			for(unsigned int j=0; j<myParameters.elementTissue.size(); j++){
				if( myParameters.elementTissue[j]== i+1 ) temp_vect.push_back( j );
			}
			if(i<5){
			temp_mat.push_back( myParameters.material_parameters[ 2*i ] );
			temp_mat.push_back( myParameters.material_parameters[ 2*i+1 ] );
			} else{ 
			temp_mat.push_back( myParameters.material_parameters[ 2 ] ) ;
			temp_mat.push_back( myParameters.material_parameters[ 2+1 ] );
			}

			new_model->addElementSet( temp_vect, myParameters.material, temp_mat); // material 
			// & material_parameters hay que cambiar. 
			// Nop, material no se cambia en este caso, porque ambos son "NH"

			temp_mat.clear();
			temp_vect.clear();
		}
	}//else return EXIT_FAILURE;
	
	// Boundary Conditions.
	// new_model->addConstraint(myParameters.constraint_bC_list, "Fix" , "0"); 
	// new_model->addConstraint(myParameters.constraint_bC_list, "Fix" , "1"); 
	// new_model->addConstraint(myParameters.constraint_bC_list, "Fix" , "2");

	new_model->addConstraint(myParameters.constraint_bC_list, 
		myParameters.constraint_type, 
		myParameters.constraint_axis);
	
	// ====================== Contact Plates ========================
	//std::list<int> unique_node_vector;
	u_node_vector.clear() ;
	for(unsigned int i_1 =0; i_1<(myParameters.nodes.size()/3); i_1++)
	{
		unique_node_vector.push_back(i_1);
	}
	// std::vector<double> bb = myParameters.boundingBox;
	bb.clear() ;  bb = myParameters.boundingBox;
	de = (((float)abs(bb[3] - (bb[2]+20.0) ))-((float) myParameters.thickness));

	std::cout << std::endl;
	std::cout << "Bounding box = [" << bb[0] <<", "<< bb[1] <<", "<< bb[2] <<", ";
	std::cout << bb[3] <<", "<< bb[4] <<", "<< bb[5] <<"] "<< std::endl;

	// 1. Plate 1
	a.clear() ; a.push_back(bb[0]-50);	a.push_back(bb[2]);	a.push_back(bb[4]-50);
	b.clear() ; b.push_back(bb[1]+50);	b.push_back(bb[2]);	b.push_back(bb[4]-50);
	c.clear() ; c.push_back(bb[0]-50);	c.push_back(bb[2]);	c.push_back(bb[5]+50);

	disp.clear() ; disp.push_back(0); disp.push_back(+20); disp.push_back(0);
	std::cout << "Displacement: [" << disp[0] << ", " << disp[1] << ", " << disp[2] << "]" << std::endl;
//	
	new_model->addContactPlate(a,c,b,disp, unique_node_vector);		
//
	// Modelo de mamograf�a !!

	// 2. Plate 2.
	a.clear() ; a.push_back(bb[0]-50);	a.push_back(bb[3]);	a.push_back(bb[4]-50);
	b.clear() ; b.push_back(bb[0]-50);	b.push_back(bb[3]);	b.push_back(bb[5]+50);
	c.clear() ; c.push_back(bb[1]+50);	c.push_back(bb[3]);	c.push_back(bb[4]-50);
//
	disp.clear() ; disp.push_back(0); disp.push_back(-de); disp.push_back(0);
	std::cout << "Displacement: [" << disp[0] << ", " << disp[1] << ", " << disp[2] << "]" << std::endl;
//
	new_model->addContactPlate(a,c,b,disp, unique_node_vector);

	//std::vector<int> node_vect;
	//	for(int i =0; i<unique_node_vector.size(); i++) node_vect.push_back( unique_node_vector[i]);
	std::vector<double> AccelerationDirection;
		AccelerationDirection.push_back( 0 );
		AccelerationDirection.push_back( -1 );
		AccelerationDirection.push_back( 0 );
	new_model->addGravity(unique_node_vector, "100", AccelerationDirection);


	// =============================== System Params ========================
	new_model->addSystemParams("<TimeStep>", myParameters.timeStep);
	new_model->addSystemParams("<TotalTime>", myParameters.totalTime);
	new_model->addSystemParams("<DampingCoeff>", myParameters.DampingCoeff);
	new_model->addSystemParams("<Density>", myParameters.Density);
	new_model->addSystemParams("<DoDeformableCollision>", myParameters.DoDeformableCollision);

	// ======================== Model Writer =========================== 
	new_model->writeModel( (char*) niftysim_filename.c_str());
	std::cout << niftysim_filename << std::endl;

	//std::cout << "ACABA DE SALIR DE ESCRIBIR EL XML" << std::endl;
	// delete new_model[];

	// ======================= Executing NiftySim
/*	std::cout<< std::endl;
	std::cout << action <<std::endl; 
	std::cout<< std::endl;
	
	sys_a = 0;
	sys_a=system( (const char*) action.c_str() );

	delete new_model;

	std::cout<< std::endl;
	std::cout << "NiftySim return value: " << sys_a << std::endl;  // 1= mal acabada la simulacion
	std::cout<< std::endl;

	return sys_a;
	//return 0;
	*/
}

void NiftySimEjecutable::createUnstructuredGrid( std::vector<double> point_vector, vtkCellArray* cells, vtkSmartPointer<vtkUnstructuredGrid> &grid)
{
	int numberofpoints = point_vector.size()/3;
	//int numberofcells = cells_vector.size()/4;

	// ========================================================
	vtkPoints * points = vtkPoints::New();
	// vtkCellArray* cells = vtkCellArray::New();

	points->SetNumberOfPoints( numberofpoints );
	//double auxPoint[3] = {0.0, 0.0, 0.0};
		auxPoint[0] = 0.0;
		auxPoint[1] = 0.0;
		auxPoint[2] = 0.0;
	for(int pointId = 0; pointId < numberofpoints; pointId++)
	{
		auxPoint[0] = point_vector[ pointId*3 ];
		auxPoint[1] = point_vector[ pointId*3 +1];
		auxPoint[2] = point_vector[ pointId*3 +2];

		points->SetPoint(pointId, auxPoint);
	}

	grid->SetPoints(points);
	grid->SetCells(VTK_TETRA, cells);

	points->Delete();
}

void NiftySimEjecutable::getPoints_compressedBreast( 
					vtkSmartPointer<vtkUnstructuredGrid> niftysim_mesh,
					std::vector<double> &f_points)
{
	//vtkPoints* initial_points = niftysim_mesh->GetPoints();
	//vtkDataArray* disp = niftysim_mesh->GetPointData()->GetArray(0);
	initial_points = niftysim_mesh->GetPoints();
	vtk_disp = niftysim_mesh->GetPointData()->GetArray(0);

	//double pt[3] = {0.0, 0.0, 0.0};
	//double d[3] = {0.0, 0.0, 0.0};
	//double temp[3] = {0.0, 0.0, 0.0};
		pt[0] = 0.0; pt[1] = 0.0; pt[2] = 0.0;
		d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;
		temp[0] = 0.0; temp[1] = 0.0; temp[2] = 0.0;

	for(int i=0; i<initial_points->GetNumberOfPoints(); i++)
	{
		initial_points->GetPoint(i, pt);
		vtk_disp->GetTuple(i, d);

		temp[0] = pt[0] + d[0];
		temp[1] = pt[1] + d[1];
		temp[2] = pt[2] + d[2];

		f_points[ 3*i ] = temp[0];
		f_points[ 3*i +1 ] = temp[1];
		f_points[ 3*i +2 ] = temp[2];
	}

	// Cells
	// vtkCellArray* cells = niftysim_mesh->GetCells();
	cells = niftysim_mesh->GetCells();

	// Writen the compressed breast !
	//vtkSmartPointer<vtkUnstructuredGrid> comprrGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	createUnstructuredGrid( f_points, cells, comprrGrid);

	compressedgridname = m_outputDir + "\\" + m_parameters.position + "_compressedGrid.vtk";
	//vtkSmartPointer<vtkUnstructuredGridWriter> UGridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		UGridWriter->SetInputData( comprrGrid );
		UGridWriter->SetFileName( compressedgridname.c_str() );
		UGridWriter->Update();
		UGridWriter->Write();
}

// public
int NiftySimEjecutable::Update()
{
	 // #ifndef _GPU_
	m_i_points.clear() ;	m_i_points = m_parameters.nodes;
	for(unsigned int xx = 0; xx< m_i_points.size(); xx++){	m_finalpoints.push_back( 0 ); }

	//CompressionFunction( m_parameters, m_outputDir );
	
	CompressionFunction( m_parameters , m_outputDir );
/*		if( sys_a!=0) return 1;

	vtkSmartPointer< ErrorObserver > errorObserver = vtkSmartPointer< ErrorObserver >::New();

	vtk_niftysim_filename = m_outputDir + "\\" + m_parameters.position + "_NiftySim-meshgrid.vtk";
	//vtkSmartPointer<vtkUnstructuredGridReader> niftysimOutput_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
		niftysimOutput_reader->SetFileName( (const char*) vtk_niftysim_filename.c_str() );
		niftysimOutput_reader->AddObserver( vtkCommand::ErrorEvent, errorObserver);
		niftysimOutput_reader->AddObserver( vtkCommand::WarningEvent, errorObserver);
		niftysimOutput_reader->Update();
		if( (errorObserver->GetError()) ){
			std::cout << "Vtk Error Event!" << std::endl;
			std::cout << errorObserver->GetErrorMessage() << std::endl;
			return 1;
		}
		if( (errorObserver->GetWarning()) ){
			std::cout << "Vtk Warning Event!" << std::endl;
			std::cout << errorObserver->GetWarningMessage() << std::endl;
			return 1;
		}
	//vtkSmartPointer<vtkUnstructuredGrid> niftygrid = niftysimOutput_reader->GetOutput();
	//	niftygrid = niftysimOutput_reader->GetOutput();

	getPoints_compressedBreast( niftysimOutput_reader->GetOutput(), m_finalpoints);
	// #endif
*/
	return 0;
}

#endif
