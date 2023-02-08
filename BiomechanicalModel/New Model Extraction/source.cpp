#include <iostream>
#include <fstream>

// Itk
#include "itkImage.h"
#include "itkImageFileReader.h"

typedef itk::Image<char,3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;

// Itk-Vtk Glue
#include "itkImageToVTKImageFilter.h"
typedef itk::ImageToVTKImageFilter<ImageType> ItktoVtkType;

// Vtk
#include "vtkImageData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"

#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"

#include "vtkTetra.h"
	
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/read_vtk_image_data.h>

#include <boost/lexical_cast.hpp>
#include <boost/functional.hpp>

typedef short Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Mesh Extraction
typedef CGAL::Triangulation_cell_base_3<Tr>	Cell_Base;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>::Vertices_in_complex_iterator Complex_Vertex_Iterator;

typedef typename C3t3::Cell_iterator Cell_iterator;
typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
typedef typename Tr::Vertex_handle Vertex_handle;
typedef typename Tr::Point Point_3;

void usage()
{
    std::cout << std::endl;
    std::cout << "inputfilename outputfilename arguments" << std::endl;
    std::cout << std::endl;
    std::cout << "- inputfilename : Name of the original volume (3D image) in itk format -i.e. .nrrd or .mhd" << std::endl;
    std::cout << "- outputfilename :  Name of the mesh file with extesion .inp or .vtk"  << std::endl;
    std::cout << "- arguments :" << std::endl;
    std::cout << " \t --value_outside : (default 0)" << std::endl;
    std::cout << " \t --facet_angle : (default 30)" << std::endl;
    std::cout << " \t --facet_size : (default 5)" << std::endl;
    std::cout << " \t --facet_distance : (default 1)" << std::endl;
    std::cout << " \t --cell_radius_edge_ratio : (default 1)" << std::endl;
    std::cout << " \t --cell_size : (default 5)" << std::endl;
};

struct argum{
    argum() {
        value_outside = 0;
        facet_angle = 30;
        facet_size = 5;
        facet_distance = 1;
        cell_radius_edge_ratio = 1;
        cell_size = 5;
    };
    ~argum(){};

    int value_outside;
    double facet_angle;
    double facet_size;
    double facet_distance;
    double cell_radius_edge_ratio;
    double cell_size;
};

argum argParser(int argc, char* argv[])
{
    argum argumentos;
    for(int i=3; i<argc; i++)
    {
        if( strcmp(argv[i],"--value_outside")==0 ) argumentos.value_outside = atof(argv[i+1]);
        if( strcmp(argv[i],"--facet_angle")==0 ) argumentos.facet_angle = atof(argv[i+1]);
        if( strcmp(argv[i],"--facet_size")==0 ) argumentos.facet_size = atof(argv[i+1]);
        if( strcmp(argv[i],"--facet_distance")==0 ) argumentos.facet_distance = atof(argv[i+1]);
        if( strcmp(argv[i],"--cell_radius_edge_ratio")==0 ) argumentos.cell_radius_edge_ratio = atof(argv[i+1]);
        if( strcmp(argv[i],"--cell_size")==0 ) argumentos.cell_size = atof(argv[i+1]);
    }
    return argumentos;
}


void writeInpFile(std::string outputfilename, C3t3 c3t3)
{
    std::ofstream file;
    file.open( outputfilename.c_str() );

    file << "*HEADING \n";
    file << "E.Garcia - O.Diaz Universitat de Girona, 2018 \n";
    file << "Units: Millimetres (mm) \n";
    file << "**============================================================================== \n";

    // Material definitions
    file << "**MATERIAL DEFINITIONS BEGIN \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    file << "**MATERIAL PLACEHOLDER BEGIN \n";
    file << "**Part: MASK 1 \n";
    file << "**Part Material: MASK 1 \n";
    file << "*SOLID SECTION, ELSET=PT_MASK_1, MATERIAL=PM_MASK_1 \n";
    file << "*MATERIAL, NAME=PM_MASK_1 \n";
    file << "**MATERIAL PLACEHOLDER END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    file << "**MATERIAL DEFINITIONS END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    // Contact definitions
    file << "**CONTACT DEFINITIONS BEGIN \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    file << "**CONTACT DEFINITIONS END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    
    // Node data
    file << "**NODE DATA BEGIN \n";
    file << "**============================================================================== \n";
    file << "*NODE \n";
    file << "**------------------------------------------------------------------------------ \n";
    file << "**INDEX, X COORD, Y COORD, Z COORD \n";
    file << "**------------------------------------------------------------------------------ \n";

    Tr tr = c3t3.triangulation();
    int i = 0;
    /*
    // the actual trick: create a map of Point_3 (Vertex_Handle map did not work for me...)
    // std::map<Point_3, int> V;
    for (Tr::All_vertices_iterator it=t.all_vertices_begin(); it != t.all_vertices_end(); ++it)
    {
	    // add the current Point_3 to the map with its current index
	    // V[it->point()] = i;
        if( i>0 ){	
            file << i <<", "<<it->point()[0]<<", " <<it->point()[1]<<", " <<it->point()[2]<<"\n";
            //file << i <<", "<<it->point()<<"\n";
        };
    	++i;
    }
    */

    std::map<Point_3, int> V;
    for (Tr::All_vertices_iterator it=tr.all_vertices_begin(); it != tr.all_vertices_end(); ++it)
    {   
        if( i>0 ){	
            file << i <<", "<<it->point()[0]<<", " <<it->point()[1]<<", " <<it->point()[2]<<"\n";

            V[it->point()] = i;
        };
	    ++i;
    }


    file << "**============================================================================== \n";
    file << "**NODE DATA END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    // Element data
    file << "**SOLID ELEMENT DATA BEGIN \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    
    file << "**ELEMENTS - Part: Mask 1 BEGIN \n";
    file << "**Part Material: MASK 1 \n";
    file << "**Material ID: PM_MASK_1 \n";
    file << "**Material Index: 1 \n";
    
    file << "**============================================================================== \n";
    file << "**ELEMENTS (TETRAHEDRA) - Part: Mask 1 BEGIN \n";
    file << "*ELEMENT, TYPE=C3D4, ELSET=PT_MASK_1 \n";

    i=1;
    /*
    for (Tr::All_cells_iterator it = t.all_cells_begin(); it != t.all_cells_end(); ++it)
    {
        if(i>0){
    	//std::cout<<"Complex #"<<i<<" : "<<it->subdomain_index()<<std::endl;
	    const Tr::Cell c(*it);
    	const Tr::Vertex_handle v0 = c.vertex(0);
    	const Tr::Vertex_handle v1 = c.vertex(1);
    	const Tr::Vertex_handle v2 = c.vertex(2);
    	const Tr::Vertex_handle v3 = c.vertex(3);
    	//std::cout<<v0->point()<<" ; "<<v1->point()<<" ; "<<v2->point()<<" ; "<<v3->point()<<std::endl;
    	//std::cout<<"Indices: "<<std::distance(t.all_vertices_begin(), v0)<<"/"<<std::distance(t.all_vertices_begin(), v1)<<"/"<<std::distance(t.all_vertices_begin(), v2)<<"/"<<std::distance(t.all_vertices_begin(), v3)<<std::endl;

        file << i << "," << std::distance(t.all_vertices_begin(), v0)+1 << "," << std::distance(t.all_vertices_begin(), v1)+1 << "," << std::distance(t.all_vertices_begin(), v2)+1 << "," << std::distance(t.all_vertices_begin(), v3)+1<<"\n";
        }
	    ++i;
    }
    */	
    C3t3::Vertex_handle v0;
	C3t3::Vertex_handle v1;
	C3t3::Vertex_handle v2;
	C3t3::Vertex_handle v3;

    for (C3t3::Cell_iterator it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it)
    {
	    const Tr::Cell c(*it);
	    v0 = c.vertex(0);
	    v1 = c.vertex(1);
	    v2 = c.vertex(2);
	    v3 = c.vertex(3);
	    // std::cout<<"Point #"<<i<<" : "<<v0->point()<<" ; "<<v1->point()<<" ; "<<v2->point()<<" ; "<<v3->point()<<std::endl;
	    // std::cout<<"Indices: "<<V[v0->point()]<<"/"<<V[v1->point()]<<"/"<<V[v2->point()]<<"/"<<V[v3->point()]<<std::endl;

        file << i << "," << V[v0->point()] << "," << V[v1->point()] << "," << V[v2->point()] << "," << V[v3->point()] <<"\n";
	    ++i;
    }

    file << "**ELEMENTS (TETRAHEDRA) - Part: Mask 1 END \n";
    file << "**============================================================================== \n";
    file << "**ELEMENTS - Part: Mask 1 END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    file << "**SOLID ELEMENT DATA END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    // Shell Element data
    file << "**SHELL ELEMENT DATA BEGIN \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    file << "**SHELL ELEMENT DATA END \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";

    // Contact Surface Data
    file << "**CONTACT SURFACE DATA BEGIN \n";
    file << "**============================================================================== \n";
    file << "**============================================================================== \n";
    file << "**CONTACT SURFACE DATA END \n";
    file << "**============================================================================== \n";

    file.close();

}


void writeVTKmesh( std::string outputfilename, C3t3 c3t3)
{
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    Tr tr = c3t3.triangulation();
	//std::cout << "createUnstructuredGrid" << std::endl;
    int numberofpoints = tr.number_of_vertices();
    //std::cout << "number of points " << numberofpoints << std::endl;
    int numberofcells = c3t3.number_of_cells();
    //std::cout << "number of cells " << numberofcells << std::endl;

	// ========================================================
    // Points!
	vtkPoints * points = vtkPoints::New();
	vtkCellArray * cells = vtkCellArray::New();

	points->SetNumberOfPoints( numberofpoints );
	double * auxPoint = new double[3];
		auxPoint[0] = 0.0;
		auxPoint[1] = 0.0;
		auxPoint[2] = 0.0;

    int i = 0;
    std::map<Point_3, int> V;
    for (Tr::All_vertices_iterator it=tr.all_vertices_begin(); it != tr.all_vertices_end(); ++it)
    {   
        if( i>0 ){	
            auxPoint[0] = it->point()[0];
            auxPoint[1] = it->point()[1];
            auxPoint[2] = it->point()[2];

            points->SetPoint(i-1, auxPoint);

            V[it->point()] = i-1;
        };
	    ++i;
    }

    // Cell!
	vtkSmartPointer<vtkIdTypeArray> idCells = vtkSmartPointer<vtkIdTypeArray>::New();

	idCells->SetNumberOfComponents(4);
	idCells->SetNumberOfTuples( numberofcells );

	vtkIdType* tuple = new vtkIdType[4];

    i = 0;
	C3t3::Vertex_handle v0;
	C3t3::Vertex_handle v1;
	C3t3::Vertex_handle v2;
	C3t3::Vertex_handle v3;

    for (C3t3::Cell_iterator it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it)
    {
	    const Tr::Cell c(*it);
	    v0 = c.vertex(0);
	    v1 = c.vertex(1);
	    v2 = c.vertex(2);
	    v3 = c.vertex(3);
	    // std::cout<<"Point #"<<i<<" : "<<v0->point()<<" ; "<<v1->point()<<" ; "<<v2->point()<<" ; "<<v3->point()<<std::endl;
	    // std::cout<<"Indices: "<<V[v0->point()]<<"/"<<V[v1->point()]<<"/"<<V[v2->point()]<<"/"<<V[v3->point()]<<std::endl;

        //tuple[0] = 4; 
		tuple[0] = V[v0->point()];
		tuple[1] = V[v1->point()];
		tuple[2] = V[v2->point()];
		tuple[3] = V[v3->point()];

        cells->InsertNextCell(4,tuple);

	    ++i;
    }

    /*
	vtkIntArray * data = vtkIntArray::New();
		data->SetName("PectoralBC");
	bool is_inside = false;
	for(int bc = 0; bc< numberofpoints; bc++ )	{
		is_inside = false;
		for(int j=0; j<bC_list.size(); j++){
			if (bc==bC_list[j]) 
			{
					is_inside = true;
					std::cout << "Is INSIDE " << bc << std::endl;
			}
		}

		if(is_inside) data->InsertNextValue(1);
		else data->InsertNextValue(0);
	}
	
	//vtkIdType * id_tissue;
	vtkIntArray * tissues = vtkIntArray::New();
		tissues->SetName("Tissue");
	for( int t=0; t<tissue_list.size(); t++){
		tissues->InsertNextValue(  tissue_list[t] );
	}
    */

	grid->SetPoints(points); // Points
	grid->SetCells(VTK_TETRA, cells); // Cells

	// grid->GetPointData()->SetScalars( data ); // BoundaryConditions
	// grid->GetCellData()->SetScalars( tissues ); // Tissues

	vtkSmartPointer<vtkUnstructuredGridWriter> initialGridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		initialGridWriter->SetInputData( grid );
		initialGridWriter->SetFileName( outputfilename.c_str() );
		initialGridWriter->Update();
		initialGridWriter->Write();

}


void writeMesh(std::string outputfilename, C3t3 c3t3)
{
    if( outputfilename.compare(outputfilename.size()-4,4,".inp")==0 )  {
        std::cout << "Writing .inp file..." << std::endl;
        writeInpFile(outputfilename, c3t3);
    } else if ( outputfilename.compare(outputfilename.size()-4,4,".vtk")==0 ) {
        std::cout << "Writing .vtk file..." << std::endl;
        writeVTKmesh(outputfilename, c3t3);
    } else {
        std::cout << "File name does not correspond to neither inp nor vtk file format" << std::endl;
        std::cout << "Writing .mesh file..." << std::endl;
        outputfilename = outputfilename + ".mesh";
        std::ofstream medit_file( outputfilename.c_str() );
        c3t3.output_to_medit( medit_file );
        std::cout << "Use medit to visualize the mesh!" << std::endl;
    }
}


int main(int argc, char* argv[])
{
    if(argc==1 || strcmp(argv[1],"--h")==0 || strcmp(argv[1],"--help")==0 ) {
        usage();
        return EXIT_SUCCESS;
        }
    
    std::string inputfilename = argv[1];
    std::string outputfilename = argv[2];
    argum argumentos;

    if(argc>3) argumentos = argParser(argc, argv);
    
    std::cout << std::endl;
    std::cout << " Value Outside : " << argumentos.value_outside << std::endl;
    std::cout << " Facet Angle : " << argumentos.facet_angle << std::endl;
    std::cout << " Facet Size : " << argumentos.facet_size << std::endl;
    std::cout << " Facet Distance : " << argumentos.facet_distance << std::endl;
    std::cout << " Cell Radius-Edge Ratio : " << argumentos.cell_radius_edge_ratio << std::endl;
    std::cout << " Cell Size : " << argumentos.cell_size << std::endl;
    std::cout << std::endl;

    ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( inputfilename );
        try{
            reader->Update();
        } catch( itk::ExceptionObject & excp){
            std::cout << "Reader Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }

    // Itk to Vtk
    ItktoVtkType::Pointer itkToVtk = ItktoVtkType::New();
        itkToVtk->SetInput( reader->GetOutput() );
        try{
            itkToVtk->Update();
        } catch( itk::ExceptionObject & excp) {
            std::cout << "Itk To Vtk Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }

    // Vtk to CGAL
    CGAL::Image_3 image = CGAL::read_vtk_image_data( itkToVtk->GetOutput() );
    if( image.image() == 0 ){
        std::cout << "CGAL conversion Exception!"<< std::endl;
        return EXIT_FAILURE;
    }

    // Labeled image to Mesh
    using namespace CGAL::parameters;
    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(
        image,
        value_outside = argumentos.value_outside
    );

    // Mesh Criteria
    /*
    Mesh_criteria criteria(
        facet_angle=30,
        facet_size = 5,
        facet_distance = 1,
        cell_radius_edge_ratio = 30,
        cell_size = 5
    );
    */
    Mesh_criteria criteria(
        facet_angle = argumentos.facet_angle,
        facet_size = argumentos.facet_size,
        facet_distance = argumentos.facet_distance,
        cell_radius_edge_ratio = argumentos.cell_radius_edge_ratio,
        cell_size = argumentos.cell_size
    );

    // Mesh creation 
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    // Mesh Extraction??
    std::cout << std::endl;
    Tr tr = c3t3.triangulation();
    std::cout << "number of points " << tr.number_of_vertices() << std::endl;
    std::cout << "number of elements " << c3t3.number_of_cells() << std::endl;
    std::cout << std::endl;

    // Write mesh
    writeMesh( outputfilename, c3t3);

    return EXIT_SUCCESS;
}



