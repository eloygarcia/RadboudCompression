#include <iostream>
#include <fstream>

// Itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkFlipImageFilter.h"


typedef itk::Image<unsigned short,3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
typedef itk::BinaryBallStructuringElement<unsigned short, 3> StructuringElementType;
typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> DilateType;
typedef itk::BinaryMorphologicalClosingImageFilter< ImageType, ImageType, StructuringElementType> ClosingType;
typedef itk::BinaryMorphologicalOpeningImageFilter< ImageType, ImageType, StructuringElementType> OpeningType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef itk::FlipImageFilter< ImageType > FlipType;

#include "itkPoint.h"
typedef itk::Point<double, ImageType::ImageDimension> PointType;

// Itk-Vtk Glue
#include "itkImageToVTKImageFilter.h"
typedef itk::ImageToVTKImageFilter<ImageType> ItktoVtkType;

// Vtk
#include "vtkImageData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"

#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkCellIterator.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkCellData.h"

#include "vtkTetra.h"

#include "vtkStaticCleanUnstructuredGrid.h"

#include "vtkCellCenters.h"
#include "vtkCellDataToPointData.h"

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
//#include <CGAL/generate_label_weights.h>
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

struct argum{
    argum() {
        value_outside = 0;
        facet_angle = 30;
        facet_size = 2;
        facet_distance = 1;
        cell_radius_edge_ratio = 2;
        cell_size = 2;
        bc_thickness = 0;

        pixel_spacing.push_back(0.273); // X size
        pixel_spacing.push_back(0.273); // Y size
        pixel_spacing.push_back(0.273); // Z size

        bc = "VICTRE";
        side = "R";
    };
    ~argum(){};

    int value_outside;
    double facet_angle;
    double facet_size;
    double facet_distance;
    double cell_radius_edge_ratio;
    double cell_size;
    double bc_thickness;
    std::vector<float> pixel_spacing;
    std::string bc;
    std::string side;
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
        if( strcmp(argv[i], "--pixel_spacing")==0 ){
            argumentos.pixel_spacing[0] = atof(argv[i+1]);
            argumentos.pixel_spacing[1] = atof(argv[i+2]);
            argumentos.pixel_spacing[2] = atof(argv[i+3]);
        };
        if(strcmp(argv[i],"--bc")==0){
            argumentos.bc = argv[i+1];
        };
        if( strcmp(argv[i], "--bc_thickness")==0 ) argumentos.bc_thickness = atof(argv[i+1]);
        if( strcmp(argv[i], "--breast_side")==0 ) argumentos.side = argv[i+1];
    }
    return argumentos;
}

void usage()
{
    argum argumentos;
    std::cout << std::endl;
    std::cout << "inputfilename outputfilename arguments" << std::endl;
    std::cout << std::endl;
    std::cout << "- inputfilename : Name of the original volume (3D image) in itk format -i.e. .nrrd or .mhd" << std::endl;
    std::cout << "- outputfilename :  Name of the mesh file with extesion .inp or .vtk"  << std::endl;
    std::cout << "- arguments :" << std::endl;
    std::cout << " \t --value_outside : (default "<< argumentos.value_outside << ")" << std::endl;
    std::cout << " \t --facet_angle : (default "<< argumentos.facet_angle << ")" << std::endl;
    std::cout << " \t --facet_size : (default "<< argumentos.facet_size << ")" << std::endl;
    std::cout << " \t --facet_distance : (default "<< argumentos.facet_distance << ")" << std::endl;
    std::cout << " \t --cell_radius_edge_ratio : (default "<< argumentos.cell_radius_edge_ratio << ")" << std::endl;
    std::cout << " \t --cell_size : (default "<< argumentos.cell_size << ")" << std::endl;
    std::cout << " \t --pixel_spacing : (default " << argumentos.pixel_spacing[0] << ", " << argumentos.pixel_spacing[1] << ", " << argumentos.pixel_spacing[0] << ")" << std::endl;
    std::cout << " \t --bc : (default VICTRE -i.e. no pectoral muscle-) or  CT" << std::endl;
    std::cout << " \t --bc_thickness : (default "<< argumentos.bc_thickness << ")" << std::endl;
    std::cout << " \t --breast_side : (default "<< argumentos.side << ")" << std::endl;
};

int updatePhysicalInformation(ImageType::Pointer & image, std::vector<float> pixel_spacing, std::string inputfilename)
{
    ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    ImageType::PointType  origin;
        origin[0] = 0.0;
        origin[1] = 0.0;
        origin[2] = 0.0;
    image->SetOrigin( origin );
    ImageType::DirectionType direction;
        direction.SetIdentity();
    image->SetDirection( direction );
    ImageType::SpacingType spacing;
        spacing[0] = pixel_spacing[0];
        spacing[1] = pixel_spacing[1];
        spacing[2] = pixel_spacing[2];
    image->SetSpacing( spacing );

    std::cout << "Writing " << inputfilename.substr(0, inputfilename.length()-5) << "-resampled.nrrd image" << std::endl;

    WriterType::Pointer writer = WriterType::New();
        //writer->SetFileName( inputfilename.substr(0, inputfilename.length()-5) +"-resampled.nrrd" );
        writer->SetFileName( inputfilename );
        writer->SetInput( image );
        writer->UseCompressionOn();
        try{
            writer->Update();
        } catch( itk::ExceptionObject & excp ) {
            std::cout << "Writer Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        };
    return 0;
}

void vtkMeshWithBC(C3t3 c3t3, vtkSmartPointer<vtkUnstructuredGrid>& vtk_umesh)
{
    Tr tr = c3t3.triangulation();
	//std::cout << "createUnstructuredGrid" << std::endl;
    int numberofpoints = tr.number_of_vertices();
    //std::cout << "number of points " << numberofpoints << std::endl;
    int numberofcells = c3t3.number_of_cells();
    //std::cout << "number of cells " << numberofcells << std::endl;

    // temporal vectors for boundary conditions and element labelling
    std::vector<float> i_points;
    std::vector<int> elements;
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

            // temporal nodes vector
            i_points.push_back(auxPoint[0]);
            i_points.push_back(auxPoint[1]);
            i_points.push_back(auxPoint[2]);
            //

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

		// temporal element vector
		elements.push_back( tuple[0] );
		elements.push_back( tuple[1] );
		elements.push_back( tuple[2] );
		elements.push_back( tuple[3] );
		///

        cells->InsertNextCell(4,tuple);

	    ++i;
    }


	vtk_umesh->SetPoints(points); // Points
	vtk_umesh->SetCells(VTK_TETRA, cells); // Cells

    vtkSmartPointer< vtkStaticCleanUnstructuredGrid > cleaningMesh = vtkSmartPointer< vtkStaticCleanUnstructuredGrid>::New();
        cleaningMesh->SetInputData(vtk_umesh);
        cleaningMesh->ToleranceIsAbsoluteOn();
        cleaningMesh->SetTolerance(0.0);
        cleaningMesh->RemoveUnusedPointsOn();
        cleaningMesh->Update();

    std::cout << std::endl;
    std::cout << "Final number of points: "<< cleaningMesh->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "Final number of cells: "<< cleaningMesh->GetOutput()->GetNumberOfCells() << std::endl;

    vtk_umesh = cleaningMesh->GetOutput();
}

void labelingElements(vtkSmartPointer<vtkUnstructuredGrid>& u_grid, ImageType::Pointer image, std::string bc, int bc_thickness)
{
    /* Labeling elements */
    // Set material description in the elements:
    // 0 - background
    // 1 - Fatty
    // 2 - Glandular
    // 3 - skin
    //double tempCenter[3] = {0.0,0.0,0.0};

    vtkSmartPointer<vtkCellCenters> centers = vtkSmartPointer<vtkCellCenters>::New();
        centers->SetInputData(u_grid);
        centers->Update();

     vtkPointSet* pointSet = centers->GetOutput();
     //std::cout << pointSet->GetNumberOfPoints() << std::endl;
     ImageType::IndexType idx;

    double * pt; //[3] ={0.0,0.0,0.0};
    PointType itkpt;
    ImageType::PixelType pixelValue;

    vtkIntArray * data = vtkIntArray::New();
	data->SetName("materials");

/*
// https://discourse.vtk.org/t/get-centroids-in-python-script/5073/4
    cell_Iterator = data.NewCellIterator()
    while cell_Iterator:
        cellID = cell_Iterator.GetCellId()
        cellData = data.GetCell(cellID)
        centroid = cellData.GetCentroid()
        cell_Iterator.GoToNextCell()
  */
     for( int i=0; i<pointSet->GetNumberOfPoints(); i++){
        pt = pointSet->GetPoint(i);
            itkpt[0] = pt[0];
            itkpt[1] = pt[1];
            itkpt[2] = pt[2];
        image->TransformPhysicalPointToIndex( itkpt, idx);
            pixelValue = image->GetPixel(idx);
            if(pixelValue==0) pixelValue=3;
        data->InsertNextValue( pixelValue );

        ImageType::SizeType ss = image->GetLargestPossibleRegion().GetSize();

        //std::cout << "Size: [" << ss[0] <<", "<< ss[1] <<", "<< ss[2] <<"]"<< std::endl;
        //std::cout << "[" << idx[0] << ", "<< idx[1] << ", "<< idx[2] << "] " << std::endl;
        //std::cout << "Pixel Value: " << pixelValue << std::endl;

     };

     u_grid->GetCellData()->SetScalars( data ); // Tissues

    // initialize point data array
    vtkIntArray * bcdata = vtkIntArray::New();
    bcdata->SetName("boundaryConditions");

    if( strcmp(bc.c_str(), "VICTRE")==0){
        std::cout << "Using VICTRE labeling for BC" << std::endl;
          // Boundary conditions:
          // This is performed using point data.
          //
         vtkSmartPointer<vtkCellDataToPointData> cell2point = vtkSmartPointer<vtkCellDataToPointData>::New();
            cell2point->SetInputData(u_grid);
            cell2point->Update();

        // std::cout << cell2point->GetOutput()->GetCellData()->GetNumberOfArrays() << std::endl;
        // std::cout << cell2point->GetOutput()->GetPointData()->GetNumberOfArrays() << std::endl;
        // std::cout << cell2point->GetOutput()->GetPointData()->GetArray(0)->GetNumberOfTuples() << std::endl;
        // std::cout << cell2point->GetOutput()->GetPointData()->GetArray(0)->GetName() << std::endl;
        for( int i=0; i<cell2point->GetOutput()->GetPointData()->GetArray(0)->GetNumberOfTuples(); i++){
            if( cell2point->GetOutput()->GetPointData()->GetArray(0)->GetTuple1(i) > 3) bcdata->InsertNextValue(1);
            else bcdata->InsertNextValue(0);
            //std::cout << cell2point->GetOutput()->GetPointData()->GetArray(0)->GetTuple1(i) << std::endl;
        }
    }else if( strcmp(bc.c_str(), "CT")==0){
        std::cout << "Using CT labeling for BC" << std::endl;
        vtkSmartPointer<vtkPoints> bcpoints = u_grid->GetPoints();
        int numberofpoints = bcpoints->GetNumberOfPoints();
        double * bbox = bcpoints->GetBounds();
        double * point;
        for(int i=0; i<numberofpoints; i++){
            pt = bcpoints->GetPoint(i);
            // std::cout << "[" << pt[0] << ", " << pt[1] << ", " << pt[2] << "]" << std::endl;
            /*
            This code extract the boundary conditions using CT images, that do not contains pectoral muscle.
            In this case, the bc_thickess y define the volume in milimeters where the bc have to be located.
            If there are no bc nodes in the final mesh check (and change if necessary) the axis and limits of the
             */
            if( pt[0] < (bbox[0]+bc_thickness)) bcdata->InsertNextValue(1); // for radboud CT images pt[0] < (bbox[0]+bc_thickness)
                                                                            // for UC Davis (Sarno paper) pt[2] < (bbox[4]+bc_thickness)
            else bcdata->InsertNextValue(0);
        }
    }
    u_grid->GetPointData()->SetScalars( bcdata );

}

int main(int argc, char* argv[])
{
    // arguments!
    argum argumentos;

    if(argc==1 || strcmp(argv[1],"--h")==0 || strcmp(argv[1],"--help")==0 ) {
        usage();
        return EXIT_SUCCESS;
        }
    
    std::string inputfilename = argv[1];
    std::string outputfilename = argv[2];

    if(argc>3) argumentos = argParser(argc, argv);
    
    std::cout << std::endl;
    std::cout << " Value Outside : " << argumentos.value_outside << std::endl;
    std::cout << " Facet Angle : " << argumentos.facet_angle << std::endl;
    std::cout << " Facet Size : " << argumentos.facet_size << std::endl;
    std::cout << " Facet Distance : " << argumentos.facet_distance << std::endl;
    std::cout << " Cell Radius-Edge Ratio : " << argumentos.cell_radius_edge_ratio << std::endl;
    std::cout << " Cell Size : " << argumentos.cell_size << std::endl;
    std::cout << " Image spacing : [" << argumentos.pixel_spacing[0] << ", " << argumentos.pixel_spacing[1] << ", " << argumentos.pixel_spacing[0] << "]" << std::endl;
    std::cout << " Boundary Condition Type : " << argumentos.bc << std::endl;
    if(strcmp(argumentos.bc.c_str(), "CT")==0)    std::cout << " Boundary Condition Thickness : " << argumentos.bc_thickness << std::endl;
    std::cout << " Breast Side: " << argumentos.side << std::endl;
    std::cout << std::endl;

    if((strcmp(argumentos.bc.c_str(), "CT")==0) && (argumentos.bc_thickness==0)){
        std::cout << "Check BC conditions. BC thickness cannot be 0 for CT images." << std::endl;
        return EXIT_FAILURE;
    }

    // Itk Image
    ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( inputfilename );
        try{
            reader->Update();
        } catch( itk::ExceptionObject & excp){
            std::cout << "Reader Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }

    ImageType::Pointer img = reader->GetOutput();
    ImageType::SpacingType sp = img->GetSpacing();

    // Breast Side: Flipping to the R view
    if(strcmp(argumentos.side.c_str(), "R")!=0 || strcmp(argumentos.side.c_str(), "L")==0){
        FlipType::Pointer flipFilter = FlipType::New();
        FlipType::FlipAxesArrayType flipAxes;
            flipAxes[0] = false;
            flipAxes[1] = false;
            flipAxes[2] = true;

        flipFilter->SetInput(img);
        flipFilter->SetFlipAxes(flipAxes);
        flipFilter->Update();

        img=flipFilter->GetOutput();
    }


    // Physical information
    std::string subname = inputfilename.substr( inputfilename.length()-5, inputfilename.length());
    //std::cout << subname.c_str() << std::endl;

    // if((strcmp(subname.c_str(),".nrrd")!=0 && strcmp(subname.c_str(),".mhd")!=0) ||
    //    (sp[0]!=0.273 && sp[1]!=0.273 && sp[2]!=0.273)
    //    )
    //{
        int e = updatePhysicalInformation(img, argumentos.pixel_spacing, inputfilename);
        if(e!=0) return EXIT_FAILURE;
    //};


    // Thresholding
    ThresholdType::Pointer threshold = ThresholdType::New();
        threshold->SetInput( img );
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(1);
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(255);
        threshold->Update();
   ImageType::SizeType ss = threshold->GetOutput()->GetLargestPossibleRegion().GetSize();
   std::cout << "Size: [" << ss[0] <<", "<< ss[1] <<", "<< ss[2] <<"]"<< std::endl;

    // Closing
    std::cout << "Closing element!!" << std::endl;
    StructuringElementType structuringElement2;
        // structuringElement.SetRadius(1); // 3x3 structuring element
        structuringElement2.SetRadius(2);
        structuringElement2.CreateStructuringElement();

    //ClosingType::Pointer morphofilter = ClosingType::New();
    OpeningType::Pointer morphofilter = OpeningType::New();
        morphofilter->SetInput(threshold->GetOutput());
        morphofilter->SetKernel(structuringElement2);
        morphofilter->SetBackgroundValue(0);
        morphofilter->SetForegroundValue(1);
        morphofilter->Update();

    // Dilate
    std::cout << "Including Dilate ball element!!!" << std::endl;
    StructuringElementType structuringElement;
        // structuringElement.SetRadius(1); // 3x3 structuring element
        structuringElement.SetRadius(1);
        structuringElement.CreateStructuringElement();

    DilateType::Pointer binaryDilate = DilateType::New();
        binaryDilate->SetKernel(structuringElement);
        binaryDilate->SetInput(morphofilter->GetOutput());
        binaryDilate->SetDilateValue(1);
        binaryDilate->Update();




    std::cout << "itk to vtk" << std::endl;
    // Itk to Vtk
    ItktoVtkType::Pointer itkToVtk = ItktoVtkType::New();
        //itkToVtk->SetInput( threshold->GetOutput() );
        itkToVtk->SetInput( binaryDilate->GetOutput() );
        //itkToVtk->SetInput( morphofilter->GetOutput() );
        try{
            itkToVtk->Update();
        } catch( itk::ExceptionObject & excp) {
            std::cout << "Itk To Vtk Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }
     //std::cout << "vtk output tuples:" << itkToVtk->GetOutput()->GetPointData()->GetScalars()->GetNumberOfTuples() << std::endl;
     //std::cout << "itk voxles: [" << ss[0] * ss[1] * ss[2] <<"]"<< std::endl;


    std::cout << "vtk to cgal" << std::endl;
    // Vtk to CGAL
    CGAL::Image_3 image = CGAL::read_vtk_image_data( itkToVtk->GetOutput() );
    if( image.image() == 0 ){
        std::cout << "CGAL conversion Exception!"<< std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "labeled Image to mesh" << std::endl;
    // Labeled image to Mesh
    using namespace CGAL::parameters;

    //const float sigma = 5.f;
    //CGAL::Image_3 img_weights = CGAL::Mesh_3::generate_label_weights(image, sigma);

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(
        image,
        value_outside = argumentos.value_outside,
        relative_error_bound = 1e-6
    );

    // Mesh Criteria
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
    std::cout << "Initial number of points " << tr.number_of_vertices() << std::endl;
    std::cout << "Initial number of elements " << c3t3.number_of_cells() << std::endl;
    std::cout << std::endl;

    // Write mesh
    //writeMesh( outputfilename, c3t3);
    vtkSmartPointer<vtkUnstructuredGrid> vtk_umesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkMeshWithBC(c3t3, vtk_umesh);

    // Labeling elements
    labelingElements(vtk_umesh,img, argumentos.bc, argumentos.bc_thickness);

    /// Writing the mesh
	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		//writer->SetInputData( cleaningMesh->GetOutput() );
		writer->SetInputData( vtk_umesh );
		writer->SetFileName( outputfilename.c_str() );
		writer->Update();
		writer->Write();

    return EXIT_SUCCESS;
}