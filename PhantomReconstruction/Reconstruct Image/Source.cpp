#include "RegularGrid.h"
#include "NewProjection.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::Image< unsigned char, 3> ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGridWriter.h"

#include "vtkPoints.h"
#include "vtkCellArray.h"

void extractGrid(vtkSmartPointer<vtkUnstructuredGrid> grid, 
				 std::vector<float> &initial_points,
				 std::vector<float> &final_points, 
				 std::vector<int> &elements)
{
	// =============================== EXTRACCION DE LA MALLA ===========================
	// Points
//	std::vector<double> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = grid->GetPoints();
	int NoP = grid->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};
	for (int i=0; i<NoP; i++)
	{
		i_points->GetPoint(i,temp_point);

		initial_points.push_back( temp_point[0] );
		initial_points.push_back( temp_point[1] );
		initial_points.push_back( temp_point[2] );
	}
	
	// std::vector<double> final_points;
	// for(int xx = 0; xx< initial_points.size(); xx++){	final_points.push_back( 0 ); }
	
	// Elements
	int accumm = 0;
	int minimum = 999999999;
//	std::vector<int> elements; // = input_mri->getElements_model();
	vtkCellArray * i_cells = grid->GetCells(); // LA PUTA QUE LOS PARIO !!
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

	// Displacements
	vtkSmartPointer<vtkDataArray> vtk_disp = grid->GetPointData()->GetArray(0);

	double pt[3] = {0.0, 0.0, 0.0};
	double d[3] = {0.0, 0.0, 0.0};
	double temp[3] = {0.0, 0.0, 0.0};
		pt[0] = 0.0; pt[1] = 0.0; pt[2] = 0.0;
		d[0] = 0.0; d[1] = 0.0; d[2] = 0.0;
		temp[0] = 0.0; temp[1] = 0.0; temp[2] = 0.0;

	final_points.clear();

	for(int i=0; i<i_points->GetNumberOfPoints(); i++)
	{
		i_points->GetPoint(i, pt);
		vtk_disp->GetTuple(i, d);

		temp[0] = pt[0] + d[0];
		temp[1] = pt[1] + d[1];
		temp[2] = pt[2] + d[2];

		// f_points[ 3*i ] = temp[0];
		// f_points[ 3*i +1 ] = temp[1];
		// f_points[ 3*i +2 ] = temp[2];

		final_points.push_back( temp[0] );
		final_points.push_back( temp[1] );
		final_points.push_back( temp[2] );
	}

};


int main( int argc, char* argv[])
{
    if(argc<4)
    {
        std::cout << "Usage: " << argv[0] << " niftysim_mesh itk_image output_mesh output_image" << std::endl;
        return EXIT_FAILURE;
    }

    /* VTK */
    std::string mesh = argv[1];
    std::string image = argv[2];
    std::string output_mesh = argv[3];
    std::string output_image = argv[4];

    vtkSmartPointer< vtkUnstructuredGridReader > meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( mesh.c_str());
        meshReader->Update();

/*
    std::vector<float> initial_points;
    std::vector<float> final_points;
    std::vector<int> elements;
    extractGrid( meshReader->GetOutput(), initial_points, 
                        final_points, elements)
*/
    // Initial Points
	std::vector<float> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = meshReader->GetOutput()->GetPoints();
	int NoP = meshReader->GetOutput()->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};

    // Final Points
    vtkDataArray* BoundCond = meshReader->GetOutput()->GetPointData()->GetScalars("displacements");//->GetArray(0);
    std::cout << BoundCond->GetNumberOfTuples() << std::endl;
    std::cout << meshReader->GetOutput()->GetNumberOfPoints() << std::endl;
    double * disp;
    std::vector<float> final_points;

	for (int i=0; i<NoP; i++)
	{
		i_points->GetPoint(i,temp_point);

		initial_points.push_back( temp_point[0] );
		initial_points.push_back( temp_point[1] );
		initial_points.push_back( temp_point[2] );

        disp = BoundCond->GetTuple(i);
        
        final_points.push_back(temp_point[0]+disp[0]);
        final_points.push_back(temp_point[1]+disp[1]);
        final_points.push_back(temp_point[2]+disp[2]);
	}
	
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

    // Create Unstructured Grid
    vtkSmartPointer< vtkUnstructuredGrid> u_grid = vtkSmartPointer< vtkUnstructuredGrid>::New();
    vtkPoints * points = vtkPoints::New();

    points->SetNumberOfPoints( NoP );
    for(int i=0; i<NoP; i++)
    {
        temp_point[0] = final_points[3*i];
        temp_point[1] = final_points[(3*i)+1];
        temp_point[2] = final_points[(3*i)+2];

        points->SetPoint(i, temp_point);
    }

    u_grid->SetPoints( points );
    u_grid->SetCells( VTK_TETRA, i_cells);

    // Write moved mesh
    vtkSmartPointer< vtkUnstructuredGridWriter> meshWriter = vtkSmartPointer< vtkUnstructuredGridWriter>::New();
        meshWriter->SetInputData(u_grid);
        meshWriter->SetFileName( output_mesh.c_str() );
        meshWriter->Update();



    /* ITK */
    ReaderType::Pointer imageReader = ReaderType::New();
        imageReader->SetFileName( image );
        try {
            imageReader->Update();
        } catch(itk::ExceptionObject & excp) {
            std::cout << "ImageReader exception!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }
    
    ImageType::PointType origen = imageReader->GetOutput()->GetOrigin();
	ImageType::SpacingType spacing = imageReader->GetOutput()->GetSpacing();
	ImageType::SizeType size = imageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
	ImageType::DirectionType direction = imageReader->GetOutput()->GetDirection();

    std::cout << "Origen : [" << origen[0] << ", " << origen[1] << ", " << origen[2] << "] " << std::endl;
	std::cout << "Spacing : [" << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << "] " << std::endl;
	std::cout << "Size : [" << size[0] << ", " << size[1] << ", " << size[2] << "] " << std::endl;
	std::cout << "Direction : [" << direction[0][0] << ", " << direction[0][1] << ", " << direction[0][2] << "] " << std::endl;
	std::cout << "\t" << direction[1][0] << ", " << direction[1][1] << ", " << direction[1][2] << "] " << std::endl;
	std::cout << "\t" << direction[2][0] << ", " << direction[2][1] << ", " << direction[2][2] << "] " << std::endl;

    float spa[3] = {0.0,0.0,0.0};   
        spa[0] = spacing[0];
        spa[1] = spacing[1];
        spa[2] = spacing[2];
    
    ImageType::Pointer new_image = ImageType::New();
    NewProjection * compressedImage = new NewProjection();
        compressedImage->getCompressedImage( initial_points, 
                                             final_points, 
                                             elements,
                                             imageReader->GetOutput(),
                                             new_image,
                                             spa);

	WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( output_image );
		writer->SetInput( new_image );
		writer->UseCompressionOn();
		try{
			writer->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cout << "Writer Exception Object!" << std::endl;
			std::cout << excp << std::endl;
			return EXIT_FAILURE;
		}

    return EXIT_SUCCESS;
}