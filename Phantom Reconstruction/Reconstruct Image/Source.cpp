#include "RegularGrid.h"
#include "NewProjection.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGridWriter.h"

#include "vtkPoints.h"
#include "vtkCellArray.h"

typedef itk::Image< unsigned char, 3> ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

void extractNodesAndElements(vtkUnstructuredGrid * mesh, std::vector<float> & nodes, std::vector<int> & elements)
{
    // W
    nodes.clear();
    elements.clear();

    // Nodes Extraction
    vtkPoints * i_points = mesh->GetPoints();
    int NoP = mesh->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};

    for (int i=0; i<NoP; i++)
	{
		i_points->GetPoint(i,temp_point);

		nodes.push_back( temp_point[0] );
		nodes.push_back( temp_point[1] );
		nodes.push_back( temp_point[2] );
	};


    // Elements Extraction
	int accumm = 0;
	int minimum = 999999999;
	vtkCellArray * i_cells = mesh->GetCells(); // LA PUTA QUE LOS PARIO !!
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

}

int main( int argc, char* argv[])
{
    if(argc<4)
    {
        std::cout << "Usage: " << argv[0] << " original_mesh compressed_mesh original_itk_image phantom_image" << std::endl;
        return EXIT_FAILURE;
    }

    std::string original_mesh = argv[1];
    std::string compressed_mesh = argv[2];
    std::string image = argv[3];
    std::string output_image = argv[4];

    /* VTK */
     /* original_mesh */
    vtkSmartPointer< vtkUnstructuredGridReader > originalMeshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        originalMeshReader->SetFileName( original_mesh.c_str());
        originalMeshReader->Update();

    std::vector<float> initial_points;
    std::vector<int> elements;
    extractNodesAndElements( originalMeshReader->GetOutput(), initial_points, elements);

    /* compressed_mesh*/
    vtkSmartPointer< vtkUnstructuredGridReader > compressedMeshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        compressedMeshReader->SetFileName( compressed_mesh.c_str());
        compressedMeshReader->Update();

    std::vector<float> final_points;
    std::vector<int> unused_elements;
    extractNodesAndElements( compressedMeshReader->GetOutput(), final_points, unused_elements);


    /* ITK */
    /* Read original image */
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

    /* IMAGE RECONSTRUCTION */
    float spa[3] = {0.0,0.0,0.0}; //just preparing the function
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

    /* WRITE COMPRESSED BREAST IMAGE */
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