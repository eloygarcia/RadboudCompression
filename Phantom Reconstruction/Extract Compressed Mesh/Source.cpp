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
if(argc<2)
    {
        std::cout << "Usage: " << argv[0] << " original_mesh compressed_mesh" << std::endl;
        return EXIT_FAILURE;
    }

    /* VTK */
    std::string original_mesh = argv[1];
    std::string output_mesh = argv[2];

    /* original_mesh */
    vtkSmartPointer< vtkUnstructuredGridReader > meshReader = vtkSmartPointer< vtkUnstructuredGridReader>::New();
        meshReader->SetFileName( original_mesh.c_str());
        meshReader->Update();

    // Initial Points
	std::vector<float> initial_points; // = input_mri->getPoints_model();
	vtkPoints * i_points = meshReader->GetOutput()->GetPoints();
	int NoP = meshReader->GetOutput()->GetNumberOfPoints();
	double temp_point[3] = {0.0,0.0,0.0};

    // Final Points

    // Include a try catch HERE!!!
    std::cout << "If the software fails here, please, check the node displacements" << std::endl;
    vtkDataArray* BoundCond = meshReader->GetOutput()->GetPointData()->GetScalars("displacements");//->GetArray(0);

//    std::cout << BoundCond->GetNumberOfTuples() << std::endl;
//    std::cout << meshReader->GetOutput()->GetNumberOfPoints() << std::endl;
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

//	std::cout << "here 1" << std::endl;

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

    std::cout << "Compressed mesh!" << std::endl;

    return EXIT_SUCCESS;
}