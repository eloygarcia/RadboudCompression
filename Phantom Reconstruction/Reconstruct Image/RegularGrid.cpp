#ifndef __RegularGrid_cpp
#define __RegularGrid_cpp

#include "RegularGrid.h"

RegularGrid::RegularGrid(void)
{
	is_points = false;
	is_elements = false;

	is_update = false;

	// ====== Initialize ==========
//	m_spacing = new float[3];
	is_spacing = false;

//	m_origen = new float[3];
	is_origen = false;

//	m_size = new int[3];
	is_size = false;

	is_init = false;

	// ========== Bounding Box ============
	m_boundingBox = new float[6];
	is_boundingBox = false;

	// std::cout << "Inicializando Grid" << std::endl;
}

RegularGrid::~RegularGrid(void)
{
//	delete [] m_origen;
//	delete [] m_spacing;
//	delete [] m_size;
//
//	delete [] m_points;
	delete [] m_boundingBox;
//
//	delete [] m_elements;
	delete [] m_OBBs;

	delete [] m_flags;
	delete [] m_cumsum;
	delete [] m_correspondingElement;

	// std::cout << "Grid Destruida! " << std::endl;ç

	is_points = false;
	is_elements = false;

	is_spacing = false;
	is_origen = false;
	is_size = false;

	is_init = false;
	is_update = false;

	is_boundingBox = false;
}

void RegularGrid::ComputeBoundingBox()
{
	float x_max = -99999.9; float x_min = 99999.9;
	float y_max = -99999.9; float y_min = 99999.9;
	float z_max = -99999.9; float z_min = 99999.9;

	for(int i=0; i<m_numberOfPoints; i++)
	{
		if(m_points[3*i] < x_min) x_min= m_points[3*i];
		if(m_points[3*i] > x_max) x_max= m_points[3*i];

		if(m_points[3*i+1] < y_min) y_min= m_points[3*i+1];
		if(m_points[3*i+1] > y_max) y_max= m_points[3*i+1];

		if(m_points[3*i+2] < z_min) z_min= m_points[3*i+2];
		if(m_points[3*i+2] > z_max) z_max= m_points[3*i+2];
	}

	m_boundingBox[0]=x_min;
	m_boundingBox[1]=x_max;
	m_boundingBox[2]=y_min;
	m_boundingBox[3]=y_max;
	m_boundingBox[4]=z_min;
	m_boundingBox[5]=z_max;

	is_boundingBox = true;

//	std::cout << "BoundingBox: ["<< m_boundingBox[0]<< ", "<< m_boundingBox[1] << ", "<< m_boundingBox[2] <<
//		 ", "<< m_boundingBox[3] << ", "<< m_boundingBox[4] << ", "<< m_boundingBox[5] << std::endl;
}

void RegularGrid::Initialize()
{
	if( !is_origen ){
		ComputeBoundingBox();
		m_origen = new float[3];
			m_origen[0] = m_boundingBox[0]-1;
			m_origen[1] = m_boundingBox[2]-1;
			m_origen[2] = m_boundingBox[4]-1;
			is_origen = true;
	}

	if( !is_spacing ){
		m_spacing = new float[3];
			m_spacing[0] = 1;
			m_spacing[1] = 1;
			m_spacing[2] = 1;
			is_spacing = true;
	};

	if( !is_size ){
		if( is_spacing && is_boundingBox){
			m_size = new int[3];
				m_size[0] = (int)ceil((m_boundingBox[1]+2 - m_boundingBox[0] )/m_spacing[0]);
				m_size[1] = (int)ceil((m_boundingBox[3]+2 - m_boundingBox[2] )/m_spacing[1]);
				m_size[2] = (int)ceil((m_boundingBox[5]+2 - m_boundingBox[4] )/m_spacing[2]);
				m_numberOfPixels = m_size[0] * m_size[1] *  m_size[2];
		}
	}

	is_init = true;
}

float RegularGrid::max( float * idx )
{
	float maximum = -999.999;
	for(int i=0; i<4; i++) { if ( idx[i]> maximum ) maximum = idx[i]; }
	return maximum;
}

float RegularGrid::min(float * idx)
{
	float minimum = 9999.99999;
	for(int i=0; i<4; i++) { if ( idx[i]< minimum ) minimum = idx[i]; }
	return minimum;
}

void RegularGrid::Update()
{
	if( is_points && is_elements ) {
		if( !is_spacing || !is_origen || !is_size ) Initialize();

		if( !is_boundingBox ) ComputeBoundingBox();


//		std::cout  << "Origen: [" << m_origen[0] << ", " << m_origen[1] << ", " << m_origen[2] <<"] " << std::endl;
//		std::cout  << "Size: [" << m_size[0] << ", " << m_size[1] << ", " << m_size[2] <<"] " << std::endl;
//		std::cout  << "spacing: [" << m_spacing[0] << ", " << m_spacing[1] << ", " << m_spacing[2] <<"] " << std::endl;
//		std::cout << " Number Of Pixels = " << m_numberOfPixels << std::endl;

//		std::cout << m_points[0]<< "  " << m_points[1] << "  " << m_points[2] << " .... " << m_points[3*m_numberOfPoints-2] << " " << m_points[3*m_numberOfPoints-1]  << " " << m_points[3*m_numberOfPoints] << std::endl;
//		std::cout << m_elements[0]<< "  " << m_elements[1] << "  " << m_elements[2] << " .... " << m_elements[3*m_numberOfElements-2] << " " << m_elements[3*m_numberOfElements-1]  << " " << m_elements[3*m_numberOfElements] << std::endl;
		// Computing Oriented Bounding Boxes !
			float * temp_x = new float[4];
			float * temp_y = new float[4];
			float * temp_z = new float[4];

			for( int i=0; i<m_numberOfElements; i++){
				temp_x[0] = m_points[ 3* m_elements[4*i] ];
				temp_x[1] = m_points[ 3* m_elements[4*i+1] ];
				temp_x[2] = m_points[ 3* m_elements[4*i+2] ];
				temp_x[3] = m_points[ 3* m_elements[4*i+3] ];
				m_OBBs[6*i] = min( temp_x );
				m_OBBs[6*i+1] = max( temp_x );

				temp_y[0] = m_points[ 3* m_elements[4*i] +1 ];
				temp_y[1] = m_points[ 3* m_elements[4*i+1] +1 ];
				temp_y[2] = m_points[ 3* m_elements[4*i+2] +1 ];
				temp_y[3] = m_points[ 3* m_elements[4*i+3] +1 ];
				m_OBBs[6*i+2] = min( temp_y );
				m_OBBs[6*i+3] = max( temp_y );

				temp_z[0] = m_points[ 3* m_elements[4*i] +2 ];
				temp_z[1] = m_points[ 3* m_elements[4*i+1] +2 ];
				temp_z[2] = m_points[ 3* m_elements[4*i+2] +2 ];
				temp_z[3] = m_points[ 3* m_elements[4*i+3] +2 ];
				m_OBBs[6*i+4] = min( temp_z );
				m_OBBs[6*i+5] = max( temp_z );
			}

			delete [] temp_x;
			delete [] temp_y;
			delete [] temp_z;

//		std::cout << m_OBBs[0]<< "  " << m_OBBs[1] << "  " << m_OBBs[2] << " .... " << m_OBBs[6*m_numberOfElements-2] << " " << m_OBBs[6*m_numberOfElements-1]  << " " << m_OBBs[6*m_numberOfElements] << std::endl;
		// Flags !
			m_flags = new int[ m_numberOfPixels ];
			for( int i=0; i<m_numberOfPixels; i++){	m_flags[i]=0;	}

			int sp_x_min, sp_y_min, sp_z_min;
			int sp_x_max, sp_y_max, sp_z_max;
			int idx = 0;

			for( int i=0; i<m_numberOfElements; i++){
				sp_x_min = (int) floor((m_OBBs[6*i] - m_origen[0])/m_spacing[0]);
				sp_y_min = (int) floor((m_OBBs[6*i+2] - m_origen[1])/m_spacing[1]);
				sp_z_min = (int) floor((m_OBBs[6*i+4] - m_origen[2])/m_spacing[2]);

				sp_x_max = (int) ceil((m_OBBs[6*i+1] - m_origen[0])/m_spacing[0]);
				sp_y_max = (int) ceil((m_OBBs[6*i+3] - m_origen[1])/m_spacing[1]);
				sp_z_max = (int) ceil((m_OBBs[6*i+5] - m_origen[2])/m_spacing[2]);
		
				for( int j_z=sp_z_min; j_z<sp_z_max; j_z++) {
					for( int j_y=sp_y_min; j_y<sp_y_max ; j_y++){
						for( int j_x=sp_x_min; j_x<sp_x_max; j_x++){
							idx = j_x + (m_size[0]*j_y) + (m_size[0]*m_size[1]*j_z);
							//std::cout << idx << " de " << numberOfPixels << "pixels "<< std::endl;
							m_flags[idx]= m_flags[idx]+1;
						}
					}
				}
			}

		// Accumulative sum !!
			m_cumsum  = new int[ m_numberOfPixels ];
			m_cumsum[0] = 0;
			int sum = 0;
			for( int i = 1; i<m_numberOfPixels; i++) {
				sum = sum + m_flags[ i-1 ];
				m_cumsum[i] = sum;
			}

			int maximum_elements = m_cumsum[ m_numberOfPixels-1 ];

		// Alocación de elementos !!
			m_correspondingElement = new int[ maximum_elements ];
			for( int i=0 ; i<maximum_elements; i++) m_correspondingElement[i] = 0;
			for( int i=0; i<m_numberOfPixels; i++){	m_flags[i]=0;	}

			int nim = 0;
			for( int i=0; i<m_numberOfElements; i++){
				sp_x_min = (int) floor((m_OBBs[6*i] - m_origen[0])/m_spacing[0]);
				sp_y_min = (int) floor((m_OBBs[6*i+2] - m_origen[1])/m_spacing[1]);
				sp_z_min = (int) floor((m_OBBs[6*i+4] - m_origen[2])/m_spacing[2]);

				sp_x_max = (int) ceil((m_OBBs[6*i+1] - m_origen[0])/m_spacing[0]);
				sp_y_max = (int) ceil((m_OBBs[6*i+3] - m_origen[1])/m_spacing[1]);
				sp_z_max = (int) ceil((m_OBBs[6*i+5] - m_origen[2])/m_spacing[2]);
		
				for( int j_z=sp_z_min; j_z<sp_z_max; j_z++) {
					for( int j_y=sp_y_min; j_y<sp_y_max ; j_y++){
						for( int j_x=sp_x_min; j_x<sp_x_max; j_x++){
							idx = j_x + (m_size[0]*j_y) + (m_size[0]*m_size[1]*j_z);
							nim = m_cumsum[ idx ] + m_flags[ idx ];

							m_correspondingElement[ nim ] = i;
							m_flags[idx]= m_flags[idx]+1;
						}
					}
				}
			}
			is_update = true;
	}
}

#endif