#include "cuda_runtime.h"
#include "device_launch_parameters.h"


__device__ void computealpha( const int  ind, const float origen, const float voxelSize,const float p1, const float p2, float &alpha) {	// alpha = [ X_plane(i) - X1] / [X2-X1] ;
		alpha = ((ind*voxelSize + origen)-p1)/(p2-p1); };

__device__ void Determinante( float* a0, float* a1, float* a2, float & det )
{
	float w[3] = {0.0,0.0,0.0};
		w[0] = a1[1]*a2[2] - a1[2]*a2[1];
		w[1] = -a1[0]*a2[2] + a1[2]*a2[0];
		w[2] = a1[0]*a2[1] - a1[1]*a2[0];
	
	det = a0[0]*w[0] + a0[1]*w[1] + a0[2]*w[2];
	//return det;
}

__device__ void triLinearInterpolator(const float* imagepointer, const int* size, const float* pixel, float & value) // es un pixel continuo, con lo cual, vale
{
//	float value = 0.0;
	int x0, y0, z0, x1, y1, z1;
		x0 = (int) (floorf(pixel[0])); x1=x0+1; // x1 = (int) (ceilf(pixel[0]));
		y0 = (int) (floorf(pixel[1])); y1=y0+1; // y1 = (int) (ceilf(pixel[1]));
		z0 = (int) (floorf(pixel[2])); z1=z0+1; // z1 = (int) (ceilf(pixel[2]));

	float xd = (pixel[0]-x0)/(x1-x0);
	float yd = (pixel[1]-y0)/(y1-y0);
	float zd = (pixel[2]-z0)/(z1-z0);

	int idx_1 = 0;		int idx_2 = 0;

		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z0);
	float c00 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z1);
	float c01 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z0);
	float c10 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z1);
	float c11 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );

		float c0 = c00 * (1-yd) + c10 * yd;
		float c1 = c01 * (1-yd) + c11 * yd;

		float c = c0 * (1-zd) + c1 * zd;

		value = c;
}

__device__ void computeBaricentricCoordinates(float* posicionPunto, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* baricentricCoordinates, bool & is_inside)
{
	is_inside = false;

	float vap[3] = {0.0,0.0,0.0};
		vap[0] = posicionPunto[0]-vertex_0[0];
		vap[1] = posicionPunto[1]-vertex_0[1];
		vap[2] = posicionPunto[2]-vertex_0[2];

	float vbp[3] = {0.0,0.0,0.0};
		vbp[0] = posicionPunto[0]-vertex_1[0];
		vbp[1] = posicionPunto[1]-vertex_1[1];
		vbp[2] = posicionPunto[2]-vertex_1[2];

	float vab[3] = {0.0,0.0,0.0};
		vab[0] = vertex_1[0]-vertex_0[0];
		vab[1] = vertex_1[1]-vertex_0[1];
		vab[2] = vertex_1[2]-vertex_0[2];

	float vac[3] = {0.0,0.0,0.0};
		vac[0] = vertex_2[0]-vertex_0[0];
		vac[1] = vertex_2[1]-vertex_0[1];
		vac[2] = vertex_2[2]-vertex_0[2];
	
	float vad[3] = {0.0,0.0,0.0};
		vad[0] = vertex_3[0]-vertex_0[0];
		vad[1] = vertex_3[1]-vertex_0[1];
		vad[2] = vertex_3[2]-vertex_0[2];

	float vbc[3] = {0.0,0.0,0.0};
		vbc[0] = vertex_2[0]-vertex_1[0];
		vbc[1] = vertex_2[1]-vertex_1[1];
		vbc[2] = vertex_2[2]-vertex_1[2];

	float vbd[3] = {0.0,0.0,0.0};
		vbd[0] = vertex_3[0]-vertex_1[0];
		vbd[1] = vertex_3[1]-vertex_1[1];
		vbd[2] = vertex_3[2]-vertex_1[2];

	float va6 = 0;
        Determinante(vbp, vbd, vbc, va6);
	float vb6 = 0;
        Determinante(vap, vac, vad, vb6);
	float vc6 = 0;
        Determinante(vap, vad, vab, vc6);
	float vd6 = 0;
        Determinante(vap, vab, vac, vd6);
    float v_temp = 0;
        Determinante(vab, vac, vad, v_temp);
	float v6 = 1/abs(v_temp);

	baricentricCoordinates[0] = va6*v6;
	baricentricCoordinates[1] = vb6*v6;
	baricentricCoordinates[2] = vc6*v6;
	baricentricCoordinates[3] = vd6*v6;

	
	if( (baricentricCoordinates[0]>=0 && baricentricCoordinates[0] <= 1) &&
		(baricentricCoordinates[1]>=0 && baricentricCoordinates[1] <= 1) &&
		(baricentricCoordinates[2]>=0 && baricentricCoordinates[2] <= 1) &&
		(baricentricCoordinates[3]>=0 && baricentricCoordinates[3] <= 1) &&
		baricentricCoordinates[0]+baricentricCoordinates[1]+baricentricCoordinates[2]+baricentricCoordinates[3]<1.01 &&
		baricentricCoordinates[0]+baricentricCoordinates[1]+baricentricCoordinates[2]+baricentricCoordinates[3]>0.99 )
		is_inside=true;
	
	//return is_inside;
}

__device__ void computeCartessianCoordinates(float* baricentricCoordinates, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* posicionPunto)
{
	posicionPunto[0] = baricentricCoordinates[0]*vertex_0[0] + baricentricCoordinates[1]*vertex_1[0] + baricentricCoordinates[2]*vertex_2[0] + baricentricCoordinates[3]*vertex_3[0];
	posicionPunto[1] = baricentricCoordinates[0]*vertex_0[1] + baricentricCoordinates[1]*vertex_1[1] + baricentricCoordinates[2]*vertex_2[1] + baricentricCoordinates[3]*vertex_3[1];
	posicionPunto[2] = baricentricCoordinates[0]*vertex_0[2] + baricentricCoordinates[1]*vertex_1[2] + baricentricCoordinates[2]*vertex_2[2] + baricentricCoordinates[3]*vertex_3[2];
}

/*
__global__ void kernel_projection(const int* dev_3d_size, const float* dev_3d_spacing, const float* dev_3d_origen, const float* dev_3d_imagepointer,
								  const int* dev_2d_size, const float* dev_2d_spacing, const float* dev_2d_origen, float* dev_2d_imagepointer,
								  const float* source)
{
	int i = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	int numberOfPixels2d = dev_2d_size[0] * dev_2d_size[1]; // *size2d[2];

	if( i<numberOfPixels2d){

	// Fila y columna de la imagen !
	int row = (int) floorf( i / dev_2d_size[0] ); 
	int col = (int) ( i - (row*dev_2d_size[0]) ); 
	

	// Posición del pixel del detector !!
	float x2 = dev_2d_origen[0] + (col * dev_2d_spacing[0] );
	float y2 = dev_2d_origen[1] + (row * dev_2d_spacing[1] );

	// Vector de dirección !!
	float vect[3];
		vect[0] = x2 - source[0];
		vect[1] = y2 - source[1];
		vect[2] = dev_2d_origen[2] - source[2];  // Revisar este punto.

	// Distancia de la fuente al detector !!
	float xa = pow(vect[0],2);
	float ya = pow(vect[1],2);
	float za = pow(vect[2],2);

	float dist12 = sqrt( xa + ya + za );

	// Cálculo del alpha en Z !!
	float temp = 0.0f;
	computealpha( 0, dev_3d_origen[2], dev_3d_spacing[2], source[2], dev_2d_origen[2], temp );

	// Resolver la ecuación de la recta !!
	float temp_dist[3] = {0.0f, 0.0f, 0.0f};
	float point[3] = {0.0f, 0.0f, 0.0f};
	float step = 0.0005f;
	float t=temp;

	float pixel[3] = {0.0f, 0.0f, 0.0f};

	float value = 0.0f;
	float length = 0.0f;

	// Longitud del paso !!
	temp_dist[0] = pow(step*vect[0],2);
	temp_dist[1] = pow(step*vect[1],2);
	temp_dist[2] = pow(step*vect[2],2);
	float l_step = sqrt(temp_dist[0] + temp_dist[1] + temp_dist[2]);

	while( t<1 )
	{
		// Calculo del siguiente punto en la recta
		point[0] = source[0] + t * vect[0];
		point[1] = source[1] + t * vect[1];
		point[2] = source[2] + t * vect[2];
		// Posición del voxel parcial !!
		pixel[0] = (point[0] - dev_3d_origen[0]) / dev_3d_spacing[0];
		pixel[1] = (point[1] - dev_3d_origen[1]) / dev_3d_spacing[1];
		pixel[2] = (point[2] - dev_3d_origen[2]) / dev_3d_spacing[2];
		// Interpolación trilienal
		if((pixel[0]<0 || pixel[0]>dev_3d_size[0]-1) || (pixel[1]<0 || pixel[1]>dev_3d_size[1]-1) || (pixel[2]<0 || pixel[2]>dev_3d_size[2]-1))  value = 0;
		else triLinearInterpolator( dev_3d_imagepointer, dev_3d_size, pixel, value);
		// Longitud acumulada !!
		length += (1000 * value * l_step);
		// Nuevo punto!! 
		t+=step;
	}

	if(length>0 & length < 65535 ) dev_2d_imagepointer[ i ] = length;
	else dev_2d_imagepointer[ i ] = 0.0f;
	
	}
}
*/

__global__ void kernel_projection(const int* dev_3d_size, const float* dev_3d_spacing, const float* dev_3d_origen, const unsigned char* dev_3d_imagepointer,
								  const float* dev_i_points, const float* dev_f_points, const int* dev_elements,
								  const float* dev_grid_origen, const float* dev_grid_spacing, const int* dev_grid_size,
								  const int* dev_flags, const int* dev_cumsum, const int* dev_correspondingElements,
								  const int* dev_2d_size, const float* dev_2d_spacing, const float* dev_2d_origen, unsigned char* dev_2d_imagepointer) 
{
	long int i = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	//int numberOfPixels2d = dev_2d_size[0] * dev_2d_size[1] * dev_2d_size[2];

		int index_z = (int)  (i /((long int)dev_2d_size[0]*(long int)dev_2d_size[1]));
		int index_y = (int) ((i - ((long int) index_z*(long int)dev_2d_size[0]*(long int)dev_2d_size[1]))/(long int) dev_2d_size[0]);
		int index_x = (int) ( i - ((long int) index_z*(long int)dev_2d_size[0]*(long int)dev_2d_size[1]) - ((long int) index_y*(long int)dev_2d_size[0]));

		float temp_position[3] = {0.0,0.0,0.0};
			temp_position[0] = dev_2d_origen[0] + dev_2d_spacing[0]*(index_x + 0.5);
			temp_position[1] = dev_2d_origen[1] + dev_2d_spacing[1]*(index_y + 0.5);
			temp_position[2] = dev_2d_origen[2] + dev_2d_spacing[2]*(index_z + 0.5);

		float temp_VoxelGrid[3] = {0.0,0.0,0.0};
			temp_VoxelGrid[0] = floorf((temp_position[0] - dev_grid_origen[0])/dev_grid_spacing[0]);
			temp_VoxelGrid[1] = floorf((temp_position[1] - dev_grid_origen[1])/dev_grid_spacing[1]);
			temp_VoxelGrid[2] = floorf((temp_position[2] - dev_grid_origen[2])/dev_grid_spacing[2]);
		long int correspondingIndex = 0;
			correspondingIndex = temp_VoxelGrid[0] + (dev_grid_size[0]*temp_VoxelGrid[1]) + (dev_grid_size[0]*dev_grid_size[1]*temp_VoxelGrid[2]);
			
		// El pixel le tengo. Puedo saber qué elementos hay
		int index_start = dev_cumsum[ correspondingIndex ]; // Este i no vale! "i" es para la imagen!" hay que calcular otro para la grid regular
		int number_of_elements_here = dev_flags[ correspondingIndex ];
		// int count=0;
		int element_index_number = 0;
		int element_number = 0;

		int temp_element[4] = {0,0,0,0};
		float temp_vertex_0[3] = {0.0,0.0,0.0};
		float temp_vertex_1[3] = {0.0,0.0,0.0};
		float temp_vertex_2[3] = {0.0,0.0,0.0};
		float temp_vertex_3[3] = {0.0,0.0,0.0};

		float old_vertex_0[3] = {0.0,0.0,0.0};
		float old_vertex_1[3] = {0.0,0.0,0.0};
		float old_vertex_2[3] = {0.0,0.0,0.0};
		float old_vertex_3[3] = {0.0,0.0,0.0};

		float baricentricCoordinates[4] = {0.0,0.0,0.0,0.0};
		bool is_inside = false;

		float pixel[3] ={0.0,0.0,0.0};
		float sss_position[3] = {0.0,0.0,0.0};

		long int index = 0;
		float pixel_value = 0;

		if( number_of_elements_here != 0 ){
		for( int j=0; j<number_of_elements_here; j++ ){
			element_index_number = index_start + j;
			element_number = dev_correspondingElements[ element_index_number ];
			
				temp_element[0] = dev_elements[ 4*element_number ];
				temp_element[1] = dev_elements[ 4*element_number +1 ];
				temp_element[2] = dev_elements[ 4*element_number +2 ];
				temp_element[3] = dev_elements[ 4*element_number +3 ];
		
				temp_vertex_0[0] = dev_f_points[ 3*temp_element[0]  ];
				temp_vertex_0[1] = dev_f_points[ 3*temp_element[0] +1 ];
				temp_vertex_0[2] = dev_f_points[ 3*temp_element[0] +2 ];

				temp_vertex_1[0] = dev_f_points[ 3*temp_element[1] ];
				temp_vertex_1[1] = dev_f_points[ 3*temp_element[1] +1 ];
				temp_vertex_1[2] = dev_f_points[ 3*temp_element[1] +2 ];
				
				temp_vertex_2[0] = dev_f_points[ 3*temp_element[2] ];
				temp_vertex_2[1] = dev_f_points[ 3*temp_element[2] +1 ];
				temp_vertex_2[2] = dev_f_points[ 3*temp_element[2] +2 ];
				
				temp_vertex_3[0] = dev_f_points[ 3*temp_element[3] ];
				temp_vertex_3[1] = dev_f_points[ 3*temp_element[3] +1 ];
				temp_vertex_3[2] = dev_f_points[ 3*temp_element[3] +2 ];

			computeBaricentricCoordinates(temp_position, temp_vertex_0, temp_vertex_1, temp_vertex_2, temp_vertex_3, baricentricCoordinates, is_inside);

			if( is_inside ){
				// dev_2d_imagepointer[i] = 0; // hasta aquí llega!!!
				//if(dev_2d_imagepointer[i] == 0 ) 
				//{
					// Primer código !
//					 final_image[i] = element_number+1;
					// Segundo código !
					// 1. Hay que acalcular primero la posición cartesiana en el espacio de la imagen de entrada
				//	std::cout << "Temp. Element: [" << temp_element[0] << ", " << temp_element[1] << ", " << temp_element[2] << ", " << temp_element[3] << "] " << std::endl;
					
						old_vertex_0[0] = dev_i_points[ 3*temp_element[0]  ];
						old_vertex_0[1] = dev_i_points[ 3*temp_element[0] +1 ];
						old_vertex_0[2] = dev_i_points[ 3*temp_element[0] +2 ];
		
						old_vertex_1[0] = dev_i_points[ 3*temp_element[1] ];
						old_vertex_1[1] = dev_i_points[ 3*temp_element[1] +1 ];
						old_vertex_1[2] = dev_i_points[ 3*temp_element[1] +2 ];
					
						old_vertex_2[0] = dev_i_points[ 3*temp_element[2] ];
						old_vertex_2[1] = dev_i_points[ 3*temp_element[2] +1 ];
						old_vertex_2[2] = dev_i_points[ 3*temp_element[2] +2 ];
									
						old_vertex_3[0] = dev_i_points[ 3*temp_element[3] ];
						old_vertex_3[1] = dev_i_points[ 3*temp_element[3] +1 ];
						old_vertex_3[2] = dev_i_points[ 3*temp_element[3] +2 ];

					computeCartessianCoordinates( baricentricCoordinates, old_vertex_0, old_vertex_1, old_vertex_2, old_vertex_3, sss_position);

						// 2. Localizamos el pixel continuo e interpolamos el valor. 
					
					//	pixel[0] = (sss_position[0]-origenCuda[0])/spacingCuda[0];
					//	pixel[1] = (sss_position[1]-origenCuda[1])/spacingCuda[1];
					//	pixel[2] = (sss_position[2]-origenCuda[2])/spacingCuda[2];
						

						pixel[0] = (sss_position[0]-dev_3d_origen[0])/dev_3d_spacing[0];
						pixel[1] = (sss_position[1]-dev_3d_origen[1])/dev_3d_spacing[1];
						pixel[2] = (sss_position[2]-dev_3d_origen[2])/dev_3d_spacing[2];

												//Checkeando la información ! 
//					std::cout << "Temporal position. [" << temp_position[0] <<", " << temp_position[1] << ", " << temp_position[2] << "] " << std::endl;
//					std::cout << "    coordenadas baricentricas: [" << baricentricCoordinates[0] << ", " << baricentricCoordinates[1] << ", " << baricentricCoordinates[2] << ", " << baricentricCoordinates[3] << "] " << std::endl;
//					std::cout << "Reconstructed position. [" << sss_position[0] <<", " << sss_position[1] << ", " << sss_position[2] << "] " << std::endl;
//					std::cout << " El pixel es [" << pixel[0] <<", " << pixel[1] << ", " << pixel[2] << "]  de size = [" << sizeCuda[0] <<", " << sizeCuda[1] << ", "<< sizeCuda[2] << "] " << std::endl;
//					std::cout << std::endl;
					// Sleep(250);

					index = (long int)((long int)floorf(pixel[0]) + ((long int)dev_3d_size[0] * (long int)floorf(pixel[1]) ) + ((long int)dev_3d_size[0]*(long int)dev_3d_size[1]*(long int)floorf(pixel[2])));
					//dev_2d_imagepointer[i] = (unsigned char)index;
					if(index>0) dev_2d_imagepointer[i]=dev_3d_imagepointer[index];

  					
					// dev_2d_imagepointer[i] = (unsigned char)dev_3d_imagepointer[index];
					 
					//triLinearInterpolator(dev_3d_imagepointer, dev_3d_size, pixel, pixel_value);
//						dev_2d_imagepointer[i] = pixel_value;
					//else std::cout << "Fuera de image!" << std::endl;


//						trilinearinterpolation( image_input , sizeCuda, pixel, value);
//						nearestNeighborInterpolation( image_input , sizeCuda, pixel, value);
						// Hay que cambiar esta interpolacion por nearest neigbors

//						final_image[i] = value;
				/*	*/

				//}else{
				//	count++;
				//	std::cout << i << ".  Más de un elemento en el voxel !! " << final_image[i] << " and " << element_number << std::endl;
				//	std::cout << "    coordenadas baricentricas: [" << barcor[0] << ", " << barcor[1] << ", " << barcor[2] << ", " << barcor[3] << "] " << std::endl;
				//	std::cout << "    coordenadas baricentricas: [" << baricentricCoordinates[0] << ", " << baricentricCoordinates[1] << ", " << baricentricCoordinates[2] << ", " << baricentricCoordinates[3] << "] " << std::endl;
				//	std::cout << std::endl;
				//	Sleep(500);

				// } // endif( final_image[i]==0 )
			} // endif( is_inside == true )		
		} // endfor
		} // endif( number_of_elements_here!=0 )		

		// dev_2d_imagepointer[i] = (char)number_of_elements_here;
		//dev_2d_imagepointer[i] = 1.0f;
}


__global__ void fill_dos( char * imagepointer)
{
	long unsigned int i = blockDim.x * blockIdx.x +threadIdx.x;
	imagepointer[i]=2;
}