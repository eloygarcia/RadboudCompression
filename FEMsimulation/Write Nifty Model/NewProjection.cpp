#ifndef __NewProjection_cpp
#define __NewProjection_cpp

#include "NewProjection.h"

NewProjection::NewProjection(void)
{
	importfilter = ImportFilterType::New();

	storing = new int[8];
	best = 0;
	temp_value = 0;

	// function getCompresedImage
	bounding_box = new float[6];
	for(unsigned int i=0; i<6; i++){
		bounding_box[i]=0;
	}
	n = 0;
	m = 0;
	
	m_origen = new float[3];
		m_origen[0] = 0.0;
		m_origen[1] = 0.0;
		m_origen[2] = 0.0;
	m_spacing = new float[3];
		m_spacing[0]=1;
		m_spacing[1]=1;
		m_spacing[2]=1;
	m_size = new unsigned int[3];
		m_size[0] = 1.0;
		m_size[1] = 1.0;
		m_size[2] = 1.0;
	numberOfPixels = 0;

	origenCuda = new float[3];
		origenCuda[0] = 0.0;
		origenCuda[1] = 0.0;
		origenCuda[2] = 0.0;
	sizeCuda = new unsigned int[3];
		sizeCuda[0] = 0;
		sizeCuda[1] = 0;
		sizeCuda[2] = 0;
	spacingCuda = new float[3];
		spacingCuda[0] = 0.0;
		spacingCuda[1] = 0.0;
		spacingCuda[2] = 0.0;

	temp_x = new float[4];
		temp_x[0] = 0.0;
		temp_x[1] = 0.0;
		temp_x[2] = 0.0;
		temp_x[3] = 0.0;
	temp_y = new float[4];
		temp_y[0] = 0.0;
		temp_y[1] = 0.0;
		temp_y[2] = 0.0;
		temp_y[3] = 0.0;
	temp_z = new float[4];
		temp_z[0] = 0.0;
		temp_z[1] = 0.0;
		temp_z[2] = 0.0;
		temp_z[3] = 0.0;

	sp_x_min = 0;
	sp_y_min = 0;
	sp_z_min = 0;
	
	sp_x_max = 0;
	sp_y_max = 0;
	sp_z_max = 0;

	idx = 0;
	nim = 0;

	voxel_position = new float[3];
		voxel_position[0] = 0.0;
		voxel_position[1] = 0.0;
		voxel_position[2] = 0.0;
	//sum = 0;
	maximum_elements = 0;

	index_z = 0;
	index_y = 0;
	index_x = 0;

	temp_position = new float[3];
		temp_position[0] = 0.0;
		temp_position[1] = 0.0;
		temp_position[2] = 0.0;
	baricentricCoordinates = new float[4];
		baricentricCoordinates[0] = 0.0;
		baricentricCoordinates[1] = 0.0;
		baricentricCoordinates[2] = 0.0;
		baricentricCoordinates[3] = 0.0;
	barcor = new float[4];
		barcor[0] = 0.0;
		barcor[1] = 0.0;
		barcor[2] = 0.0;
		barcor[3] = 0.0;

	index_start = 0;
	number_of_elements_here = 0;
	count=0;

	element_index_number = 0;
	element_number = 0;

	temp_element = new int[4];
		temp_element[0] = 0;
		temp_element[1] = 0;
		temp_element[2] = 0;
		temp_element[3] = 0;
	temp_vertex_0 = new float[3];
		temp_vertex_0[0] = 0.0;
		temp_vertex_0[1] = 0.0;
		temp_vertex_0[2] = 0.0;
	temp_vertex_1 = new float[3];
		temp_vertex_1[0] = 0.0;
		temp_vertex_1[1] = 0.0;
		temp_vertex_1[2] = 0.0;
	temp_vertex_2 = new float[3];
		temp_vertex_2[0] = 0.0;
		temp_vertex_2[1] = 0.0;
		temp_vertex_2[2] = 0.0;
	temp_vertex_3 = new float[3];
		temp_vertex_3[0] = 0.0;
		temp_vertex_3[1] = 0.0;
		temp_vertex_3[2] = 0.0;

	is_inside = false;

	old_vertex_0 = new float[3];
		old_vertex_0[0] = 0.0;
		old_vertex_0[1] = 0.0;
		old_vertex_0[2] = 0.0;
	old_vertex_1 = new float[3];
		old_vertex_1[0] = 0.0;
		old_vertex_1[1] = 0.0;
		old_vertex_1[2] = 0.0;
	old_vertex_2 = new float[3];
		old_vertex_2[0] = 0.0;
		old_vertex_2[1] = 0.0;
		old_vertex_2[2] = 0.0;
	old_vertex_3 = new float[3];
		old_vertex_3[0] = 0.0;
		old_vertex_3[1] = 0.0;
		old_vertex_3[2] = 0.0;

	sss_position = new float[3];
		sss_position[0] = 0.0;
		sss_position[1] = 0.0;
		sss_position[2] = 0.0;
	value=0.0;
	pixel = new float[3];
		pixel[0] = 0.0;
		pixel[1] = 0.0;
		pixel[2] = 0.0;

	//function Determinante
	det = 0.0;
	//function TrilinearInterpolation.
	x0=0; y0=0; z0=0; x1=0; y1=0; z1=0;
	xd = 0; yd = 0; zd =0;
	idx_1 = 0; idx_2 = 0;
	index=0;
	//funcion nearestNeighborInterpolation.
	val000=0; val001=0; val010=0; val011=0; val100=0; val101=0; val110=0; val111=0;
	sumVal=0;

	c00=0; c10=0; c01=0; c11=0;
	c0=0; c1=0;
	c=0;
	//function computeBarycentricCoordinates
	V=0;
	v1=0; v2=0; v3=0; v4=0; 
	is_inside = false;
	//function max
	maximum = -99999.9999;
	//function min
	minimum = 999999.99999;
	//function get_boundingBox
	x_max = -9999999; x_min = 99999999;
	y_max = -9999999; y_min = 99999999;
	z_max = -9999999; z_min = 99999999;
}

NewProjection::~NewProjection(void)
{
/*	delete corresponding_element;
	delete flags;
	delete temp_flags;

	delete bb;
	*/
	delete m_origen;
	delete m_spacing;
	delete m_size;

//	delete final_image;

}

// private
float NewProjection::Determinante( float* v0, float* v1, float* v2, float* v3 )
{
	//float det = 0.0;
	det  = (v0[2]*v1[1]*v2[0]*1 - v0[1]*v1[2]*v2[0]*1 - 
			v0[2]*v1[0]*v2[1]*1 + v0[0]*v1[2]*v2[1]*1 + 
			v0[1]*v1[0]*v2[2]*1 - v0[0]*v1[1]*v2[2]*1 - 
			v0[2]*v1[1]*1*v3[0] + v0[1]*v1[2]*1*v3[0] + 
			v0[2]*1*v2[1]*v3[0] - 1*v1[2]*v2[1]*v3[0] - 
			v0[1]*1*v2[2]*v3[0] + 1*v1[1]*v2[2]*v3[0] + 
			v0[2]*v1[0]*1*v3[1] - v0[0]*v1[2]*1*v3[1] - 
			v0[2]*1*v2[0]*v3[1] + 1*v1[2]*v2[0]*v3[1] + 
			v0[0]*1*v2[2]*v3[1] - 1*v1[0]*v2[2]*v3[1] - 
			v0[1]*v1[0]*1*v3[2] + v0[0]*v1[1]*1*v3[2] + 
			v0[1]*1*v2[0]*v3[2] - 1*v1[1]*v2[0]*v3[2] - 
			v0[0]*1*v2[1]*v3[2] + 1*v1[0]*v2[1]*v3[2])/6;

	return det;
}

void NewProjection::trilinearinterpolation(const pixelType* imagepointer, const int* size, const float* pixel, pixelType & value) // es un pixel continuo, con lo cual, vale
{
//	float value = 0.0;
c=0;
	// int x0, y0, z0, x1, y1, z1;
		x0 = (int) (floor(pixel[0])); x1 = x0 +1; // (int) (ceil(pixel[0]));
		y0 = (int) (floor(pixel[1])); y1 = y0 +1; //(int) (ceil(pixel[1]));
		z0 = (int) (floor(pixel[2])); z1 = z0 +1; // (int) (ceil(pixel[2]));

	//float xd = (pixel[0]-x0)/(x1-x0);
	//float yd = (pixel[1]-y0)/(y1-y0);
	//float zd = (pixel[2]-z0)/(z1-z0);
	xd = (pixel[0]-x0)/(x1-x0);
	yd = (pixel[1]-y0)/(y1-y0);
	zd = (pixel[2]-z0)/(z1-z0);

	idx_1 = 0;	idx_2 = 0;
	c00=0; c10=0; c01=0; c11=0;
	c0=0; c1=0;
	c=0;	

	if( (xd>0) && ( xd<size[0]) && (yd>0) && (yd<size[1]) && (zd>0) && (zd<size[2]))
	{
	// int idx_1 = 0;		int idx_2 = 0;
		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z0);
	//float c00 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
	c00 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z1);
	//float c01 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
	c01 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z0);
	//float c10 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
	c10 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z1);
	//float c11 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
	c11 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );

		//float c0 = c00 * (1-yd) + c10 * yd;
		//float c1 = c01 * (1-yd) + c11 * yd;
		//
		//float c = c0 * (1-zd) + c1 * zd;
		
		c0 = c00 * (1-yd) + c10 * yd;
		c1 = c01 * (1-yd) + c11 * yd;
		
		c = c0 * (1-zd) + c1 * zd;
	}

		value = c;
}

void NewProjection::nearestNeighborInterpolation( const pixelType * imagepointer, const int* s_size, const float* pixel, pixelType & value)
{
	valido = 1;
	index = -1;

	sumVal=0;

	x0 = (int) (floor(pixel[0])); x1 = x0 +1; // (int) (ceil(pixel[0]));
	y0 = (int) (floor(pixel[1])); y1 = y0 +1; //(int) (ceil(pixel[1]));
	x0 = (int) (floor(pixel[2])); z1 = z0 +1; // (int) (ceil(pixel[2]));

//	for( int i=0; i<8; i++) storing[i] = 0;
	
	// Valor 000
		index = x0 + (s_size[0] * y0) + (s_size[0]*s_size[1]*z0);
		val000 = imagepointer[ index ];
	// Valor 001
		index = x1 + (s_size[0] * y0) + (s_size[0]*s_size[1]*z0);
		val001 = imagepointer[ index ];
	// Valor 010
		index = x0 + (s_size[0] * y1) + (s_size[0]*s_size[1]*z0);
		val010 = imagepointer[ index ];
	// Valor 011
		index = x1 + (s_size[0] * y1) + (s_size[0]*s_size[1]*z0);
		val011 = imagepointer[ index ];
	// Valor 100
		index = x0 + (s_size[0] * y0) + (s_size[0]*s_size[1]*z1);
		val100 = imagepointer[ index ];
	// Valor 101
		index = x1 + (s_size[0] * y0) + (s_size[0]*s_size[1]*z1);
		val101 = imagepointer[ index ];
	// Valor 110
		index = x0 + (s_size[0] * y1) + (s_size[0]*s_size[1]*z1);
		val110 = imagepointer[ index ];
	// Valor 111
		index = x1 + (s_size[0] * y1) + (s_size[0]*s_size[1]*z1);
		val111 = imagepointer[ index ];

	sumVal = val000 + val001 + val010 + val011 + val100 + val101 + val110 + val111;

	if( sumVal<4 ) valido = 0;
	else if( (sumVal>3) && (sumVal<12) ) valido = 1;
	else if( sumVal>11 ) valido = 2;


		// if( (index<(s_size[0]*s_size[1]*s_size[2])) && (0<index)) valido = imagepointer[index];
		// else std::cout << "Fuera de image!" << std::endl;

/*			if(valido<8) storing[valido]++;
		index = x1 + (size[0] * y0) + (size[0]*size[1]*z0);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x0 + (size[0] * y0) + (size[0]*size[1]*z1);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x1 + (size[0] * y0) + (size[0]*size[1]*z1);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x0 + (size[0] * y1) + (size[0]*size[1]*z0);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x1 + (size[0] * y1) + (size[0]*size[1]*z0);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x0 + (size[0] * y1) + (size[0]*size[1]*z1);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;
		index = x1 + (size[0] * y1) + (size[0]*size[1]*z1);
			valido = imagepointer[index];
			if(valido<8) storing[valido]++;

			best = 0;
			temp_value = 0;
	for( int i=0; i<8; i++) {
		if(storing[i]>best){
				best=storing[i];
				temp_value = (unsigned short)i;
		}
	}

	value = temp_value;
	*/
			value = valido;
}

bool NewProjection::computeBaricentricCoordinates(float* posicionPunto, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* baricentricCoordinates)
{
	baricentricCoordinates[0] = 0.0;
	baricentricCoordinates[1] = 0.0;
	baricentricCoordinates[2] = 0.0;
	baricentricCoordinates[3] = 0.0;

	//float V = 0;
	V = abs(Determinante( vertex_0, vertex_1, vertex_2, vertex_3 ));

	//float v1 = 0; 
	v1 = Determinante( posicionPunto, vertex_1, vertex_2, vertex_3);
	//float v2 = 0;
	v2 = Determinante( vertex_0, posicionPunto, vertex_2, vertex_3);
	//float v3 = 0;
	v3 = Determinante( vertex_0, vertex_1, posicionPunto, vertex_3);
	//float v4 = 0;
	v4 = Determinante( vertex_0, vertex_1, vertex_2 ,posicionPunto);
		
	//bool is_inside = false;
	is_inside = false;

	if( ((v1/V)>= -.1 && (v1/V)<= 1.1) && ((v2/V)>= -.1 && (v2/V)<= 1.1) && ((v3/V)>= -.1 && (v3/V)<=1.1) && ((v4/V)>= -.1 && (v4/V)<= 1.1) )
	//if( ((v1/V)>=0 && (v1/V)<= 1) && ((v2/V)>=0 && (v2/V)<= 1) && ((v3/V)>=0 && (v3/V)<=1) && ((v4/V)>=0 && (v4/V)<=1) )
	{
		/*
		std::cout << std::endl;
		std::cout << " punto : [" << posicionPunto[0] << ", " << posicionPunto[1] << ", " << posicionPunto[2] << "] " << std::endl;
		std::cout << "vertex_0 : [" << vertex_0[0] << ", " << vertex_0[1] << ", " << vertex_0[2] << "] " << std::endl;
		std::cout << "vertex_1 : [" << vertex_1[0] << ", " << vertex_1[1] << ", " << vertex_1[2] << "] " << std::endl;
		std::cout << "vertex_2 : [" << vertex_2[0] << ", " << vertex_2[1] << ", " << vertex_2[2] << "] " << std::endl;
		std::cout << "vertex_3 : [" << vertex_3[0] << ", " << vertex_3[1] << ", " << vertex_3[2] << "] " << std::endl;
		std::cout << std::endl;

		std::cout << "Baricentric Coordinates : [" << v1/V <<", " << v2/V << ", " << v3/V << ", " << v4/V << "] " << std::endl;
		*/
		is_inside = true;
		/**/
		baricentricCoordinates[0] = v1/V; 
		baricentricCoordinates[1] = v2/V;
		baricentricCoordinates[2] = v3/V;
		baricentricCoordinates[3] = v4/V;

	}

	// std::cout << "Baricentric Coordinates : [" << v1/V <<", " << v2/V << ", " << v3/V << ", " << v4/V << "] " << std::endl;

	// Sleep(2000);

	return is_inside;
}

void NewProjection::computeCartessianCoordinates(float* baricentricCoordinates, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* posicionPunto)
{
	posicionPunto[0] = baricentricCoordinates[0]*vertex_0[0] + baricentricCoordinates[1]*vertex_1[0] + baricentricCoordinates[2]*vertex_2[0] + baricentricCoordinates[3]*vertex_3[0];
	posicionPunto[1] = baricentricCoordinates[0]*vertex_0[1] + baricentricCoordinates[1]*vertex_1[1] + baricentricCoordinates[2]*vertex_2[1] + baricentricCoordinates[3]*vertex_3[1];
	posicionPunto[2] = baricentricCoordinates[0]*vertex_0[2] + baricentricCoordinates[1]*vertex_1[2] + baricentricCoordinates[2]*vertex_2[2] + baricentricCoordinates[3]*vertex_3[2];
}

float NewProjection::max(float * idx)
{
	//float maximum = -99999.9999;
	maximum = -999.999;
	for(int i=0; i<4; i++) { if ( idx[i]> maximum ) maximum = idx[i]; }
	return maximum;
}

float NewProjection::min(float * idx)
{
	//float minimum = 999999.99999;
	minimum = 9999.999;
	for(int i=0; i<4; i++) { if ( idx[i]< minimum ) minimum = idx[i]; }
	return minimum;
}

void NewProjection::get_boundingBox(float* points, int numberOfPoints, float* bounding_box)
{
	//double x_max = -9999999; double x_min = 99999999;
	//double y_max = -9999999; double y_min = 99999999;
	//double z_max = -9999999; double z_min = 99999999;

	x_max = -9999999; x_min = 99999999;
	y_max = -9999999; y_min = 99999999;
	z_max = -9999999; z_min = 99999999;

	for(int i=0; i<numberOfPoints; i++)
	{
		if(points[3*i] < x_min) x_min= points[3*i];
		if(points[3*i] > x_max) x_max= points[3*i];

		if(points[3*i+1] < y_min) y_min= points[3*i+1];
		if(points[3*i+1] > y_max) y_max= points[3*i+1];

		if(points[3*i+2] < z_min) z_min= points[3*i+2];
		if(points[3*i+2] > z_max) z_max= points[3*i+2];
	}

	bounding_box[0]=x_min;
	bounding_box[1]=x_max;
	bounding_box[2]=y_min;
	bounding_box[3]=y_max;
	bounding_box[4]=z_min;
	bounding_box[5]=z_max;

//	std::cout << "BoundingBox: ["<< bounding_box[0]<< ", "<< bounding_box[1] << ", "<< bounding_box[2] <<
//		 ", "<< bounding_box[3] << ", "<< bounding_box[4] << ", "<< bounding_box[5] << std::endl;
}

// public:

void NewProjection::getCompressedImage(std::vector<double> initial_points, std::vector<double> deformed_points, 
						std::vector<int> elements, ImageType::Pointer image_original, 
						ImageType::Pointer &compressed_image, float* im_spacing)
{
	ImageType::DirectionType direction;
		direction.SetIdentity();
	image_original->SetDirection(direction);

	std::cout << std::endl;
	std::cout << "Comienza la compression por medio de coordenadas baricéntricas !" << std::endl;
	std::cout << std::endl;

	/* Point ! */
	n = initial_points.size()/3;
	float * arr = new float[ initial_points.size() ];
	for(int i=0; i<initial_points.size(); i++) { arr[i] = (float)initial_points[i]; };

	float * new_arr = new float[ deformed_points.size() ];
	for(int i=0; i<deformed_points.size(); i++) { new_arr[i] = (float)deformed_points[i]; };
	
	std::cout << " Constante n = " << n << std::endl;

	// Bounding Box
	float * bounding_box_0 = new float[6];
	for(int i=0; i<6; i++)		bounding_box_0[i]=0;

	get_boundingBox(arr,n,bounding_box_0);
	
	std::cout << std::endl;
	std::cout << "Initial BoundingBox: ["<< bounding_box_0[0]<< ", "<< bounding_box_0[1] << ", "<< bounding_box_0[2] <<
		 ", "<< bounding_box_0[3] << ", "<< bounding_box_0[4] << ", "<< bounding_box_0[5] << std::endl;
	std::cout<< std::endl;

	for(int i=0; i<6; i++)		bounding_box[i]=0;
	get_boundingBox(new_arr,n,bounding_box);

	std::cout << "BoundingBox: ["<< bounding_box[0]<< ", "<< bounding_box[1] << ", "<< bounding_box[2] <<
		 ", "<< bounding_box[3] << ", "<< bounding_box[4] << ", "<< bounding_box[5] << std::endl;
	std::cout<< std::endl;
	
	/* Elements ! */
	m = elements.size()/4;
	int min_cha = 9999;
	int * cha = new int[ elements.size() ];
	for(int i=0; i<elements.size(); i++)	{
			cha[i] = elements[i]; 
//			if(cha[i]<min_cha) min_cha = cha[i];
	}
	
	std::cout << " Constante m = " << m << std::endl;
//	std::cout << " Mincha = " << min_cha  << std::endl;

	/* Leo la imagen original que compone la malla*/
	ImageType::PointType origen_input = image_original->GetOrigin();
		origenCuda[0] = origen_input[0];
		origenCuda[1] = origen_input[1];
		origenCuda[2] = origen_input[2];
	ImageType::SizeType size_input = image_original->GetLargestPossibleRegion().GetSize();
		sizeCuda[0] = size_input[0];
		sizeCuda[1] = size_input[1];
		sizeCuda[2] = size_input[2];
	ImageType::SpacingType spacing_input = image_original->GetSpacing();
		spacingCuda[0] = spacing_input[0];
		spacingCuda[1] = spacing_input[1];
		spacingCuda[2] = spacing_input[2];
	ImageType::DirectionType directioCuda = image_original->GetDirection();
	image_input = image_original->GetBufferPointer();

	std::cout<<std::endl;
	std::cout << " -- Image Input -- " << std::endl;
	std::cout << "Origen : [" << origenCuda[0] << ", " << origenCuda[1] << ", " << origenCuda[2] << "] " << std::endl;
	std::cout << "Size : [" << sizeCuda[0] << ", " << sizeCuda[1] << ", " << sizeCuda[2] << "] " << std::endl;
	std::cout << "Spacing : [" << spacingCuda[0] << ", " << spacingCuda[1] << ", " << spacingCuda[2] << "] " << std::endl;
	std::cout << "Delimitation: \n\t" << origenCuda[0] <<", "<< origenCuda[0]+(sizeCuda[0]*spacingCuda[0]) << std::endl;
	std::cout << "\t" << origenCuda[1] <<", "<< origenCuda[1]+(sizeCuda[1]*spacingCuda[1]) << std::endl;
	std::cout << "\t" << origenCuda[2] <<", "<< origenCuda[2]+(sizeCuda[2]*spacingCuda[2]) << std::endl;
	std::cout << "Direction: \n"<< directioCuda << std::endl;
	std::cout << std::endl;
	
/* GRID ! */
	/* Compute bounding boxes */
	float* bb = new float[ 6*m ];

	// OOBB para todos los elementos 
	for( int i=0; i<m; i++){		// Puede ser que esto sea lo que pasa en la reconstrucción de la malla ??
		temp_x[0] = new_arr[ 3*cha[4*i] ];
		temp_x[1] = new_arr[ 3*cha[4*i+1] ];
		temp_x[2] = new_arr[ 3*cha[4*i+2] ];
		temp_x[3] = new_arr[ 3*cha[4*i+3] ];
		bb[6*i] = min( temp_x );
		bb[6*i+1] = max( temp_x );

		temp_y[0] = new_arr[ 3*cha[4*i] +1 ];
		temp_y[1] = new_arr[ 3*cha[4*i+1] +1 ];
		temp_y[2] = new_arr[ 3*cha[4*i+2] +1 ];
		temp_y[3] = new_arr[ 3*cha[4*i+3] +1 ];
		bb[6*i+2] = min( temp_y );
		bb[6*i+3] = max( temp_y );

		temp_z[0] = new_arr[ 3*cha[4*i] +2 ];
		temp_z[1] = new_arr[ 3*cha[4*i+1] +2 ];
		temp_z[2] = new_arr[ 3*cha[4*i+2] +2 ];
		temp_z[3] = new_arr[ 3*cha[4*i+3] +2 ];
		bb[6*i+4] = min( temp_z );
		bb[6*i+5] = max( temp_z );
	}
	
	/* ORGANIZANDO LA GRID REGULAR !*/
		m_origen[0] = bounding_box[0]-1;
		m_origen[1] = bounding_box[2]-1;
		m_origen[2] = bounding_box[4]-1;
		//m_origen[0] = 0;
		//m_origen[1] = 0;
		//m_origen[2] = bounding_box[4] - 2.;

		m_spacing[0] = 1;
		m_spacing[1] = 1;
		m_spacing[2] = 1;

		//spacing[0] = 5;
		//spacing[1] = 5;
		//spacing[2] = 5;

		m_size[0] = (int)ceil((bounding_box[1]+2-bounding_box[0])/m_spacing[0]);
		m_size[1] = (int)ceil((bounding_box[3]+2-bounding_box[2])/m_spacing[1]);
		m_size[2] = (int)ceil((bounding_box[5]+2-bounding_box[4])/m_spacing[2]);
		
		//m_size[0] = 280; // 140 mm
		//m_size[1] = 460; // 230 mm
		//m_size[2] = (int)floor((bounding_box[5]+4-bounding_box[4])/m_spacing[2]);;

	numberOfPixels = m_size[0]*m_size[1]*m_size[2];

	int * flags = new int[ numberOfPixels ];
	if( numberOfPixels > 0 ) {
		for( int i=0; i<numberOfPixels; i++){
			flags[i]=0;
		}
	}

	std::cout << std::endl;
	std::cout << "NumberOfPixels : " << numberOfPixels << std::endl;
	std::cout << "Bounding box: [" << bounding_box[0] << ", " <<bounding_box[1] << ", " <<bounding_box[2] << ", " <<bounding_box[3] << ", " <<bounding_box[4] << ", " <<bounding_box[5] << "] " << std::endl;
	std::cout << "Origen : [" << m_origen[0] << ", " << m_origen[1] << ", " << m_origen[2] << "] " << std::endl;
	std::cout << "Spacing : [" << m_spacing[0] << ", " << m_spacing[1] << ", " << m_spacing[2] << "] " << std::endl;
	std::cout << "Size : [" << m_size[0] << ", " << m_size[1] << ", " << m_size[2] << "] " << std::endl;
	std::cout << std::endl;

	for( int i=0; i<m; i++){
		sp_x_min = (int) floor((bb[6*i] - m_origen[0])/m_spacing[0]);	
		sp_y_min = (int) floor((bb[6*i+2] - m_origen[1])/m_spacing[1]);
		sp_z_min = (int) floor((bb[6*i+4] - m_origen[2])/m_spacing[2]);

		sp_x_max = (int) ceil((bb[6*i+1] - m_origen[0])/m_spacing[0]);
		sp_y_max = (int) ceil((bb[6*i+3] - m_origen[1])/m_spacing[1]);
		sp_z_max = (int) ceil((bb[6*i+5] - m_origen[2])/m_spacing[2]);

		for( int j_z=sp_z_min; j_z<sp_z_max; j_z++) {
			for( int j_y=sp_y_min; j_y<sp_y_max; j_y++){
				for( int j_x=sp_x_min; j_x<sp_x_max; j_x++){
					
					idx = j_x + (m_size[0]*j_y) + (m_size[0]*m_size[1]*j_z);
					//std::cout << idx << " de " << numberOfPixels << "pixels "<< std::endl;
					flags[idx]= flags[idx]+1;

				}
			}
		}
	}

	int * cumsum = new int[numberOfPixels+1];
		cumsum[0] = 0;
	int sum = 0;
	for( int i=1; i<numberOfPixels+1; i++){
		sum = sum + flags[i-1];
		cumsum[i] = sum;
	}

	//std::cout << "Cum sum = " << cumsum[0] <<", " << cumsum[1] <<", " << cumsum[2] << " [...] " << cumsum[numberOfPixels-1]  <<", " << cumsum[numberOfPixels] <<", " << cumsum[numberOfPixels+1] << std::endl;

	// int maximum_elements = cumsum[numberOfPixels];
	maximum_elements = cumsum[numberOfPixels-1];
	std::cout << std::endl;
	std::cout << "Maximum elements = " << maximum_elements << std::endl;
	std::cout << std::endl;

	// int* correspondingElement = new int[ maximum_elements ];
	correspondingElement = new int[ maximum_elements ];
	for( int i=0; i<maximum_elements; i++ ){
		correspondingElement[i]=0;
	}

	for( int i=0; i<numberOfPixels; i++){
		flags[i]=0;
	}

	for( int i=0; i<m; i++){
		sp_x_min = (int) floor((bb[6*i] - m_origen[0])/m_spacing[0]);
		sp_y_min = (int) floor((bb[6*i+2] - m_origen[1])/m_spacing[1]);
		sp_z_min = (int) floor((bb[6*i+4] - m_origen[2])/m_spacing[2]);

		sp_x_max = (int) ceil((bb[6*i+1] - m_origen[0])/m_spacing[0]);
		sp_y_max = (int) ceil((bb[6*i+3] - m_origen[1])/m_spacing[1]);
		sp_z_max = (int) ceil((bb[6*i+5] - m_origen[2])/m_spacing[2]);
		
		for( int j_z=sp_z_min; j_z<sp_z_max; j_z++) {
			for( int j_y=sp_y_min; j_y<sp_y_max; j_y++){
				for( int j_x=sp_x_min; j_x<sp_x_max; j_x++){

					idx = j_x + (m_size[0]*j_y) + (m_size[0]*m_size[1]*j_z);
					//int nim = cumsum[idx]+temp_flags[idx];
					nim = cumsum[idx]+flags[idx];
					
					correspondingElement[ nim ] = i;
					//temp_flags[idx]= temp_flags[idx]+1;
					flags[idx]= flags[idx]+1;

				}
			}
		}
	}


	// HASTA AQUI LLEGA !!
/* Create deformed image. */
	// size, origen, spacing.

	// el paddel esta entre Y y Z, así que habría que calcular el offset para cada uno de ellos
	// y ajustar el paddel.

	// 2.8mm thick screening compression paddle		// 180x240 mm^2
	// 1 mm thick breast support		// 240x290 mm^2

	// float y_offset = (240 - (bounding_box[3] - bounding_box[2]))/2; // Y es el lado largo 
		float y_offset = 0;
	// float z_offset = (180 - (bounding_box[5] - bounding_box_0[4]))/2; // Z es el lado corto

	// Hay que definir un nuevo spacing y size, utilizaremos el espacio delimitado por el modelo  igual que antes
	float im_origen[3] = {0.0,0.0,0.0};
		im_origen[0] = bounding_box[0]-.5;// m_origen[0]; // compression en X	--> Breast support en X
		// im_origen[1] = m_origen[1];
		im_origen[1] = m_origen[1]-.5;
		// im_origen[2] = m_origen[2];
		im_origen[2] = bounding_box[4]-.5;
	
	int im_size[3] = {0,0,0};
		im_size[0] = (int)ceil((bounding_box[1]+2-bounding_box[0])/im_spacing[0]);
		
		im_size[1] = (int)ceil((bounding_box[3]+2-bounding_box[2])/im_spacing[1]);
		//im_size[1] = (int)ceil( 240 /im_spacing[1]); // 250 mm / spacing_y

		im_size[2] = (int)ceil((bounding_box[5]+2-bounding_box[4])/im_spacing[2]);
		// im_size[2] = (int)ceil((bounding_box[5]-(bounding_box_0[4]+1))/im_spacing[2]);
		// im_size[2] = (int)ceil( 180 /im_spacing[2]); // 180 mm /spacing_z;

	int im_numberOfPixels = im_size[0] * im_size[1] * im_size[2];


	final_image = new pixelType[ im_numberOfPixels ]; // Cambiar cuando rehagas la imagen !!!  NO int SINO float !!
	for( int i=0; i<im_numberOfPixels; i++) {
		final_image[i] = 0;
	}

	float voxel_position[3] = {0.0,0.0,0.0};
	int temp_pixelGG[3] = {0.0,0.0,0.0};
	int index_GG = 0;

	for( int i=0; i<im_numberOfPixels; i++){
		index_z = (int) floor( i /(im_size[0]*im_size[1]));
		index_y = (int) floor((i - (index_z*im_size[0]*im_size[1]))/im_size[0]);
		index_x = (int) ( i - (index_z*im_size[0]*im_size[1]) - (index_y*im_size[0]));

		//float temp_position[3] = {0.0,0.0,0.0};
			temp_position[0] = im_origen[0] + im_spacing[0]*(index_x + 0.5);
			temp_position[1] = im_origen[1] + im_spacing[1]*(index_y + 0.5);
			temp_position[2] = im_origen[2] + im_spacing[2]*(index_z + 0.5);

			// Tengo la posición en la grid pequeña, y ahora tengo que hallar la posicion en la grid grande!
			temp_pixelGG[0] = (int) floor((temp_position[0] - m_origen[0]) / m_spacing[0] );
			temp_pixelGG[1] = (int) floor((temp_position[1] - m_origen[1]) / m_spacing[1] );
			temp_pixelGG[2] = (int) floor((temp_position[2] - m_origen[2]) / m_spacing[2] );

			if( (temp_pixelGG[0]>0) && (temp_pixelGG[0]<m_size[0]) &&
				(temp_pixelGG[1]>0) && (temp_pixelGG[1]<m_size[1]) &&
				(temp_pixelGG[2]>0) && (temp_pixelGG[2]<m_size[2]) )
			{

				index_GG = (temp_pixelGG[2]*m_size[0]*m_size[1]) + (temp_pixelGG[1]*m_size[0]) + temp_pixelGG[0];

			// El pixel le tengo. Puedo saber qué elementos hay
			index_start = cumsum[ index_GG ];
			number_of_elements_here = flags[ index_GG ];
			count=0;

			if( number_of_elements_here != 0 ){
			for( int j=0; j<number_of_elements_here; j++ ){
				element_index_number = index_start + j;
				element_number = correspondingElement[ element_index_number ];
				// int temp_element[4] = {0,0,0,0};
					temp_element[0] = cha[ 4*element_number ];
					temp_element[1] = cha[ 4*element_number +1 ];
					temp_element[2] = cha[ 4*element_number +2 ];
					temp_element[3] = cha[ 4*element_number +3 ];

				// float temp_vertex_0[3] = {0.0,0.0,0.0};
					temp_vertex_0[0] = new_arr[ 3*temp_element[0]  ];
					temp_vertex_0[1] = new_arr[ 3*temp_element[0] +1 ];
					temp_vertex_0[2] = new_arr[ 3*temp_element[0] +2 ];

				// float temp_vertex_1[3] = {0.0,0.0,0.0};
					temp_vertex_1[0] = new_arr[ 3*temp_element[1] ];
					temp_vertex_1[1] = new_arr[ 3*temp_element[1] +1 ];
					temp_vertex_1[2] = new_arr[ 3*temp_element[1] +2 ];
				
				// float temp_vertex_2[3] = {0.0,0.0,0.0};
					temp_vertex_2[0] = new_arr[ 3*temp_element[2] ];
					temp_vertex_2[1] = new_arr[ 3*temp_element[2] +1 ];
					temp_vertex_2[2] = new_arr[ 3*temp_element[2] +2 ];
				
				// float temp_vertex_3[3] = {0.0,0.0,0.0};
					temp_vertex_3[0] = new_arr[ 3*temp_element[3] ];
					temp_vertex_3[1] = new_arr[ 3*temp_element[3] +1 ];
					temp_vertex_3[2] = new_arr[ 3*temp_element[3] +2 ];

				is_inside = computeBaricentricCoordinates(temp_position, temp_vertex_0, temp_vertex_1, temp_vertex_2, temp_vertex_3, baricentricCoordinates);

				if( is_inside ){
					if(final_image[i] == 0 ) 
					{
						// Primer código !
	//					 final_image[i] = element_number+1;
						// Segundo código !
				/*	*/		// 1. Hay que acalcular primero la posición cartesiana en el espacio de la imagen de entrada
					//	std::cout << "Temp. Element: [" << temp_element[0] << ", " << temp_element[1] << ", " << temp_element[2] << ", " << temp_element[3] << "] " << std::endl;
					
							old_vertex_0[0] = arr[ 3*temp_element[0]  ];
							old_vertex_0[1] = arr[ 3*temp_element[0] +1 ];
							old_vertex_0[2] = arr[ 3*temp_element[0] +2 ];
		
							old_vertex_1[0] = arr[ 3*temp_element[1] ];
							old_vertex_1[1] = arr[ 3*temp_element[1] +1 ];
							old_vertex_1[2] = arr[ 3*temp_element[1] +2 ];
					
							old_vertex_2[0] = arr[ 3*temp_element[2] ];
							old_vertex_2[1] = arr[ 3*temp_element[2] +1 ];
							old_vertex_2[2] = arr[ 3*temp_element[2] +2 ];
									
							old_vertex_3[0] = arr[ 3*temp_element[3] ];
							old_vertex_3[1] = arr[ 3*temp_element[3] +1 ];
							old_vertex_3[2] = arr[ 3*temp_element[3] +2 ];

						computeCartessianCoordinates( baricentricCoordinates, old_vertex_0, old_vertex_1, old_vertex_2, old_vertex_3, sss_position);

							// 2. Localizamos el pixel continuo e interpolamos el valor. 
					
						//	pixel[0] = (sss_position[0]-origenCuda[0])/spacingCuda[0];
						//	pixel[1] = (sss_position[1]-origenCuda[1])/spacingCuda[1];
						//	pixel[2] = (sss_position[2]-origenCuda[2])/spacingCuda[2];
						

							pixel[0] = (sss_position[0]-origen_input[0])/spacing_input[0];
							pixel[1] = (sss_position[1]-origen_input[1])/spacing_input[1];
							pixel[2] = (sss_position[2]-origen_input[2])/spacing_input[2];

													//Checkeando la información ! 
				/*		std::cout << "Temporal position. [" << temp_position[0] <<", " << temp_position[1] << ", " << temp_position[2] << "] " << std::endl;
						std::cout << "    coordenadas baricentricas: [" << baricentricCoordinates[0] << ", " << baricentricCoordinates[1] << ", " << baricentricCoordinates[2] << ", " << baricentricCoordinates[3] << "] " << std::endl;
						std::cout << "Reconstructed position. [" << sss_position[0] <<", " << sss_position[1] << ", " << sss_position[2] << "] " << std::endl;
						std::cout << " El pixel es [" << pixel[0] <<", " << pixel[1] << ", " << pixel[2] << "]  de size = [" << sizeCuda[0] <<", " << sizeCuda[1] << ", "<< sizeCuda[2] << "] " << std::endl;
						std::cout << std::endl;
						Sleep(250);*/

						index = (unsigned int) floor(pixel[0]+0.5) + (sizeCuda[0] * (unsigned int) floor(pixel[1]+0.5) ) + (sizeCuda[0]*sizeCuda[1]*(unsigned int) floor(pixel[2]+0.5));
						if( (index<(sizeCuda[0]*sizeCuda[1]*sizeCuda[2])) && (0<index))
							final_image[i] = image_input[index];
						else{
							std::cout << "Fuera de image!" << std::endl;
						//Checkeando la información ! 
						std::cout << "Temporal position. [" << temp_position[0] <<", " << temp_position[1] << ", " << temp_position[2] << "] " << std::endl;
						std::cout << "    coordenadas baricentricas: [" << baricentricCoordinates[0] << ", " << baricentricCoordinates[1] << ", " << baricentricCoordinates[2] << ", " << baricentricCoordinates[3] << "] " << std::endl;
						std::cout << "Reconstructed position. [" << sss_position[0] <<", " << sss_position[1] << ", " << sss_position[2] << "] " << std::endl;
						std::cout << "Index: " << index << " de " << (sizeCuda[0]*sizeCuda[1]*sizeCuda[2]) << std::endl;
						std::cout << " El pixel es [" << pixel[0] <<", " << pixel[1] << ", " << pixel[2] << "]  de size = [" << sizeCuda[0] <<", " << sizeCuda[1] << ", "<< sizeCuda[2] << "] " << std::endl;
						std::cout << std::endl;
						Sleep(250);
						}

						
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


					} // endif( final_image[i]==0 )
				} // endif( is_inside == true )		
			} // endfor
			} // endif( number_of_elements_here!=0 )
			} else { final_image[i] = 0; };
	}
	
	std::cout << "LLega hasta aquí?? " << std::endl;

/* Imagen ITK */

	start[0] = 0; start[1] = 0; start[2] = 0;
	size_im[0] = im_size[0]; size_im[1] = im_size[1]; size_im[2] = im_size[2];
	ImageType::RegionType region(start, size_im);
	
	origen_im[0] = im_origen[0]; origen_im[1] = im_origen[1]; origen_im[2] = im_origen[2];
	spacing_im[0] = im_spacing[0]; spacing_im[1] = im_spacing[1]; spacing_im[2] = im_spacing[2];
	
	//ImportFilterType::Pointer importfilter = ImportFilterType::New();
		importfilter->SetOrigin( origen_im );
		importfilter->SetRegion( region );
		importfilter->SetSpacing( spacing_im );
		importfilter->SetDirection( direction );
	//const bool importImagefilterWillOwnTheBuffer = false;
		//importfilter->SetImportPointer( flags, numberOfPixels, importImagefilterWillOwnTheBuffer );
		//importfilter->SetImportPointer( final_image, im_numberOfPixels, importImagefilterWillOwnTheBuffer );
		importfilter->SetImportPointer( final_image, im_numberOfPixels, true );

		//importfilter->SetImportPointer( final_image, numberOfPixels, importImagefilterWillOwnTheBuffer );
		importfilter->Update();

	compressed_image = importfilter->GetOutput();
	std::cout << std::endl;
	std::cout << " Sale de la Compresion !! " << std::endl;
	std::cout << std::endl;

	delete[] correspondingElement;
	delete[] cumsum;
	//delete[] flags;
	//delete[] temp_flags;
	delete[] bb;
	delete[] cha;
	delete[] arr;
	delete[] new_arr;

//	delete[] origenCuda;
//	delete[] spacingCuda;
//	delete[] sizeCuda;
//	delete[] image_input[];
}

#endif