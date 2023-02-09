#pragma once

#ifndef __NewProjection_h
#define __NewProjection_h

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>

#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkImageFileWriter.h"

//typedef unsigned short pixelType;
//typedef short pixelType;
typedef char pixelType;

typedef itk::Image<pixelType, 3> ImageType;
typedef itk::ImportImageFilter<pixelType,3> ImportFilterType;
typedef itk::ImageFileWriter< ImageType > WriterType;

class NewProjection
{
public:
	NewProjection(void);
	~NewProjection(void);

private:
	pixelType * image_input;
	pixelType * final_image;

	float * bounding_box;
	int* correspondingElement;

	int n;
	int m;

	float * m_origen;
	float * m_spacing;
	unsigned int * m_size;
	int numberOfPixels;

	float * origenCuda;
	unsigned int * sizeCuda;
	float * spacingCuda;

	float * temp_x;
	float * temp_y;
	float * temp_z;
	
	int sp_x_min;
	int sp_y_min;
	int sp_z_min;
	
	int sp_x_max;
	int sp_y_max;
	int sp_z_max;

	int idx;
	int nim;

	//int sum;
	long maximum_elements;

	int index_z;
	int index_y;
	int index_x;

	float * temp_position;
	float * baricentricCoordinates;
	float * barcor;

	int index_start;
	int number_of_elements_here;
	int count;

	int element_index_number;
	int element_number;

	int * temp_element;
	float * temp_vertex_0;
	float * temp_vertex_1;
	float * temp_vertex_2;
	float * temp_vertex_3;

//	bool is_inside;
	float * voxel_position;

	float *  old_vertex_0;
	float * old_vertex_1;
	float * old_vertex_2;
	float * old_vertex_3;
	float * sss_position;
	pixelType value;
	float * pixel;


	ImageType::IndexType start;
	ImageType::SizeType size_im;
	ImageType::PointType origen_im;
	ImageType::SpacingType spacing_im;
	ImageType::DirectionType direction2d;

	//int* corresponding_element;
	//int * flags;
	//int* temp_flags;
	//float* bb;
	//int* cumsum;

	ImportFilterType::Pointer importfilter;
	
	// function Determinante
	float det;
	// function TrilinearInterpolation
	int x0, y0, z0, x1, y1, z1;
	float xd;
	float yd;
	float zd;
	int idx_1, idx_2;
	unsigned int index;
	pixelType valido;
	int * storing;

	pixelType c00, c01, c10, c11;
	pixelType c0, c1;
	pixelType c;

	// function nearestNeihborInterpolation
	pixelType val000, val001, val010, val011, val100, val101, val110, val111;
	pixelType sumVal;
	// function ComputeBaryCenticCoordinates
	float V;
	float v1, v2, v3, v4;
	bool is_inside;
	// function max
	float maximum;
	// function min
	float minimum;
	//function get_boundingbox
	double x_max, x_min, y_max, y_min, z_max, z_min;

	int best;
	int temp_value;

// 
private:
	float Determinante( float* v0, float* v1, float* v2, float* v3 );
	void trilinearinterpolation(const pixelType* imagepointer, const int* size, const float* pixel, pixelType & value); // es un pixel continuo, con lo cual, vale
	bool computeBaricentricCoordinates(float* posicionPunto, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* baricentricCoordinates);
	void computeCartessianCoordinates(float* baricentricCoordinates, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* posicionPunto);
	float max(float * idx);
	float min(float * idx);
	void get_boundingBox(float* points, int numberOfPoints, float* bounding_box);

	void nearestNeighborInterpolation( const pixelType * imagepointer, const int* s_size, const float* pixel, pixelType & value);

	//ImageType::Pointer compressed_image;

public:
	void getCompressedImage(std::vector<double> initial_points, 
		std::vector<double> deformed_points, 
		std::vector<int> elements,
		ImageType::Pointer image_original, 
		ImageType::Pointer &compressed_image, float* im_spacing);

};

// #include "NewProjection.cpp"

#endif