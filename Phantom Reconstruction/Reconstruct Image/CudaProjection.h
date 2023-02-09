#pragma once

#ifndef __CudaProjection_h
#define __CudaProjection_h

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
/*
#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkImageFileWriter.h"

#include <unistd.h>

//typedef unsigned short pixelType;
//typedef short pixelType;
typedef char pixelType;

typedef itk::Image<pixelType, 3> ImageType;
typedef itk::ImportImageFilter<pixelType,3> ImportFilterType;
typedef itk::ImageFileWriter< ImageType > WriterType;
*/
class CudaProjection
{
public:
	CudaProjection(void);
	~CudaProjection(void);

private:
    cudaError_t cudaStatus;
	unsigned int bl; // Numero de bloques a lanzar, en conjuntos de 512 threads

// HOST
    // initial image
	int * m_initial_size;
	float * m_initial_spacing;
	float * m_initial_origen;
    unsigned char *	m_initial_imagepointer;

    // Mesh 
    float * m_i_points;
    float * m_f_points;
    int * m_elements;
    
    // Grid!!
	float * m_grid_origen;
	float * m_grid_spacing;
	int  * m_grid_size;

    int * m_flags;
	int * m_cumsum;
	int * m_correspondingElements;

    // === Imagen Simulada
	int * m_simulada_size;
	float * m_simulada_spacing;
	float * m_simulada_origen;
	unsigned char * m_simulada_imagepointer;

// DEVICE 
	// === Imagen Inicial
	int * dev_initial_size;
	float * dev_initial_spacing;
	float * dev_initial_origen;
long numberOfPixels; // esta es host !
	unsigned char * dev_initial_imagepointer;
	
	// === Mesh and Grid ! ===
int numberOfPoints;
	float * dev_i_points;
	float * dev_f_points;
int numberOfElements;
	int * dev_elements;

	float * dev_grid_origen;
	float * dev_grid_spacing;
	int * dev_grid_size;
int numberOfVoxelsGrid;

	int * dev_flags;
	int * dev_cumsum;
long maximumCorrespondingElements;
	int * dev_correspondingElements;
	// =======================

	// === Imagen Simulada
	int* dev_simulada_size;
	float* dev_simulada_spacing;
	float* dev_simulada_origen;
long numberOfPixels_Simulada; // esta tambi�n es host !
	unsigned char * dev_simulada_imagepointer;

	
// METODOS !
private:
    void Initialize();

public:
    void SetInitialImage(int* imageSize, 
						 float* imageSpacing, 
						 float* imageOrigin, 
						 unsigned char* imagePointer);
    void SetMesh(int numPoints, float* i_points, float* f_points, int numElements, int* elements);
    void SetGrid(int* gridSize, float* gridSpacing, float* gridOrigin,
                int* flags, int* cumsum, int* correspondingElements);
    void SetFinalImage(int* finalSize, float* finalSpacing, float* finalOrigin, unsigned char* finalPointer);

public: 
	void Update();
	unsigned char* GetFinalImage() { return m_simulada_imagepointer;};



};

// #include "NewProjection.cpp"

#endif