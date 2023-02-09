//#pragma once

#ifndef __RegularGrid_h
#define __RegularGrid_h

#include <stdio.h>
#include <math.h>
#include <iostream>

class RegularGrid
{
public:

	RegularGrid(void);
	~RegularGrid(void);

private:

	float * m_origen;
	bool is_origen;
	
	float * m_spacing;
	bool is_spacing;

	int * m_size;
	bool is_size;
	int m_numberOfPixels;

	bool is_init;

	float * m_points;
	int m_numberOfPoints;
	bool is_points;

	float * m_boundingBox;
	bool is_boundingBox;

	int m_numberOfElements;
	bool is_elements;
	int * m_elements;

	float * m_OBBs;
	
	// 
	int * m_flags;
	int * m_cumsum;
	int * m_correspondingElement;

	bool is_update;
public:

	void Initialize();
	void ComputeBoundingBox();
	float min( float * idx );
	float max( float * idx );

	void setPoints( int numofpoints, float* points ){
		m_numberOfPoints = numofpoints;
		m_points = points;

		is_points = true;
	};

	void setElements(int numofelem, int* elements) {
		m_elements = elements;
		m_numberOfElements = numofelem;
		is_elements = true;

		m_OBBs = new float[ 6*m_numberOfElements ];
	};

	void setSpacing(float* spacing){
		m_spacing = spacing;
		is_spacing = true;
	};

	void setOrigen(float* origen){
		m_origen = origen;
		is_origen = true;
	};

	void setSize(int * size){
		m_size = size;
		is_size = true;
		
		m_numberOfPixels = size[0] * size[1] *  size[2];
	};

	void Update();

	int * getFlags(){ return m_flags;};
	int * getCumSum(){ return m_cumsum;};
	int * getElementList(){ return m_correspondingElement;};
	int * getSize(){return m_size;};
	float * getSpacing(){ return m_spacing;};
	float * getOrigen(){ return m_origen;};

	void returnGrid( int * flags, int * cumsum, int* correspElements){
		if(is_update){		
			flags = m_flags;
			cumsum = m_cumsum;
			correspElements = m_correspondingElement;
		} else{ 
			flags = 0;
			cumsum = 0;
			correspElements = 0;
		}
	}

};

#include "RegularGrid.cpp"
#endif