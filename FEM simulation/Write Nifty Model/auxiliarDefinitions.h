#ifndef __auxiliarDefinitions_h
#define __auxiliarDefinitions_h

#include <iostream>
#include <vector>
#include <list>

struct modelParameters
{
	modelParameters()
	{
		angle = 0;

		number_of_tissues = 0;
		// Nodes
		nodes_Dof = "3";
		// Elements
		type = "T4";
		//subsets_on = 0;
		
		// Set Material;
		material = "NH";
		int numParams = 2;

		material_parameters.push_back( 1491.64 );
		material_parameters.push_back( 147672.24 ); // Fatty

		material_parameters.push_back( 5050.167 );
		material_parameters.push_back( 499966.555 ); // Glandular

		material_parameters.push_back( 5050.167 );
		material_parameters.push_back( 499966.555 ); // Glandular

		material_parameters.push_back( 5050.167 );
		material_parameters.push_back( 499966.555 ); // Glandular

		// SystemParams;
		// timeStep = "0.0001";
		 timeStep = "0.00005";
		// timeStep = "0.00001";
		totalTime = "3";
		DampingCoeff = "0.065";
		Density = "1";
		DoDeformableCollision = "0";

		// Constraint;
		constraint_type = "Fix";
		constraint_axis = "0";

		// Compression;
		thickness = 0;

		// Pectoral Muscle
		pectoralAngle = 0.0;

		// Paddle; 
		paddleAngle = 0.0;

		// Gravity
		gravity = 100;
	}

	std::string name;

	int number_of_tissues;

	char* inputfilename;
	// Nodes
	std::vector<double> nodes;
	char* nodes_Dof;
	std::vector<double> boundingBox;

	// Elements
	std::vector<int> elements;
	char* type;
	std::vector<int> elementTissue;

	// Subsets
	int subsets_on;

	// System Params;
	char* timeStep;
	char* totalTime;
	char* DampingCoeff;
	char* Density;
	char* DoDeformableCollision;
	
	// Set Material;
	char* material;
	int numParams;
	std::vector<double> material_parameters; 

	// Constraint;
	std::vector<int> constraint_bC_list;
	char* constraint_type;
	char* constraint_axis;
	// ymin --> normal positiva
	// ymax --> normal negativa

	// Compression;
	std::string position;
	
	// mri data !

	// mammogram data !
	std::vector<float> cen;
	float angle;
	
	std::vector<double> mm_origen;
	std::vector<double> mm_spacing;
	std::vector<int> mm_size;

	std::string side;
	float thickness;

	double pectoralAngle;

	float paddleAngle;

	float gravity;
}; 

#endif