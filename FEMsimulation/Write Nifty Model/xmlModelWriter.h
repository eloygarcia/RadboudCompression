#ifndef __xmlModelWriter_h
#define __xmlModelWriter_h

#pragma once

//#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <sstream>
#include <iterator>

#include <list>

#include "xml.h"

class xmlModelWriter
{
public:
	xmlModelWriter(void);
	~xmlModelWriter(void);

private:
	XML* document;
	XMLElement* root;
	int count;
	int children_param;

	std::ostringstream tr_oss;
	std::ostringstream rot_oss;

	stringstream ss; 
	std::ostringstream oss;

	// stringstream ss; 
	// std::ostringstream oss;

	// stringstream ss; 
	// std::ostringstream oss;
	stringstream p_ss;
	std::ostringstream p_oss;

	std::ostringstream a_ss;
	std::ostringstream b_ss;
	std::ostringstream c_ss;
	std::ostringstream disp_ss;
	stringstream n_slaves;
	std::ostringstream slvs_ss;

public:
	void addVTKMesh(char* inputfilename, char* Type, std::vector<double> Translation, std::vector<double> Rotation, char* scaleFactor);
	void addNodes( std::vector<double> nodes, char* Dof );
	void addElements( std::vector<int> elements, char* Type );
	void addElementSet( std::vector<int> elementSet, char* Type, std::vector<double> elasticParams);
	void addConstraint(std::vector<int> nodes, char* Type, char* Dof);
	void addSystemParams(char* parameter, char* value);
	void addContactPlate(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> disp, std::list<int> slvs);

	void addGravity(std::vector<int> nodes, char* AccelerationMagnitude, std::vector<double>  AccelerationDirection);
	void addOutput(char* output, char* freq);

	void writeModel( char* outputfilename );

};

//#include "xmlModelWriter.cpp"
#endif