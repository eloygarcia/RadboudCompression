#ifndef __auxiliarFunctions_Transformations_h
#define __auxiliarFunctions_Transformations_h

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <iostream> 
#include <vector>
#include <math.h>
#include <cmath>

extern float grad2rad( float degree)
{
	float aa =((degree*(float)M_PI)/180);
	return aa;
}

extern std::vector<double> get_boundingBox(std::vector<double> points)
{
	double x_max = -9999999; double x_min = 99999999;
	double y_max = -9999999; double y_min = 99999999;
	double z_max = -9999999; double z_min = 99999999;

	for(unsigned int i=0; i<(points.size()/3); i++)
	{
		if(points[3*i] < x_min) x_min= points[3*i];
		if(points[3*i] > x_max) x_max= points[3*i];

		if(points[3*i+1] < y_min) y_min= points[3*i+1];
		if(points[3*i+1] > y_max) y_max= points[3*i+1];

		if(points[3*i+2] < z_min) z_min= points[3*i+2];
		if(points[3*i+2] > z_max) z_max= points[3*i+2];
	}

	std::vector<double> temp;
		temp.push_back(x_min); 
		temp.push_back(x_max);
		temp.push_back(y_min);
		temp.push_back(y_max);
		temp.push_back(z_min);
		temp.push_back(z_max);

	std::cout << "BoundingBox: ["<< temp[0]<< ", "<< temp[1] << ", "<< temp[2] <<
		 ", "<< temp[3] << ", "<< temp[4] << ", "<< temp[5] << std::endl;
	return temp;
}

extern std::vector<float> get_centerofmass(std::vector<double> points)
{
	float sum_x = 0.0;
	float sum_y = 0.0;
	float sum_z = 0.0;
	int n_points = points.size()/3;

	for(int i=0; i<(n_points); i++)
	{
		sum_x += (float) points[3*i];
		sum_y += (float) points[3*i+1];
		sum_z += (float) points[3*i+2];
	}

	std::vector<float> centerofmass;
	centerofmass.push_back( sum_x/n_points );
	centerofmass.push_back( sum_y/n_points );
	centerofmass.push_back( sum_z/n_points );

	std::cout << std::endl;
	std::cout << "Center of mass: " << centerofmass[0] << ", " << 
		centerofmass[1] << ", " << centerofmass[2] << "]" << std::endl;

	return centerofmass;
}

// ============  Translation & Rotation =======================
extern void set_translation(std::vector<double> i_points,  std::vector<double> translation, std::vector<double> &f_points)
{
		for(unsigned int i=0; i<(i_points.size()/3); i++)
		{
			f_points[3*i]= i_points[3*i] + translation[0];
			f_points[3*i+1] = i_points[3*i+1] + translation[1];
			f_points[3*i+2] = i_points[3*i+2] + translation[2];
		}
}

extern void set_rotation(std::vector<double> i_points, std::vector<double> &f_points, std::vector<float> center, float RotX, float RotY, float RotZ)
{
	double temp_point[3] = {0.0, 0.0, 0.0};
	double temp_rot[3] = {0.0, 0.0, 0.0};
	std::vector<double> temp_vector;
	for(unsigned int i=0; i<(i_points.size() /3); i++)
	{
		temp_point[0] = i_points[3*i];
		temp_point[1] = i_points[3*i+1];
		temp_point[2] = i_points[3*i+2];

		// ================== Rotation ====================
		// --------- Tranlation: center of mass to 0,0,0 ----------
		temp_point[0] -= center[0];
		temp_point[1] -= center[1];
		temp_point[2] -= center[2];
		
		// ------------------- Rotation ------------------
			if(RotX != 0)
			{
				// RotX -- radians
				temp_rot[0] = temp_point[0];
				temp_rot[1] = temp_point[1] * cos(RotX) - temp_point[2] * sin(RotX);
				temp_rot[2] = temp_point[1] * sin(RotX) + temp_point[2] * cos(RotX);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotY != 0)
			{
				// RotY -- radians
				temp_rot[0] = temp_point[0] * cos(RotY) + temp_point[2] * sin(RotY);
				temp_rot[1] = temp_point[1];
				temp_rot[2] = - temp_point[0] * sin(RotY) + temp_point[2] * cos(RotY);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotZ != 0)
			{
				// RotZ -- radians
				temp_rot[0] = temp_point[0] * cos(RotZ) - temp_point[1] * sin(RotZ);
				temp_rot[1] = temp_point[0] * sin(RotZ) + temp_point[1] * cos(RotZ);
				temp_rot[2] = temp_point[2];
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
		// ------------ Back ---------------
		temp_point[0] += center[0];
		temp_point[1] += center[1];
		temp_point[2] += center[2];


		// ========== Compose vector ============
		temp_vector.push_back( temp_point[0] );
		temp_vector.push_back( temp_point[1] );
		temp_vector.push_back( temp_point[2] );
	}

	f_points.clear();
	f_points = temp_vector;
}

extern void set_transform(std::vector<double> i_points, std::vector<double> &f_points, std::vector<double> translation, std::vector<float> center, float RotX, float RotY, float RotZ)
{
	f_points = i_points;
/*
	ROTACION!
*/
	if( RotX!=0 || RotY!=0 ||  RotZ!=0) {
	double temp_point[3] = {0.0, 0.0, 0.0};
	double temp_rot[3] = {0.0, 0.0, 0.0};
	std::vector<double> temp_vector;
	for(unsigned int i=0; i<(f_points.size() /3); i++)
	{
		temp_point[0] = f_points[3*i];
		temp_point[1] = f_points[3*i+1];
		temp_point[2] = f_points[3*i+2];

		// ================== Rotation ====================
		// --------- Translation: center of mass to 0,0,0 ----------
		temp_point[0] = temp_point[0] - center[0];
		temp_point[1] = temp_point[1] - center[1];
		temp_point[2] = temp_point[2] - center[2];
		
		// ------------------- Rotation ------------------
			if(RotX != 0)
			{
				// RotX -- radians
				temp_rot[0] = temp_point[0];
				temp_rot[1] = temp_point[1] * cos(RotX) - temp_point[2] * sin(RotX);
				temp_rot[2] = temp_point[1] * sin(RotX) + temp_point[2] * cos(RotX);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotY != 0)
			{
				// RotY -- radians
				temp_rot[0] = temp_point[0] * cos(RotY) + temp_point[2] * sin(RotY);
				temp_rot[1] = temp_point[1];
				temp_rot[2] = - temp_point[0] * sin(RotY) + temp_point[2] * cos(RotY);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotZ != 0)
			{
				// RotZ -- radians
				temp_rot[0] = temp_point[0] * cos(RotZ) - temp_point[1] * sin(RotZ);
				temp_rot[1] = temp_point[0] * sin(RotZ) + temp_point[1] * cos(RotZ);
				temp_rot[2] = temp_point[2];
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
		// ------------ Back ---------------
		temp_point[0] = temp_point[0] + center[0];
		temp_point[1] = temp_point[1] + center[1];
		temp_point[2] = temp_point[2] + center[2];


		// ========== Compose vector ============
		temp_vector.push_back( temp_point[0] );
		temp_vector.push_back( temp_point[1] );
		temp_vector.push_back( temp_point[2] );
	}

	f_points.clear();
	f_points = temp_vector;
	}

/*
	TRANSLACION!!
*/
		for(unsigned int i=0; i<(f_points.size()/3); i++)
		{
			f_points[3*i]= f_points[3*i] + translation[0];
			f_points[3*i+1] = f_points[3*i+1] + translation[1];
			f_points[3*i+2] = f_points[3*i+2] + translation[2];
		}
}

#endif