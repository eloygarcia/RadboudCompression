#pragma once

#ifndef __niftysimEjecutable_h
#define __niftysimEjecutable_h

#include "auxiliarDefinitions.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"

#include "vtkDataArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCell.h"

#include "vtkCommand.h"

#include "xmlModelWriter.h"

class ErrorObserver : public vtkCommand
{
public:
  ErrorObserver():
    Error(false),
    Warning(false),
    ErrorMessage(""),
    WarningMessage("") {}
  static ErrorObserver *New()
  {
  return new ErrorObserver;
  }
  bool GetError() const
  {
  return this->Error;
  }
  bool GetWarning() const
  {
  return this->Warning;
  }
  void Clear()
  {
  this->Error = false;
  this->Warning = false;
  this->ErrorMessage = "";
  this->WarningMessage = "";
  }
  virtual void Execute(vtkObject *vtkNotUsed(caller),
                       unsigned long event,
                       void *calldata)
  {
  switch(event)
    {
    case vtkCommand::ErrorEvent:
      ErrorMessage = static_cast<char *>(calldata);
      this->Error = true;
      break;
    case vtkCommand::WarningEvent:
      WarningMessage = static_cast<char *>(calldata);
      this->Warning = true;
      break;
    }
  }
  std::string GetErrorMessage()
  {
  return ErrorMessage;
  }
std::string GetWarningMessage()
  {
  return WarningMessage;
  }
private:
  bool        Error;
  bool        Warning;
  std::string ErrorMessage;
  std::string WarningMessage;
};

class NiftySimEjecutable
{
public:
	NiftySimEjecutable(void);
	~NiftySimEjecutable(void);
//
private:
	//xmlModelWriter * new_model;
	int sys_a;
	float de;
		
	std::vector<double> origen;
	std::vector<double> spacing;
	std::vector<int> size;

	// ==== De: Compression Function =======
	std::string niftysim_fn;
	std::string princ;
	std::string niftysim_filename;
	std::string action;

	std::vector<double> temp_mat;
	std::list<int> unique_node_vector;

	std::vector<int> unique_element_vector;
	std::vector<int> temp_vect;

	std::list<int> u_node_vector;
	std::vector<double> bb;

	// ========= De:   ====================
	std::vector<double> m_i_points;
	std::string compressedgridname;
	std::string vtk_niftysim_filename;

	modelParameters m_parameters;
	std::string m_outputDir;
	std::vector<double> m_finalpoints;

	double * pt;
	double * d;
	double * temp;
	double * auxPoint;

	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> c;
	std::vector<double> disp;

		vtkSmartPointer<vtkUnstructuredGridReader> niftysimOutput_reader;
		vtkSmartPointer<vtkUnstructuredGrid> niftygrid;
		vtkSmartPointer<vtkUnstructuredGridWriter> UGridWriter;
		vtkSmartPointer<vtkUnstructuredGrid> comprrGrid;
		vtkPoints* initial_points;
		vtkDataArray* vtk_disp;
		vtkCellArray* cells;
	
	void CompressionFunction(modelParameters myParameters, std::string output_dir);
	void getPoints_compressedBreast( vtkSmartPointer<vtkUnstructuredGrid> niftysim_mesh,
									std::vector<double> &f_points);
	void createUnstructuredGrid( std::vector<double> point_vector, 
								 vtkCellArray* cells, 
								 vtkSmartPointer<vtkUnstructuredGrid> &grid);

public:
	void SetModelParameters( modelParameters myParameters ) { m_parameters = myParameters; };
	void SetOutputDir( std::string output_dir ){ m_outputDir = output_dir;};
	
	int Update();

	std::vector<double> GetFinalPoints(){ return m_finalpoints;};
	
};


//#include "NiftySimEjecutable.cpp"

#endif