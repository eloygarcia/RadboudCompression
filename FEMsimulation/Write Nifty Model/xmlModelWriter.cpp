#ifndef __xmlModelWriter_cpp
#define __xmlModelWriter_cpp

#include "xmlModelWriter.h"

xmlModelWriter::xmlModelWriter(void)
{	
	count=1;
	document = new XML;
	XML* x=0;

	// create root node.
	root = new XMLElement(0,"<Model>");
	document->SetRootElement(root);

	// create system params node.
	root->AddElement("<SystemParams>");
	children_param = 0;
	
}

xmlModelWriter::~xmlModelWriter(void)
{
	delete document;
}

// *****************************************************

void xmlModelWriter::addVTKMesh(char* inputfilename, char* Type, std::vector<double> Translation, std::vector<double> Rotation, char* scaleFactor)
{
	//
	// std::ostringstream tr_oss;	
	std::copy(Translation.begin(), Translation.end(), std::ostream_iterator<double>(tr_oss, " "));
	// std::ostringstream rot_oss;
	std::copy(Rotation.begin(), Rotation.end(), std::ostream_iterator<double>(rot_oss, " "));
	//
	this->root->AddElement("<VTKMesh>");
	XMLElement* VtkMesh = this->document->GetRootElement()->GetChildren()[ this->count];
		VtkMesh->AddVariable( new XMLVariable("Type", Type) );
		VtkMesh->AddContent( inputfilename,0 );
		VtkMesh->AddElement("<Translation>");
			VtkMesh->GetChildren()[0]->AddContent( (char*)tr_oss.str().c_str(), 0 );
		VtkMesh->AddElement("<Rotation>");
			VtkMesh->GetChildren()[1]->AddContent( (char*)rot_oss.str().c_str(), 0 );
		VtkMesh->AddElement("<ScaleFactor>");
			VtkMesh->GetChildren()[2]->AddContent( scaleFactor , 0 );
	this->count += 1;

	// == clear ==
	tr_oss.str("");
	rot_oss.str("");
}

void xmlModelWriter::addNodes(std::vector<double> nodes, char* Dof)
{
	// stringstream ss; 
	ss << (nodes.size()/3);
	// std::ostringstream oss;
	std::copy(nodes.begin(), nodes.end(), std::ostream_iterator<double>(oss, " "));

	this->root->AddElement("<Nodes>");
	XMLElement* NodesElement = this->document->GetRootElement()->GetChildren()[ this->count ];
		NodesElement->AddVariable( new XMLVariable("DOF", Dof ) );
		NodesElement->AddVariable( new XMLVariable("NumNodes", (char*)ss.str().c_str() ) );
		NodesElement->AddContent( (char*)(oss.str().c_str() ),0); 
		this->count +=1;

	// == clear ==
	ss.str(""); //ss.clear();
	oss.str(""); //oss.clear();
}

void xmlModelWriter::addElements(std::vector<int> elements, char* Type)
{
	// stringstream ss; 
	ss << (elements.size()/4);

	// std::ostringstream oss;
	std::copy(elements.begin(), elements.end(), std::ostream_iterator<double>(oss, " "));

	this->root->AddElement("<Elements>");
	XMLElement* ElementNode = this->document->GetRootElement()->GetChildren()[ this->count ];
		ElementNode->AddVariable( new XMLVariable("NumEls", (char*)(ss.str().c_str())) );
		ElementNode->AddVariable( new XMLVariable("Type", Type) );
		ElementNode->AddContent( (char*)(oss.str().c_str()) , 0);
		this->count += 1;

	// == clear ==
	ss.str(""); //ss.clear();
	oss.str(""); //oss.clear();
}

void xmlModelWriter::addElementSet(std::vector<int> elementSet, char* Type, std::vector<double> elasticParams)
{ 
	// Element Set
	// stringstream ss; 
	ss << elementSet.size();
	// std::ostringstream oss;
	std::copy(elementSet.begin(), elementSet.end(), std::ostream_iterator<double>(oss, " "));

	// Elastic Params
	// stringstream p_ss;
	p_ss << elasticParams.size();
	// std::ostringstream p_oss;
	std::copy(elasticParams.begin(), elasticParams.end(), std::ostream_iterator<double>(p_oss, " "));
	
	// Dado que s�lo vamos a tener un tipo de material, se puede quedar. En otros casos NO!
	// Creo que si se puede!! Pero hay que preparar los elementos aparte

	this->root->AddElement("<ElementSet>");
	XMLElement* elementSetNode = this->document->GetRootElement()->GetChildren()[ this->count ];
		elementSetNode->AddVariable( new XMLVariable( "Size", (char*)ss.str().c_str() ) );
		elementSetNode->AddElement("<Material>");
				elementSetNode->GetChildren()[0]->AddVariable( new XMLVariable("Type", Type));
				elementSetNode->GetChildren()[0]->AddElement("<ElasticParams>");
					elementSetNode->GetChildren()[0]->GetChildren()[0]->AddVariable("NumParams", (char*)p_ss.str().c_str() );
					elementSetNode->GetChildren()[0]->GetChildren()[0]->AddContent((char*)p_oss.str().c_str(), 0 );
		elementSetNode->AddContent((char*)oss.str().c_str(),0);
		this->count +=1;

	// == clear ==
	ss.str(""); //ss.clear();
	oss.str(""); //oss.clear();
	
	p_ss.str("");
	p_oss.str("");
}

void xmlModelWriter::addConstraint(std::vector<int> nodes, char* Type, char* Dof)
{
	// stringstream n_ss; 
	ss << nodes.size();
	// std::ostringstream nodes_ss;
	std::copy(nodes.begin(), nodes.end(), std::ostream_iterator<int>(oss, " "));
	//
	this->root->AddElement("<Constraint>");
	XMLElement* Constraint = this->document->GetRootElement()->GetChildren()[ this->count ];
		Constraint->AddVariable( new XMLVariable("DOF", Dof));
		Constraint->AddVariable( new XMLVariable("NumNodes", (char*)ss.str().c_str() ));
		Constraint->AddVariable( new XMLVariable("Type", Type));
		Constraint->AddElement("<Nodes>");
			Constraint->GetChildren()[0]->AddContent( (char*)oss.str().c_str() ,0);
	this->count += 1;

	// == clear ==
	ss.str(""); //ss.clear();
	oss.str(""); //oss.clear();
}

void xmlModelWriter::addContactPlate(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> disp, std::list<int> slvs)
{
	// vector to char array
	// std::ostringstream a_ss;
	std::copy(a.begin(), a.end(), std::ostream_iterator<double>(a_ss, " "));
	// std::ostringstream b_ss;
	std::copy(b.begin(), b.end(), std::ostream_iterator<double>(b_ss, " "));
	// std::ostringstream c_ss;
	std::copy(c.begin(), c.end(), std::ostream_iterator<double>(c_ss, " "));
	// std::ostringstream disp_ss;
	std::copy(disp.begin(), disp.end(), std::ostream_iterator<double>(disp_ss, " "));
	// stringstream n_slaves;
	n_slaves << slvs.size();
	// std::ostringstream slvs_ss;
	std::copy(slvs.begin(), slvs.end(), std::ostream_iterator<int>(slvs_ss, " "));
	//
	this->root->AddElement("<ContactPlate>");
	XMLElement* contactPlate = this->document->GetRootElement()->GetChildren()[ this->count ];
		contactPlate->AddElement("<a>");
		contactPlate->GetChildren()[0]->AddContent( (char*)a_ss.str().c_str(),0);
		contactPlate->AddElement("<b>");
		contactPlate->GetChildren()[1]->AddContent( (char*)b_ss.str().c_str(),0);
		contactPlate->AddElement("<c>");
		contactPlate->GetChildren()[2]->AddContent( (char*)c_ss.str().c_str(),0);
		contactPlate->AddElement("<Disp>");
		contactPlate->GetChildren()[3]->AddContent( (char*)disp_ss.str().c_str(),0);
		contactPlate->AddElement("<SlvNodes>");
		contactPlate->GetChildren()[4]->AddVariable( new XMLVariable("NumNodes", (char*) n_slaves.str().c_str()));
		contactPlate->GetChildren()[4]->AddContent( (char*)slvs_ss.str().c_str(),0);

	this->count += 1;

	// == clear ==
	a_ss.str(""); //a_ss.clear() ;
	b_ss.str(""); //b_ss.clear() ;
	c_ss.str(""); //c_ss.clear() ;
	  //
	disp_ss.str(""); //disp_ss.clear() ;
	n_slaves.str(""); //n_slaves.clear() ;
	slvs_ss.str("");//slvs_ss.clear();
}

void xmlModelWriter::addSystemParams(char* parameter, char* value)
{ 
	XMLElement* systemParams = this->document->GetRootElement()->GetChildren()[ 0 ];
		systemParams->AddElement( parameter );
		systemParams->GetChildren()[ this->children_param ]->AddContent(value, 0);
	
		this->children_param +=1;
}

void xmlModelWriter::addGravity(std::list<int> nodes, char* AccelerationMagnitude, std::vector<double> AccelerationDirection)
{
		// stringstream n_ss; 
	 ss << nodes.size();
	// std::ostringstream nodes_ss;
	 std::copy(nodes.begin(), nodes.end(), std::ostream_iterator<int>(oss, " "));
	//

	// stringstream p_ss;
	// p_ss << AccelerationMagnitude;
	// std::ostringstream p_oss;
	std::copy(AccelerationDirection.begin(), AccelerationDirection.end(), std::ostream_iterator<double>(p_oss, " "));

	this->root->AddElement("<Constraint>");
	XMLElement* Constraint = this->document->GetRootElement()->GetChildren()[ this->count ];
		// Constraint->AddVariable( new XMLVariable("DOF", Dof));
		Constraint->AddVariable( new XMLVariable("LoadShape", "POLY345"));
		Constraint->AddVariable( new XMLVariable("NumNodes", (char*)ss.str().c_str() ));
		Constraint->AddVariable( new XMLVariable("Type", "Gravity"));
		Constraint->AddElement("<Nodes>");
			Constraint->GetChildren()[0]->AddContent( (char*)oss.str().c_str() ,0);
		Constraint->AddElement("<AccelerationMagnitude>");
			Constraint->GetChildren()[1]->AddContent( (char*) AccelerationMagnitude ,0);
		Constraint->AddElement("<AccelerationDirection>");
			Constraint->GetChildren()[2]->AddContent( (char*)p_oss.str().c_str() ,0);
	this->count += 1;

	// == clear ==
	ss.str(""); //ss.clear();
	oss.str(""); //oss.clear();

	p_oss.str("");
}

void xmlModelWriter::addOutput(char* output, char* freq)
{
	this->root->AddElement("<Output>");
	XMLElement* Output =  this->document->GetRootElement()->GetChildren()[ this->count ];
		Output->AddVariable( new XMLVariable( "Freq", freq));
		Output->AddElement("<Variable>");
			Output->GetChildren()[0]->AddContent( output, 0);

	this->count += 1;
}

void xmlModelWriter::writeModel( char* outputfilename)
{
	this->document->Save( outputfilename );
}

#endif