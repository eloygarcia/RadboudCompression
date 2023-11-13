// From breastCompress.cxx
// lines 1388-1994
#include <stdio.h>

#include <iostream>
#include <fstream>
using namespace std;

int main( int argc, char* argv[]){

    std::string febioFilename = argv[1];

    std::cout << "Writing FEBio simulation file..." << std<::endl;
    ofstream febFile;

    //char febioFilename[256];
    sprintf(febioFilename.c_str(), "%s/febio_%d_%d.feb", workDir.c_str(), seed, meshCount);
    febFile.open(febioFilename.c_str());

    // header
    febFile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    febFile << "<febio_spec version=\"2.0\">\n";
    febFile << "\t<Module type=\"solid\"/>\n";

    // control
    int timeSteps;
    double dtMin;
//    if (autoRemesh){
//      timeSteps=20;
//      dtMin = 0.005;
//    } else {
      timeSteps = (int)(ceil(100/numMeshes)); // Que es numMeshes??
      dtMin = 0.001;
//    }
    double stepSize = 1.0/((double)timeSteps);

    febFile << "\t<Control>\n";
    febFile << "\t\t<time_steps>" << timeSteps << "</time_steps>\n";
    febFile << "\t\t<step_size>" << stepSize << "</step_size>\n";
    //no sé que es esto
    febFile << "\t\t<max_refs>10</max_refs>\n";
    febFile << "\t\t<max_ups>10</max_ups>\n";

    //tolerancias ??
    febFile << "\t\t<dtol>0.001</dtol>\n";
    febFile << "\t\t<etol>0.01</etol>\n";
    febFile << "\t\t<rtol>0</rtol>\n";
    febFile << "\t\t<lstol>0.9</lstol>\n";

    febFile << "\t\t<time_stepper>\n";
    febFile << "\t\t\t<dtmin>" << dtMin << "</dtmin>\n";
    febFile << "\t\t\t<dtmax>" << stepSize << "</dtmax>\n";

    // can lower max_retries to either pass or fail faster.
//    if(autoRemesh){
//      febFile << "\t\t\t<max_retries>7</max_retries>\n";
//    } else {
      febFile << "\t\t\t<max_retries>50</max_retries>\n";
//    }

    febFile << "\t\t\t<opt_iter>10</opt_iter>\n";
    febFile << "\t\t</time_stepper>\n";

     // logs??
    febFile << "\t\t<print_level>PRINT_MAJOR_ITRS</print_level>\n";
    febFile << "\t\t<plot_level>PLOT_MAJOR_ITRS</plot_level>\n";
    febFile << "\t\t<analysis type=\"static\"/>\n";
    febFile << "\t</Control>\n";

    // global variables
    febFile << "\t<Globals>\n";
    febFile << "\t\t<Constants>\n";
    febFile << "\t\t\t<T>0</T>\n";
    febFile << "\t\t\t<R>0</R>\n";
    febFile << "\t\t\t<Fc>0</Fc>\n";
    febFile << "\t\t</Constants>\n";
    febFile << "\t</Globals>\n";

    // Materials
    febFile << "\t<Material>\n";
    febFile << "\t\t<material id=\"1\" name=\"BreastFat\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << fatDensity << "</density>\n";
    febFile << "\t\t\t<E>" << fatModulus << "</E>\n";
    febFile << "\t\t\t<v>" << fatPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"2\" name=\"BreastGland\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << glandDensity << "</density>\n";
    febFile << "\t\t\t<E>" << glandModulus << "</E>\n";
    febFile << "\t\t\t<v>" << glandPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"3\" name=\"BottomPaddleMaterial\" type=\"rigid body\">\n";
    febFile << "\t\t\t<density>1e-05</density>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"4\" name=\"TopPaddleMaterial\" type=\"rigid body\">\n";
    febFile << "\t\t\t<density>1e-05</density>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"5\" name=\"BreastMuscle\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << fatDensity << "</density>\n";
    febFile << "\t\t\t<E>" << fatModulus << "</E>\n";
    febFile << "\t\t\t<v>" << fatPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"6\" name=\"BreastMass\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << massDensity << "</density>\n";
    febFile << "\t\t\t<E>" << massModulus << "</E>\n";
    febFile << "\t\t\t<v>" << massPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t</Material>\n";

    // Geometry
    febFile << "\t<Geometry>\n";

    // nodes
    febFile << "\t\t<Nodes>\n";

    int numBottomPaddleNodes = bottomPaddleNodes->GetNumberOfTuples();

    // bottom paddle nodes
    for(int i=0; i<numBottomPaddleNodes; i++){
      double loc[3];
      bottomPaddleNodes->GetTuple(i,loc);
      febFile << "\t\t\t<node id=\"" << i+1 << "\"> " << loc[0] << ", " <<
	loc[1] << ", " << loc[2] << "</node>\n";
    }

    int numTopPaddleNodes = topPaddleNodes->GetNumberOfTuples();

    // top paddle nodes
    for(int i=0; i<numTopPaddleNodes; i++){
      double loc[3];
      topPaddleNodes->GetTuple(i,loc);
      febFile << "\t\t\t<node id=\"" << i+1+numBottomPaddleNodes << "\"> " << loc[0] << ", " <<
	loc[1] << ", " << loc[2] << "</node>\n";
    }

    // breast nodes
    REAL *breastPt = tetout.pointlist;

    for(int i=0; i<tetout.numberofpoints; i++){
      febFile << "\t\t\t<node id=\"" << i+1+numBottomPaddleNodes+numTopPaddleNodes << "\"> " << *breastPt;
      breastPt++;
      febFile << ", " << *breastPt;
      breastPt++;
      febFile << ", " << *breastPt << "</node>\n";
      breastPt++;
    }

    febFile << "\t\t</Nodes>\n";

    // elements
    int elementCount = 0;

    // bottom paddle hex elements
    febFile << "\t\t<Elements type=\"hex8\" mat=\"3\" elset=\"bottomPaddleHex\">\n";
    for(int i=0; i<bottomPaddleHexElements->GetNumberOfTuples(); i++){
      double *els;
      els = bottomPaddleHexElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i+1+elementCount << "\"> ";
      for(int j=0; j<7; j++){
        febFile << (int)(els[j]) << ", ";
      }
      febFile << (int)(els[7]) << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += bottomPaddleHexElements->GetNumberOfTuples();

    // bottom paddle pent elements
    febFile << "\t\t<Elements type=\"penta6\" mat=\"3\" elset=\"bottomPaddlePent\">\n";
    for(int i=0; i<bottomPaddlePentElements->GetNumberOfTuples(); i++){
      double *els;
      els = bottomPaddlePentElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i+1+elementCount << "\"> ";
      for(int j=0; j<5; j++){
        febFile << (int)(els[j]) << ", ";
      }
      febFile << (int)(els[5]) << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += bottomPaddlePentElements->GetNumberOfTuples();

    // top paddle hex elements
    febFile << "\t\t<Elements type=\"hex8\" mat=\"4\" elset=\"topPaddleHex\">\n";
    for(int i=0; i<topPaddleHexElements->GetNumberOfTuples(); i++){
      double *els;
      els = topPaddleHexElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i+1+elementCount << "\"> ";
      for(int j=0; j<7; j++){
        febFile << (int)(els[j])+numBottomPaddleNodes << ", ";
      }
      febFile << (int)(els[7])+numBottomPaddleNodes << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += topPaddleHexElements->GetNumberOfTuples();

    // top paddle pent elements
    febFile << "\t\t<Elements type=\"penta6\" mat=\"4\" elset=\"topPaddlePent\">\n";
    for(int i=0; i<topPaddlePentElements->GetNumberOfTuples(); i++){
      double *els;
      els = topPaddlePentElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i+1+elementCount << "\"> ";
      for(int j=0; j<5; j++){
        febFile << (int)(els[j])+numBottomPaddleNodes << ", ";
      }
      febFile << (int)(els[5])+numBottomPaddleNodes << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += topPaddlePentElements->GetNumberOfTuples();

    // determine tissue type for breast tets and set muscle node list
    vtkSmartPointer<vtkIntArray> tetTypes =
      vtkSmartPointer<vtkIntArray>::New();

    tetTypes->SetNumberOfComponents(1);
    tetTypes->SetNumberOfTuples(tetout.numberoftetrahedra);

    vtkSmartPointer<vtkIdList> muscleNodes =
      vtkSmartPointer<vtkIdList>::New();

    bool fatFound = false;
    bool muscleFound = false;
    bool glandFound = false;
    bool massFound = false;

    for(int i=0; i<tetout.numberoftetrahedra; i++){
      double coords[4][3];
      double center[3];

      for(int j=0; j<4; j++){
        int idx = tetout.tetrahedronlist[4*i+j];
        for(int k=0; k<3; k++){
          coords[j][k] = tetout.pointlist[idx*3+k];
        }
      }

      // calculate center of tet
      vtkTetra::TetraCenter(coords[0], coords[1], coords[2], coords[3], center);

      // find tissue type
      //vtkIdType loc = input->FindPoint(center);
      double pcoords[3];
      int ijk[3];
      int inFOV = input->ComputeStructuredCoordinates(center, ijk, pcoords);

      if(inFOV){
        unsigned char* p = static_cast<unsigned char*>(input->GetScalarPointer(ijk));

        if (p[0] == tissue.fat || p[0] == tissue.skin){
          // use fat material
          tetTypes->SetTuple1(i, 1);
          fatFound=true;
        } else if (p[0] == tissue.muscle) {
          // use fat and add nodes to fixed boundary
	  tetTypes->SetTuple1(i, 5);
          muscleFound=true;
          for(int j=0; j<4; j++){
            muscleNodes->InsertUniqueId(tetout.tetrahedronlist[4*i+j]);
          }
        } else if(p[0] == tissue.mass) {
          tetTypes->SetTuple1(i, 6);
          massFound=true;
        } else if (p[0] == tissue.gland || p[0] == tissue.TDLU ||
		   p[0] == tissue.duct || p[0] == tissue.nipple || p[0] == tissue.cooper
		   || p[0] == tissue.artery || p[0] == tissue.vein){
          // use gland
          tetTypes->SetTuple1(i, 2);
          glandFound=true;
        } else if (p[0] == tissue.bg || p[0] == tissue.paddle){
	  tetTypes->SetTuple1(i, 1);
        } else {
          cerr << "Breast tetrahedron center material = " << p[0] << " unrecognized\n";
          return(1);
        }
      } else {
        // center outside FOV, make it fat
        tetTypes->SetTuple1(i, 1);
      }
    }

    // breast fat elements
    if(fatFound){
      febFile << "\t\t<Elements type=\"tet4\" mat=\"1\" elset=\"breastFatTet\">\n";
      for(int i=0; i<tetout.numberoftetrahedra; i++){
        if((int)tetTypes->GetTuple1(i) == 1){
          febFile << "\t\t\t<elem id=\"" << elementCount+1 << "\"> ";
          int n;
          for(int j=0; j<3; j++){
            febFile << tetout.tetrahedronlist[4*i+j]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4*i+3]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // breast gland elements
    if(glandFound){
      febFile << "\t\t<Elements type=\"tet4\" mat=\"2\" elset=\"breastGlandTet\">\n";
      for(int i=0; i<tetout.numberoftetrahedra; i++){
        if((int)tetTypes->GetTuple1(i) == 2){
          febFile << "\t\t\t<elem id=\"" << elementCount+1 << "\"> ";
          int n;
          for(int j=0; j<3; j++){
            febFile << tetout.tetrahedronlist[4*i+j]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4*i+3]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // mass elements
    if(massFound){
      febFile << "\t\t<Elements type=\"tet4\" mat=\"6\" elset=\"breastMassTet\">\n";
      for(int i=0; i<tetout.numberoftetrahedra; i++){
        if((int)tetTypes->GetTuple1(i) == 6){
          febFile << "\t\t\t<elem id=\"" << elementCount+1 << "\"> ";
          int n;
          for(int j=0; j<3; j++){
            febFile << tetout.tetrahedronlist[4*i+j]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4*i+3]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // muscle elements
    if(muscleFound){
      febFile << "\t\t<Elements type=\"tet4\" mat=\"5\" elset=\"breastMuscleTet\">\n";
      for(int i=0; i<tetout.numberoftetrahedra; i++){
        if((int)tetTypes->GetTuple1(i) == 5){
          febFile << "\t\t\t<elem id=\"" << elementCount+1 << "\"> ";
          int n;
          for(int j=0; j<3; j++){
            febFile << tetout.tetrahedronlist[4*i+j]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4*i+3]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // end geometry
    febFile << "\t</Geometry>\n";

    // boundary - muscle tissue fixed
    febFile << "\t<Boundary>\n";

    febFile << "\t\t<fix bc=\"xyz\">\n";
    for(int i=0; i<muscleNodes->GetNumberOfIds(); i++){
      vtkIdType n;
      n = muscleNodes->GetId(i);
      febFile << "\t\t\t<node id=\"" << (n+numBottomPaddleNodes+numTopPaddleNodes+1) << "\"/>\n";
    }
    febFile << "\t\t</fix>\n";

    febFile << "\t\t<fix bc=\"uvw\">\n";
    for(int i=0; i<muscleNodes->GetNumberOfIds(); i++){
      vtkIdType n;
      n = muscleNodes->GetId(i);
      febFile << "\t\t\t<node id=\"" << (n+numBottomPaddleNodes+numTopPaddleNodes+1) << "\"/>\n";
    }
    febFile << "\t\t</fix>\n";

    febFile << "\t</Boundary>\n";

    // Gravity
    febFile << "\t<Loads>\n";
    febFile << "\t\t<body_load type=\"const\">\n";
    febFile << "\t\t\t<z lc=\"1\">4905</z>\n"; //gravity in kg*mm/s^2 assuming 0.5 kg breast
    febFile << "\t\t</body_load>\n";
    febFile << "\t</Loads>\n";

    // Contacts

    // bottom contact
    febFile << "\t<Contact>\n";
    febFile << "\t\t<contact type=\"facet-to-facet sliding\">\n";
    febFile << "\t\t\t<laugon>0</laugon>\n";

    febFile << "\t\t\t<tolerance>1.0</tolerance>\n";
    febFile << "\t\t\t<penalty>25</penalty>\n";
    febFile << "\t\t\t<two_pass>0</two_pass>\n";
    febFile << "\t\t\t<auto_penalty>1</auto_penalty>\n";
    febFile << "\t\t\t<fric_coeff>0</fric_coeff>\n";
    febFile << "\t\t\t<fric_penalty>0</fric_penalty>\n";
    febFile << "\t\t\t<search_tol>0.02</search_tol>\n";
    febFile << "\t\t\t<minaug>0</minaug>\n";
    febFile << "\t\t\t<maxaug>10</maxaug>\n";
    febFile << "\t\t\t<gaptol>0</gaptol>\n";
    febFile << "\t\t\t<seg_up>0</seg_up>\n";
    febFile << "\t\t\t<surface type=\"master\">\n";
    for(int i=0; i<bottomPaddleContactSurface->GetNumberOfTuples(); i++){
      double *n;
      n = bottomPaddleContactSurface->GetTuple4(i);
      febFile << "\t\t\t\t<quad4 id=\"" << i+1 << "\"> " <<
	(int)(n[0]) << ", " <<
	(int)(n[1]) << ", " <<
	(int)(n[2]) << ", " <<
	(int)(n[3]) << "</quad4>\n";
    }
    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t\t<surface type=\"slave\">\n";
    int triCount = 1;

    double zThresh = nipplePos[2];

    for(int i=0; i<tetout.numberoftrifaces; i++){
      // test if triface is bottom contact
      double meanZ = 0.0;
      double meanX = 0.0;
      double meanY = 0.0;
      for(int j=0; j<3; j++){
        meanZ += tetout.pointlist[3*(tetout.trifacelist[3*i+j])+2];
        meanY += tetout.pointlist[3*(tetout.trifacelist[3*i+j])+1];
        meanX += tetout.pointlist[3*(tetout.trifacelist[3*i+j])];
      }
      meanZ = meanZ/3.0;
      meanY = meanY/3.0;
      meanX = meanX/3.0;

      if(rotate){
        meanZ = -sin(angle)*meanY + cos(angle)*meanZ;
        zThresh = -sin(angle)*nipplePos[1] + cos(angle)*nipplePos[2];
      }

      if(meanZ < zThresh){
        febFile << "\t\t\t\t<tri3 id=\"" << triCount << "\"> " <<
	  tetout.trifacelist[3*i]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", " <<
	  tetout.trifacelist[3*i+1]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", " <<
	  tetout.trifacelist[3*i+2]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</tri3>\n";
        triCount += 1;
      }
    }

    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t</contact>\n";

    // top contact
    febFile << "\t\t<contact type=\"facet-to-facet sliding\">\n";
    febFile << "\t\t\t<laugon>0</laugon>\n";
    febFile << "\t\t\t<tolerance>1.0</tolerance>\n";
    febFile << "\t\t\t<penalty>25</penalty>\n";
    febFile << "\t\t\t<two_pass>0</two_pass>\n";
    febFile << "\t\t\t<auto_penalty>1</auto_penalty>\n";
    febFile << "\t\t\t<fric_coeff>0</fric_coeff>\n";
    febFile << "\t\t\t<fric_penalty>0</fric_penalty>\n";
    febFile << "\t\t\t<search_tol>0.02</search_tol>\n";
    febFile << "\t\t\t<minaug>0</minaug>\n";
    febFile << "\t\t\t<maxaug>10</maxaug>\n";
    febFile << "\t\t\t<gaptol>0</gaptol>\n";
    febFile << "\t\t\t<seg_up>0</seg_up>\n";
    febFile << "\t\t\t<surface type=\"master\">\n";
    for(int i=0; i<topPaddleContactSurface->GetNumberOfTuples(); i++){
      double *n;
      n = topPaddleContactSurface->GetTuple4(i);
      febFile << "\t\t\t\t<quad4 id=\"" << i+1 << "\"> " <<
              (int)(n[0])+numBottomPaddleNodes << ", " <<
              (int)(n[1])+numBottomPaddleNodes << ", " <<
              (int)(n[2])+numBottomPaddleNodes << ", " <<
              (int)(n[3])+numBottomPaddleNodes << "</quad4>\n";
    }
    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t\t<surface type=\"slave\">\n";
    triCount = 1;

    for(int i=0; i<tetout.numberoftrifaces; i++){
      // test if triface is top contact
      double meanZ = 0.0;
      double meanX = 0.0;
      double meanY = 0.0;
      for(int j=0; j<3; j++){
        meanZ += tetout.pointlist[3*(tetout.trifacelist[3*i+j])+2];
        meanY += tetout.pointlist[3*(tetout.trifacelist[3*i+j])+1];
        meanX += tetout.pointlist[3*(tetout.trifacelist[3*i+j])];
      }
      meanZ = meanZ/3.0;
      meanY = meanY/3.0;
      meanX = meanX/3.0;

      if(rotate){
        meanZ = -sin(angle)*meanY + cos(angle)*meanZ;
        zThresh = -sin(angle)*nipplePos[1] + cos(angle)*nipplePos[2];
      }

      if(meanZ >= zThresh){
        febFile << "\t\t\t\t<tri3 id=\"" << triCount << "\"> " <<
                tetout.trifacelist[3*i]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", " <<
                tetout.trifacelist[3*i+1]+numBottomPaddleNodes+numTopPaddleNodes+1 << ", " <<
                tetout.trifacelist[3*i+2]+numBottomPaddleNodes+numTopPaddleNodes+1 << "</tri3>\n";
        triCount += 1;
      }
    }

    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t</contact>\n";

    febFile << "\t</Contact>\n";

    // constraints
    febFile << "\t<Constraints>\n";
    febFile << "\t\t<rigid_body mat=\"4\">\n";
    febFile << "\t\t\t<fixed bc=\"x\"/>\n";
    if(!rotate){
      febFile << "\t\t\t<fixed bc=\"y\"/>\n";
    }
    febFile << "\t\t\t<fixed bc=\"Rx\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Ry\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Rz\"/>\n";
    febFile << "\t\t</rigid_body>\n";
    febFile << "\t\t<rigid_body mat=\"3\">\n";
    febFile << "\t\t\t<fixed bc=\"x\"/>\n";
    if(!rotate){
      febFile << "\t\t\t<fixed bc=\"y\"/>\n";
    }
    febFile << "\t\t\t<fixed bc=\"Rx\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Ry\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Rz\"/>\n";
    febFile << "\t\t</rigid_body>\n";
    febFile << "\t\t<rigid_body mat=\"4\">\n";
    febFile << "\t\t\t<prescribed bc=\"z\" lc=\"2\">1</prescribed>\n";
    febFile << "\t\t</rigid_body>\n";
    if(rotate){
      febFile << "\t\t<rigid_body mat=\"4\">\n";
      febFile << "\t\t\t<prescribed bc=\"y\" lc=\"4\">1</prescribed>\n";
      febFile << "\t\t</rigid_body>\n";
    }
    febFile << "\t\t<rigid_body mat=\"3\">\n";
    febFile << "\t\t\t<prescribed bc=\"z\" lc=\"3\">1</prescribed>\n";
    febFile << "\t\t</rigid_body>\n";
    if(rotate){
      febFile << "\t\t<rigid_body mat=\"3\">\n";
      febFile << "\t\t\t<prescribed bc=\"y\" lc=\"5\">1</prescribed>\n";
      febFile << "\t\t</rigid_body>\n";
    }
    febFile << "\t</Constraints>\n";

    // load curves
    febFile << "\t<LoadData>\n";
    febFile << "\t\t<loadcurve id=\"1\" type=\"linear\">\n";
    febFile << "\t\t\t<point>0,0</point>\n";

    // currently have kept gravity very low but non-zero for stability.
    //    if(meshCount == 0)
    //    {
    //      febFile << "\t\t\t<point>0.3,1.0</point>\n";
    //      febFile << "\t\t\t<point>1,1.0</point>\n";
    //    }
    //    else
    //    {

    //      febFile << "\t\t\t<point>0.3,0</point>\n";
    //      febFile << "\t\t\t<point>0.4,1.0</point>\n";

    // Left a miniscule force so that the breast wouldn't be underconstrained
    // (causing problem to diverge).
    febFile << "\t\t\t<point>1,0.001</point>\n";
    //    }

    febFile << "\t\t</loadcurve>\n";
    if(rotate){
      febFile << "\t\t<loadcurve id=\"2\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << -paddleDist*cos(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"3\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << paddleDist*cos(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"4\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << paddleDist*sin(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"5\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << -paddleDist*sin(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
    } else {
      // meet in center, each paddle travels an equal distance
      febFile << "\t\t<loadcurve id=\"2\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1, " << -paddleDist << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"3\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1, " << paddleDist << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
    }
    febFile << "\t</LoadData>\n";

    // output
    febFile << "\t<Output>\n";
    febFile << "\t\t<plotfile type=\"febio\">\n";
    febFile << "\t\t\t<var type=\"displacement\"/>\n";
    febFile << "\t\t\t<compression>0</compression>\n";
    febFile << "\t\t</plotfile>\n";
    febFile << "\t</Output>\n";

    // end feb file
    febFile << "</febio_spec>\n";
    febFile.close();

    /************
     * run FEBio
     ***********/

}