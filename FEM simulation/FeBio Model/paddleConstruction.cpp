// from breastCompress.cxx
// lines 610 - 1388
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
//#include <vtkIdType.h>
#include <vtkIntArray.h>
#include <vtkMath.h>

#include <vtkUnstructuredGrid.h>
//#include "vtkUnstructuredGridWriter.h"
//#include <vtkUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
//#include <vtk
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>

#include <vtkCamera.h>
#include <vtkCellType.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>



int main(int argc, char* argv[]){
    // Eloy
    double phantomBounds[6] = {0.0, 100.0, 0.0, 100.0, 0.0, 100.0};
    double origPaddleWidth = 150.0; // Y
    double origPaddleLength = 150.0; // X
    double nipplePos[3] = {50.0, 50.0, 50.0};
    int meshCount = 0;
    int numMeshes  = 0;
    bool rotate = false;
    bool autoRemesh = false;
    double prevTopPaddleCenterDerotated[3] = {0.0,0.0,0.0};
    double prevBotPaddleCenterDerotated[3] = {0.0,0.0,0.0};
    double paddleDist = 0.0;
    double thickness = 50.;
    double angle = 0.0;
    double pi = vtkMath::Pi();
    int paddleNh = 18;
    int paddleNl = 25;
    int paddleNw = 25;

    // fin Eloy

    double topPaddleCenter[3]={0.0,0.0,0.0};
    double bottomPaddleCenter[3]={0.0,0.0,0.0};
    double paddleHeight = 1.0;

    if(meshCount == 0){
      origPaddleWidth = (phantomBounds[3]-phantomBounds[2])*1.75;
      origPaddleLength = (phantomBounds[1]-phantomBounds[0])*1.75;
    }

    double paddleWidth = origPaddleWidth;
    double paddleLength = origPaddleLength;
    double paddleEdgeRadius = 3.0;
    double paddleLipHeight = 15.0; // must be greater than paddleEdgeRadius

    // paddle retreat after remeshing
    double retreatDist = 10.0*paddleHeight;

    if(meshCount == 0){
      topPaddleCenter[0] = paddleLength/2.0 + (phantomBounds[1]-phantomBounds[0])*0.125;
      topPaddleCenter[1] = nipplePos[1];
      if(!rotate){
        topPaddleCenter[2] = phantomBounds[5] + 4*paddleHeight/2.0;
      } else {
        topPaddleCenter[2] = phantomBounds[5] + 30*paddleHeight/2.0;
      }

      bottomPaddleCenter[0] = paddleLength/2.0  + (phantomBounds[1]-phantomBounds[0])*0.1;
      bottomPaddleCenter[1] = nipplePos[1];
      if(!rotate){
        bottomPaddleCenter[2] = phantomBounds[4] - 4*paddleHeight/2.0;
      } else {
        bottomPaddleCenter[2] = phantomBounds[4] - 30*paddleHeight/2.0;
      }
    } else {
      // use previous interation's paddle position plus a retreat
      for(int i=0; i<3; i++){
        topPaddleCenter[i] = prevTopPaddleCenterDerotated[i];
        bottomPaddleCenter[i] = prevBotPaddleCenterDerotated[i];
      }
      topPaddleCenter[2] += retreatDist;
      bottomPaddleCenter[2] -= retreatDist;
    }

    if(autoRemesh){
      // try to complete full compression
      paddleDist = ((topPaddleCenter[2]-bottomPaddleCenter[2] -
                     paddleHeight - thickness)/2.0);
    } else {
      if(numMeshes >= 3)
      {
        // paddle distance must be greater for every iteration past the first one to
        // account for the distancing done when resuming from the last iteration.
        if(meshCount == 0){
          paddleDist = ((topPaddleCenter[2]-bottomPaddleCenter[2] -
                         paddleHeight - thickness)/2.0)*0.4;
        } else if(meshCount == 1) {
          paddleDist = paddleDist/2;
          paddleDist += retreatDist;
        } else {
          paddleDist -= retreatDist;
          paddleDist = (paddleDist*2)/(numMeshes-2);
          paddleDist += retreatDist;
        }
      } else {
	// 1 or 2 meshes
        if(meshCount==0) {
          paddleDist = ((topPaddleCenter[2]-bottomPaddleCenter[2] -
                         paddleHeight - thickness)/2.0)/numMeshes;
        } else if (meshCount==1) {
          paddleDist+=retreatDist;
        }
      }
    }

/*
    if(paddleDist < 0.0){
      cerr << "Paddle travel distance error.\n";
      return(1);
    }
*/

    //vtkSmartPointer<vtkDoubleArray> bottomPaddleNodes = vtkSmartPointer<vtkDoubleArray>::New();
    vtkNew<vtkPoints> bottomPaddleNodes;
        //bottomPaddleNodes->SetNumberOfComponents(3);

    //vtkSmartPointer<vtkIntArray> bottomPaddleHexElements = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkCellArray> bottomPaddleHexElements = vtkSmartPointer<vtkCellArray>::New();
        //bottomPaddleHexElements->SetNumberOfComponents(8);

    // vtkSmartPointer<vtkIntArray> bottomPaddlePentElements = vtkSmartPointer<vtkIntArray>::New();
//    vtkSmartPointer<vtkCellArray> bottomPaddlePentElements = vtkSmartPointer<vtkCellArray>::New(); // skip
        //bottomPaddlePentElements->SetNumberOfComponents(6);

    // vtkSmartPointer<vtkIntArray> bottomPaddleContactSurface = vtkSmartPointer<vtkIntArray>::New();
//     vtkSmartPointer<vtkCellArray> bottomPaddleContactSurface = vtkSmartPointer<vtkCellArray>::New(); // skip
        //bottomPaddleContactSurface->SetNumberOfComponents(4);

//    vtkSmartPointer<vtkIntArray> breastBackNodes = vtkSmartPointer<vtkIntArray>::New(); // skip
//        breastBackNodes->SetNumberOfComponents(1);

    double paddleDy = paddleWidth/paddleNw;
    double paddleDx = (paddleLength-2.0*paddleEdgeRadius)/(paddleNl-1);
    vtkIdType nodes[8];
    vtkIdType pnodes[6];

    // create main part of bottom paddle
    for(int j=0; j<paddleNl-1; j++){
      double xpos = bottomPaddleCenter[0]+paddleLength/2.0 - j*paddleDx;
      for(int i=0; i<paddleNw; i++){
        double ypos = bottomPaddleCenter[1]-paddleWidth/2.0 + i*paddleDy;
        double loc[3]={0.0,0.0,0.0};
        // add nodes
        if(i==0 && j==0){ //j e i son cero
          // add all 8 nodes
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0; // inferior (menos)
          std::cout << "point {" << loc[0] << ","<< loc[1] << ","<< loc[2] << "}" << std::endl;
          //nodes[0] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0; //(menos)
          //nodes[1] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0; // menos
          //nodes[2] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0; //menos
          //nodes[3] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0; // más
          //nodes[4] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0; // mas
          //nodes[5] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0; //mas
          //nodes[6] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0; //más
          //nodes[7] =
          bottomPaddleNodes->InsertNextPoint(loc);
        } else if(j==0){
          // add 4 nodes
          //nodes[0] = nodes[1];
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0;
          //nodes[1] =
          bottomPaddleNodes->InsertNextPoint(loc);

          //nodes[3] = nodes[2];
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0;
          //nodes[2] =
          bottomPaddleNodes->InsertNextPoint(loc);

          //nodes[4] = nodes[5];
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
          //nodes[5] =
          bottomPaddleNodes->InsertNextPoint(loc);

          //nodes[7] = nodes[6];
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
          //nodes[6] =
          bottomPaddleNodes->InsertNextPoint(loc);
        } else if(i==0){
          // add 4 nodes
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0;
          //nodes[2] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0;
          //nodes[3] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
          //nodes[6] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
          //nodes[7] =
          bottomPaddleNodes->InsertNextPoint(loc);

          // get left neighbor nodes
          // double leftN[8];
          // bottomPaddleHexElements->GetCell((j-1)*(paddleNw), leftN);
          //nodes[0] = (int)leftN[3]-1;
          //nodes[1] = (int)leftN[2]-1;
          //nodes[4] = (int)leftN[7]-1;
          //nodes[5] = (int)leftN[6]-1;
        } else {
          // add 2 nodes
          //nodes[0] = nodes[1];
          //nodes[3] = nodes[2];
          //nodes[4] = nodes[5];
          //nodes[7] = nodes[6];

          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]-paddleHeight/2.0;
          //nodes[2] =
          bottomPaddleNodes->InsertNextPoint(loc);

          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
          //nodes[6] =
          bottomPaddleNodes->InsertNextPoint(loc);

          // get left neighbor nodes
          //double *leftN;
          //leftN = bottomPaddleHexElements->GetCell((j-1)*(paddleNw)+i);
          //nodes[1] = (int)leftN[2]-1;
          //nodes[5] = (int)leftN[6]-1;
        }
*/
        // add element
/*        double nodesD[8];
        for(int a=0; a<8; a++){
          nodesD[a] = (double)(nodes[a]+1);
        }
        //bottomPaddleHexElements->InsertNextCell(nodesD);
        bottomPaddleHexElements->InsertNextCell(nodesD);

        // add surface
 //       bottomPaddleContactSurface->InsertNextTuple4(nodesD[4],nodesD[5],nodesD[6],nodesD[7]);
      }
    }
*/
}
}
    int saveTop;

/*
    // curved section
    for(int i=0; i<paddleNw; i++){
      double ypos = bottomPaddleCenter[1]-paddleWidth/2.0 + i*paddleDy;
      double xpos = bottomPaddleCenter[0]+paddleLength/2.0 - (paddleNl-1)*paddleDx;

      // middle pentahedral
      double leftN[8];
      double loc[3];
      bottomPaddleHexElements->GetCell((paddleNl-2)*(paddleNw)+i, leftN);

      if(i==0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius;
        pnodes[5] = bottomPaddleNodes->InsertNextPoint(loc);
      } else {
        pnodes[5] = pnodes[0];
      }
      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius;
      pnodes[2] = bottomPaddleNodes->InsertNextPoint(loc);
      pnodes[3] = (int)leftN[3]-1;
      pnodes[4] = (int)leftN[7]-1;
      pnodes[0] = (int)leftN[2]-1;
      pnodes[1] = (int)leftN[6]-1;

      // save 2 nodes for top pent
      int nodSave[2];
      nodSave[0] = pnodes[1];
      nodSave[1] = pnodes[4];

      double pnodesD[6];
      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }

      bottomPaddlePentElements->InsertNextCell(pnodesD);

      // bottom pentahedral
      // 0 and 3 the same
      pnodes[1] = pnodes[2];
      pnodes[4] = pnodes[5];
      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleLipHeight;
      pnodes[2] = bottomPaddleNodes->InsertNextPoint(loc);
      if(i == 0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleLipHeight;
        pnodes[5] = bottomPaddleNodes->InsertNextPoint(loc);
      } else {
        // from previous edge hex
        pnodes[5] = nodes[1];
      }

      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }
      bottomPaddlePentElements->InsertNextCell(pnodesD);

      // bottom edge hex
      nodes[0] = pnodes[5];
      nodes[4] = pnodes[4];
      nodes[1] = pnodes[2];
      nodes[5] = pnodes[1];
      if(i==0) {
        loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleLipHeight;
        nodes[3] = bottomPaddleNodes->InsertNextPoint(loc);

        loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius;
        nodes[7] = bottomPaddleNodes->InsertNextPoint(loc);
      } else {
        nodes[3] = nodes[2];
        nodes[7] = nodes[6];
      }

      loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleLipHeight;
      nodes[2] = bottomPaddleNodes->InsertNextPoint(loc);

      loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius;
      nodes[6] = bottomPaddleNodes->InsertNextPoint(loc);

      double nodesD[8];
      for(int a=0; a<8; a++){
        nodesD[a] = (double)(nodes[a]+1);
      }
      bottomPaddleHexElements->InsertNextCell(nodesD);

      // the edge
      bottomPaddleContactSurface->InsertNextCell4(nodesD[3],nodesD[7],nodesD[6],nodesD[2]);
      // add the bottom edge as well
      bottomPaddleContactSurface->InsertNextCell4(nodesD[3],nodesD[2],nodesD[1],nodesD[0]);

      // top pentahedral
      // pnodes 1,4 same
      pnodes[2] = nodSave[0];
      pnodes[5] = nodSave[1];

      if(i==0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
        pnodes[3] =  bottomPaddleNodes->InsertNextPoint(loc);
      } else {
        pnodes[3] = saveTop;
      }

      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0;
      pnodes[0] =  bottomPaddleNodes->InsertNextPoint(loc);
      saveTop = pnodes[0];

      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }
      bottomPaddlePentElements->InsertNextCell(pnodesD);

      // the edge
      bottomPaddleContactSurface->InsertNextCell4(pnodesD[5],pnodesD[2],pnodesD[0],pnodesD[3]);

      // the fan
      for(int j=0; j<paddleNh; j++){
        if(j==0){
          // from top pent
          pnodes[1] = pnodes[0];
          pnodes[4] = pnodes[3];
          pnodes[0] = nodes[5];
          pnodes[3] = nodes[4];
        } else {
          pnodes[1] = pnodes[2];
          pnodes[4] = pnodes[5];
        }

        double theta = 3.141592654/2.0/(paddleNh)*(j+1);
        if(j == paddleNh-1){
          pnodes[2] = nodes[6];
        } else {
          loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius*sin(theta);
	  loc[1] = ypos+paddleDy;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius + paddleEdgeRadius*cos(theta);
          pnodes[2] = bottomPaddleNodes->InsertNextPoint(loc);
        }
        if(i==0 && j<paddleNh-1){
          loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius*sin(theta);
	  loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2]+paddleHeight/2.0-paddleEdgeRadius + paddleEdgeRadius*cos(theta);
          pnodes[5] = bottomPaddleNodes->InsertNextPoint(loc);
        } else if(j<paddleNh-1){
          double *leftP;

          leftP = bottomPaddlePentElements->GetCell((i-1)*(3+paddleNh)+3+j);
          pnodes[5] = (int)leftP[2]-1;
        } else {
          pnodes[5] = nodes[7];
        }

        for(int a=0; a<6; a++){
          pnodesD[a] = (double)(pnodes[a]+1);
        }
        bottomPaddlePentElements->InsertNextCell(pnodesD);

        // the edge
        bottomPaddleContactSurface->InsertNextCell4(pnodesD[5],pnodesD[4],pnodesD[1],pnodesD[2]);
      }
    }

    // create top paddle nodes
    // vtkSmartPointer<vtkDoubleArray> topPaddleNodes = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPointsArray> topPaddleNodes = vtkSmartPointer<vtkPointsArray>::New();
        topPaddleNodes->SetNumberOfComponents(3);

    //vtkSmartPointer<vtkIntArray> topPaddleHexElements = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkCellArray> topPaddleHexElements = vtkSmartPointer<vtkCellArray>::New();
        //topPaddleHexElements->SetNumberOfComponents(8);

    //vtkSmartPointer<vtkIntArray> topPaddlePentElements = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkCellArray> topPaddlePentElements = vtkSmartPointer<vtkCellArray>::New();
        topPaddlePentElements->SetNumberOfComponents(6);

    vtkSmartPointer<vtkIntArray> topPaddleContactSurface = vtkSmartPointer<vtkIntArray>::New();
        topPaddleContactSurface->SetNumberOfComponents(4);

    // create main part of top paddle
    for(int j=0; j<paddleNl-1; j++){
      double xpos = topPaddleCenter[0]+paddleLength/2.0 - j*paddleDx;
      for(int i=0; i<paddleNw; i++){
        double ypos = topPaddleCenter[1]-paddleWidth/2.0 + i*paddleDy;
        double loc[3];
        // add nodes
        if(i==0 && j==0){
          // add all 8 nodes
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[0] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[1] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[2] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[3] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[4] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[5] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[6] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[7] = topPaddleNodes->InsertNextPoint(loc);
        } else if(j==0){
          // add 4 nodes
          nodes[0] = nodes[1];
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[1] = topPaddleNodes->InsertNextPoint(loc);
          nodes[3] = nodes[2];
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[2] = topPaddleNodes->InsertNextPoint(loc);
          nodes[4] = nodes[5];
          loc[0] = xpos;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[5] = topPaddleNodes->InsertNextPoint(loc);
          nodes[7] = nodes[6];
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[6] = topPaddleNodes->InsertNextPoint(loc);
        } else if(i==0){
          // add 4 nodes
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[2] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[3] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[6] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[7] = topPaddleNodes->InsertNextPoint(loc);
          // get left neighbor nodes
          double leftN[8];
          topPaddleHexElements->GetCell((j-1)*(paddleNw), leftN);
          nodes[0] = (int)leftN[3]-1;
          nodes[1] = (int)leftN[2]-1;
          nodes[4] = (int)leftN[7]-1;
          nodes[5] = (int)leftN[6]-1;
        } else {
          // add 2 nodes
          nodes[0] = nodes[1];
          nodes[3] = nodes[2];
          nodes[4] = nodes[5];
          nodes[7] = nodes[6];
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
          nodes[2] = topPaddleNodes->InsertNextPoint(loc);
          loc[0] = xpos-paddleDx;
          loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]+paddleHeight/2.0;
          nodes[6] = topPaddleNodes->InsertNextPoint(loc);
          // get left neighbor nodes
          double *leftN;
          leftN = topPaddleHexElements->GetCell((j-1)*(paddleNw)+i);
          nodes[1] = (int)leftN[2]-1;
          nodes[5] = (int)leftN[6]-1;
        }

        // add element
        double nodesD[8];
        for(int a=0; a<8; a++){
          nodesD[a] = (double)(nodes[a]+1);
        }
        topPaddleHexElements->InsertNextCell(nodesD);

        // add surface
        topPaddleContactSurface->InsertNextTuple4(nodesD[0],nodesD[3],nodesD[2],nodesD[1]);
      }
    }

    int saveBottom;

    // curved section
    for(int i=0; i<paddleNw; i++){
      double ypos = topPaddleCenter[1]-paddleWidth/2.0 + i*paddleDy;
      double xpos = topPaddleCenter[0]+paddleLength/2.0 - (paddleNl-1)*paddleDx;

      // middle pentahedral
      double leftN[8];
      double loc[3];
      topPaddleHexElements->GetCell((paddleNl-2)*(paddleNw)+i, leftN);

      if(i==0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius;
        pnodes[5] = topPaddleNodes->InsertNextPoint(loc);
      } else {
        pnodes[5] = pnodes[0];
      }

      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius;
      pnodes[2] = topPaddleNodes->InsertNextPoint(loc);
      pnodes[0] = (int)leftN[2]-1;
      pnodes[1] = (int)leftN[6]-1;
      pnodes[3] = (int)leftN[3]-1;
      pnodes[4] = (int)leftN[7]-1;

      // save 2 nodes for bottom pent
      int nodSave[2];
      nodSave[0] = pnodes[0];
      nodSave[1] = pnodes[3];

      double pnodesD[6];
      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }

      topPaddlePentElements->InsertNextCell(pnodesD);

      // top pentahedral
      int tempNode;
      tempNode = pnodes[2];
      pnodes[2] = pnodes[1];
      pnodes[1] = tempNode;
      tempNode = pnodes[5];
      pnodes[5] = pnodes[4];
      pnodes[4] = tempNode;

      if(i==0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleLipHeight;
	pnodes[3] =  topPaddleNodes->InsertNextPoint(loc);
      } else {
        pnodes[3] = saveTop;
      }

      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleLipHeight;
      pnodes[0] =  topPaddleNodes->InsertNextPoint(loc);
      saveTop = pnodes[0];

      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }
      topPaddlePentElements->InsertNextCell(pnodesD);

      // top edge hex
      nodes[0] = pnodes[4];
      nodes[1] = pnodes[1];
      nodes[4] = pnodes[3];
      nodes[5] = pnodes[0];

      if(i==0){
        loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius;
        nodes[3] = topPaddleNodes->InsertNextPoint(loc);
        loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleLipHeight;
	nodes[7] = topPaddleNodes->InsertNextPoint(loc);
      } else {
        nodes[3] = nodes[2];
        nodes[7] = nodes[6];
      }

      loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius;
      nodes[2] = topPaddleNodes->InsertNextPoint(loc);
      loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleLipHeight;
      nodes[6] = topPaddleNodes->InsertNextPoint(loc);

      double nodesD[8];
      for(int a=0; a<8; a++){
        nodesD[a] = (double)(nodes[a]+1);
      }
      topPaddleHexElements->InsertNextCell(nodesD);

      // the edge
      topPaddleContactSurface->InsertNextCell4(nodesD[3],nodesD[7],nodesD[6],nodesD[2]);
      // add the bottom edge as well
      topPaddleContactSurface->InsertNextCell4(nodesD[7],nodesD[4],nodesD[5],nodesD[6]);

      // bottom pentahedral
      // 1 and 4 the same as top pent
      pnodes[0] = nodSave[0];
      pnodes[3] = nodSave[1];

      if(i == 0){
        loc[0] = xpos-paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
        pnodes[5] = topPaddleNodes->InsertNextPoint(loc);
      } else {
        // from previous edge hex
        pnodes[5] = saveBottom;
      }

      loc[0] = xpos-paddleEdgeRadius;
      loc[1] = ypos+paddleDy;
      loc[2] = topPaddleCenter[2]-paddleHeight/2.0;
      pnodes[2] = topPaddleNodes->InsertNextPoint(loc);
      saveBottom = pnodes[2];

      for(int a=0; a<6; a++){
        pnodesD[a] = (double)(pnodes[a]+1);
      }
      topPaddlePentElements->InsertNextCell(pnodesD);

      // the edge
      topPaddleContactSurface->InsertNextCell4(pnodesD[5],pnodesD[2],pnodesD[0],pnodesD[3]);

      // the fan
      for(int j=0; j<paddleNh; j++){
        if(j==0){
          // from bottom pent
          pnodes[0] = pnodes[1];
          pnodes[3] = pnodes[4];
          // pnodes 2,5 the same
        } else {
          // pnodes 0,3 stay same
          pnodes[2] = pnodes[1];
          pnodes[5] = pnodes[4];
        }

        // need to set 1 and 4
        double theta = pi/2.0/(paddleNh)*(j+1);
        if(j == paddleNh-1){
          pnodes[1] = nodes[2];
        } else {
          loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius*sin(theta);
	  loc[1] = ypos+paddleDy;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius - paddleEdgeRadius*cos(theta);
          pnodes[1] = topPaddleNodes->InsertNextPoint(loc);
        }
        if(i==0 && j<paddleNh-1){
          loc[0] = xpos-paddleEdgeRadius-paddleEdgeRadius*sin(theta);
	  loc[1] = ypos;
          loc[2] = topPaddleCenter[2]-paddleHeight/2.0+paddleEdgeRadius - paddleEdgeRadius*cos(theta);
          pnodes[4] = topPaddleNodes->InsertNextPoint(loc);
        } else if(j<paddleNh-1){
          double *leftP;

          leftP = topPaddlePentElements->GetCell((i-1)*(3+paddleNh)+3+j);
          pnodes[4] = (int)leftP[1]-1;
        } else {

          pnodes[4] = nodes[3];
        }

        for(int a=0; a<6; a++){
          pnodesD[a] = (double)(pnodes[a]+1);
        }
        topPaddlePentElements->InsertNextCell(pnodesD);

        // the edge
        topPaddleContactSurface->InsertNextCell4(pnodesD[5],pnodesD[4],pnodesD[1],pnodesD[2]);
      }
    }

    // Initial value before rotation is unchanged.
    double topPaddleCenterRot[3], bottomPaddleCenterRot[3];
    for(int i=0; i<3; i++){
      topPaddleCenterRot[i]=topPaddleCenter[i];
      bottomPaddleCenterRot[i]=bottomPaddleCenter[i];
    }

    // rotate paddles
    if(rotate){
      double cRot = cos(angle);
      double sRot = sin(angle);

      for(int i=0; i<bottomPaddleNodes->GetNumberOfPoints(); i++){
        double t[3];
        double tnew[3];
        bottomPaddleNodes->GetCell(i, t);
        tnew[0] = t[0];
        tnew[1] = cRot*t[1] - sRot*t[2];
        tnew[2] = sRot*t[1] + cRot*t[2];
        bottomPaddleNodes->SetTuple(i, tnew);
      }

      for(int i=0; i<topPaddleNodes->GetNumberOfCells(); i++){
        double t[3];
        double tnew[3];
        topPaddleNodes->GetCell(i, t);
        tnew[0] = t[0];
        tnew[1] = cRot*t[1] - sRot*t[2];
        tnew[2] = sRot*t[1] + cRot*t[2];
        topPaddleNodes->SetTuple(i, tnew);
      }

      topPaddleCenterRot[0] = topPaddleCenter[0];
      topPaddleCenterRot[1] = cRot*topPaddleCenter[1] - sRot*topPaddleCenter[2];
      topPaddleCenterRot[2] = sRot*topPaddleCenter[1] + cRot*topPaddleCenter[2];

      bottomPaddleCenterRot[0] = bottomPaddleCenter[0];
      bottomPaddleCenterRot[1] = cRot*bottomPaddleCenter[1] - sRot*bottomPaddleCenter[2];
      bottomPaddleCenterRot[2] = sRot*bottomPaddleCenter[1] + cRot*bottomPaddleCenter[2];
    }

*/

std::cout << "Number of points " << bottomPaddleNodes->GetNumberOfPoints() << std::endl;
 vtkNew<vtkUnstructuredGrid> ugrid;
    ugrid->SetPoints(bottomPaddleNodes);
    //ugrid->SetCells(bottomPaddleHexElements);
    std::cout << "Number of points " << ugrid->GetPoints()->GetNumberOfPoints() << std::endl;

  vtkNew<vtkDataSetMapper> ugridMapper;
  ugridMapper->SetInputData(ugrid);
  ugridMapper->Update();

vtkNew<vtkNamedColors> colors;

  vtkNew<vtkActor> ugridActor;
  ugridActor->SetMapper(ugridMapper);
  ugridActor->GetProperty()->SetColor(colors->GetColor3d("Peacock").GetData());
  ugridActor->GetProperty()->EdgeVisibilityOn();

  vtkNew<vtkRenderer> renderer;
  renderer->AddActor(ugridActor);
 renderer->SetBackground(colors->GetColor3d("Beige").GetData());
  renderer->ResetCamera();
  renderer->GetActiveCamera()->Elevation(60.0);
  renderer->GetActiveCamera()->Azimuth(30.0);
 renderer->GetActiveCamera()->Dolly(1.2);

  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(renderer);

  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  renWin->SetSize(640, 480);
  renWin->SetWindowName("UGrid");

  // interact with data
  renWin->Render();
  iren->Start();


/*vtkSmartPointer<vtkUnstructuredGridWriter> ugridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    ugridWriter->SetInputData(ugrid);
    ugridWriter->SetFile
*/
  return EXIT_SUCCESS;
}
