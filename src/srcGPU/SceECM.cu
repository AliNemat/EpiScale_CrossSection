#include "SceECM.h"
#include "SceCells.h"
//# define debugModeECM 
// task: frequency of plotting the ECM should be imported. Right now is given explicitly
// bending stiffness is given inside the code. It should be given as in input from a txt file.
//isInitPhase bool variable is not active anymore.
__constant__ double sceInterCell_ECM[5]; 
//__constant__ double wLCPara_ECM[4]; 
__constant__ double restLenECMAdhSpringGPU  ;  
__constant__ double maxLenECMAdhSpringGPU ;
__constant__ double kAdhECMGPU ;
__constant__ double stiffnessECMBasalGPU ;
__constant__ double stiffnessECMBCGPU ;
__constant__ double stiffnessECMPeripGPU ;
__constant__ double lknotECMBasalGPU ;
__constant__ double lknotECMBCGPU ;
__constant__ double lknotECMPeripGPU ;


namespace patch{
	template <typename  T> std::string to_string (const T& n) 
	{
	std:: ostringstream stm ; 
	stm << n ; 
	return stm.str() ; 
	}
}



__device__
void DefineECMStiffnessAndLknot ( EType nodeType, double & stiffness, double & sponLen) {

	if (nodeType==excm) {
		stiffness=stiffnessECMBasalGPU ;  //* stiffness ; 
		sponLen=lknotECMBasalGPU  ; 
	}
	if (nodeType==perip) {
		stiffness=stiffnessECMPeripGPU ; //*stiffness ; 
		sponLen=lknotECMPeripGPU ; // 0.1 ; 
	}

	if (nodeType==bc2) {
		stiffness=stiffnessECMBCGPU;  // *stiffness ; 
		sponLen=lknotECMBCGPU ;// _sponLen ; 
	}

}

__device__
double calMorse_ECM(const double& linkLength ) {

	double forceValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterCell_ECM[0] / sceInterCell_ECM[2]
				* exp(-linkLength / sceInterCell_ECM[2])
				+ sceInterCell_ECM[1] / sceInterCell_ECM[3]
						* exp(-linkLength / sceInterCell_ECM[3]);
//		if (forceValue > 0) {
//			forceValue = 0;
//		}
	}

	return (forceValue) ; 
}

__device__
double calMorseEnergy_ECM(const double& linkLength ) {

	double energyValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		energyValue = 0;
	} else {
		energyValue = sceInterCell_ECM[0]* exp(-linkLength / sceInterCell_ECM[2])
				    - sceInterCell_ECM[1]* exp(-linkLength / sceInterCell_ECM[3]);
	}

	return (energyValue) ; 
}

/*
__device__
double calWLC_ECM(const double& linkLength ) {

	double x=linkLength/wLCPara_ECM[0] ;
	return (wLCPara_ECM[1]*( 6*x+ ( x*x*(3.0-2*x))/( (1-x)*(1-x) ) )
	       -wLCPara_ECM[2]/pow(linkLength,wLCPara_ECM[3]) ) ; 	
}
*/
__device__
bool IsValidAdhPair(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU  && dist < maxLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}
__device__
bool IsValidAdhPairForNotInitPhase(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}


__device__
double  CalAdhECM(const double& dist ) {
	return (kAdhECMGPU*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}

__device__
double  CalAdhEnergy(const double& dist ) {
	return (0.5*kAdhECMGPU*(dist-restLenECMAdhSpringGPU)*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}


EType SceECM:: ConvertStringToEType(string eNodeRead) {
	if (eNodeRead=="perip")  {
		return perip ; 
	}
	else if (eNodeRead=="bc2") {
		return bc2 ; 
	}
	else if (eNodeRead=="excm") {
		return excm ; 
	}
	else {
		cout << "Error in defining type of external nodes" << endl ; 
		return excm ;// To just return something to avoid compiler complain 
	}
} 


void SceECM::Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM, int freqPlotData, string uniqueSymbolOutput) {

maxAllNodePerCell=maxAllNodePerCellECM ; 
maxMembrNodePerCell= maxMembrNodePerCellECM ; 
maxTotalNodes=maxTotalNodesECM ; //Ali 
this->freqPlotData=freqPlotData ; 
this->uniqueSymbolOutput=uniqueSymbolOutput ; 


std::fstream readCoord_ECM ;
std::fstream readInput_ECM ;
int numberNodes_ECM ; 
double tmpPosX_ECM,tmpPosY_ECM ; 
vector<double> posXIni_ECM,posYIni_ECM ;
vector <EType> eNodeVec ;  
readCoord_ECM.open("./resources/coordinate_ECM21.txt") ;
if (readCoord_ECM.is_open()) {
	cout << "ECM coordinates file opened successfully" <<endl ; 
}
else {
	cout << "ECM coordinates file is not opened successfully" << endl ; 
}


string inputInfoText ; 
string eNodeRead ; 
readCoord_ECM>>numberNodes_ECM ;
for (int i=0 ; i<numberNodes_ECM ; i++){
	readCoord_ECM>>tmpPosX_ECM>>tmpPosY_ECM>>eNodeRead  ; 
	posXIni_ECM.push_back(tmpPosX_ECM) ; 
	posYIni_ECM.push_back(tmpPosY_ECM) ; 
	EType eNode=ConvertStringToEType(eNodeRead) ; 
	eNodeVec.push_back(eNode) ; 
}

readInput_ECM.open("./resources/ECM_input.txt") ;
if (readInput_ECM.is_open()) {
	cout << "ECM Mech input opened successfully" <<endl ; 
}
else {
	cout << "ECM Mech input is not opened successfully" << endl ; 
}

 
 readInput_ECM>> inputInfoText ; 
 for (int i=0 ; i<5; i++) {
 	readInput_ECM>> mechPara_ECM.sceInterCellCPU_ECM[i] ; //=39.0 ; 
 }
 
// readInput_ECM>>restLenECMSpring ;
// readInput_ECM>>eCMLinSpringStiff ;    
 readInput_ECM>>restLenECMAdhSpring  ;  
 readInput_ECM>>maxLenECMAdhSpring ;
 readInput_ECM>>kAdhECM ;
 //for ( int i=0 ; i<4 ; i++) {
//	readInput_ECM>>mechPara_ECM.wLCParaCPU_ECM[i] ;
// }    



std::fstream secondInput_ECM ; 
std:: string secondInputInfo ;  //dummy 
std::string secondInputFileName = "./resources/ECM_" + uniqueSymbolOutput + "input.cfg";

secondInput_ECM.open(secondInputFileName.c_str()) ;
if (secondInput_ECM.is_open()) {
	cout << "Second ECM Mech input opened successfully" <<endl ; 
}
else {
	cout << "Second ECM Mech input is not opened successfully" << endl ; 
}

 secondInput_ECM>>secondInputInfo ;  // just for information no use in the code
 secondInput_ECM>>stiffnessECMBasal ;
 secondInput_ECM>>stiffnessECMBC ;
 secondInput_ECM>>stiffnessECMPerip ;
 secondInput_ECM>>lknotECMBasal ;
 secondInput_ECM>>lknotECMBC ;
 secondInput_ECM>>lknotECMPerip ;

 cout <<" stiffness of ECM at the basal side is="<<stiffnessECMBasal <<endl ;   

 cout <<" stiffness of ECM at boundary is="<<stiffnessECMBC<<endl ; 

 cout <<" stiffness of ECM peripodial side is="<<stiffnessECMPerip<<endl ;  
 cout <<" rest len basal ECM is="<<lknotECMBasal<<endl ;
 cout <<" rest len boundary ECM is= "<<lknotECMBC<<endl ;
 cout << "rest len peripodial ECM is=" <<lknotECMPerip <<endl ; 




cout<< "number of ECM nodes is"<< numberNodes_ECM <<endl ; 
//for (int i=0 ; i<posXIni_ECM.size() ; i++){
//	cout << "ECM nodes read in cpu"<<endl; 
//	cout << posXIni_ECM[i] <<", "<<posYIni_ECM[i]<<", " << eNodeVec[i] <<endl;  ; 
//}




 for (int i=0 ; i<5; i++) {
	cout <<"Morse parameter number"<<i<<" is " <<mechPara_ECM.sceInterCellCPU_ECM[i]<<endl ; 

} 

 //cout <<"rest length of ECM spring is "<<restLenECMSpring<<endl ;   

// cout <<"ECM spring stiffness is "<<eCMLinSpringStiff<<endl ; 

 cout <<"ECM Membrane neutral adhesion length is "<<restLenECMAdhSpring<<endl ;  
 cout <<"ECM Membrane max adhesion length is "<<maxLenECMAdhSpring<<endl ;
 cout <<"ECM Membrane adhesion stiffness is "<<kAdhECM<<endl ;
 cout << "ECM only applies adhesvie force" << endl ; 

//for ( int i=0 ; i<4 ; i++) {
//	cout<<"wLC parameter "<< i << " is "<<mechPara_ECM.wLCParaCPU_ECM[i]<<endl ;  ;
//}    

cudaMemcpyToSymbol(sceInterCell_ECM,mechPara_ECM.sceInterCellCPU_ECM
			,5*sizeof(double)); 
//cudaMemcpyToSymbol(wLCPara_ECM,mechPara_ECM.wLCParaCPU_ECM
//			,4*sizeof(double)); 

cudaMemcpyToSymbol(restLenECMAdhSpringGPU, &restLenECMAdhSpring,sizeof(double));

cudaMemcpyToSymbol(maxLenECMAdhSpringGPU, &maxLenECMAdhSpring,sizeof(double));
cudaMemcpyToSymbol(kAdhECMGPU, &kAdhECM,sizeof(double));

cudaMemcpyToSymbol(stiffnessECMPeripGPU, &stiffnessECMPerip,sizeof(double));
cudaMemcpyToSymbol(stiffnessECMBCGPU, &stiffnessECMBC,sizeof(double));
cudaMemcpyToSymbol(stiffnessECMBasalGPU, &stiffnessECMBasal,sizeof(double));

cudaMemcpyToSymbol(lknotECMPeripGPU, & lknotECMPerip,sizeof(double));
cudaMemcpyToSymbol(lknotECMBCGPU, & lknotECMBC,sizeof(double));
cudaMemcpyToSymbol(lknotECMBasalGPU, & lknotECMBasal,sizeof(double));


counter=100000 ; //large number
lastPrintECM=1000000 ; // large number 
outputFrameECM=0 ; 
numNodesECM= numberNodes_ECM ; //(eCMMaxX-eCMMinX)/eCMMinDist ; 



indexECM.resize(numNodesECM,0) ;
peripORexcm.resize(numNodesECM,perip) ;

nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,0.0) ;


cellNeighborId.resize(numNodesECM,-1) ;

stiffLevel.resize(numNodesECM) ;
sponLen.resize(numNodesECM) ;

linSpringForceECMX.resize(numNodesECM,0.0); 
linSpringForceECMY.resize(numNodesECM,0.0); 
linSpringAvgTension.resize(numNodesECM,0.0); 
linSpringEnergy.resize(numNodesECM,0.0); 
morseEnergy.resize(numNodesECM,0.0); 
adhEnergy.resize(numNodesECM,0.0); 


bendSpringForceECMX.resize(numNodesECM,0.0); 
bendSpringForceECMY.resize(numNodesECM,0.0); 
 
memMorseForceECMX.resize(numNodesECM,0.0); 
memMorseForceECMY.resize(numNodesECM,0.0);
 
fBendCenterX.resize(numNodesECM,0.0); 
fBendCenterY.resize(numNodesECM,0.0); 
fBendLeftX.resize(numNodesECM,0.0); 
fBendLeftY.resize(numNodesECM,0.0); 
fBendRightX.resize(numNodesECM,0.0); 
fBendRightY.resize(numNodesECM,0.0); 
 
totalForceECMX.resize(numNodesECM,0.0); 
totalForceECMY.resize(numNodesECM,0.0);

//memNodeType.resize(maxTotalNodes,notAssigned1) ; 

thrust::sequence (indexECM.begin(),indexECM.begin()+numNodesECM);
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+604,excm );
//thrust::fill (peripORexcm.begin()+1166,peripORexcm.end(),excm );
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+int(numNodesECM/2),excm );
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+int(numNodesECM/2),excm );
 
thrust::copy(posXIni_ECM.begin(),posXIni_ECM.end(),nodeECMLocX.begin()) ; 
thrust::copy(posYIni_ECM.begin(),posYIni_ECM.end(),nodeECMLocY.begin()) ; 
thrust::copy(eNodeVec.begin(),eNodeVec.end(),peripORexcm.begin()) ;

//cout << "GPU level initial coordinates and type of external nodes are: " << endl ; 
//for (int i=0;  i<nodeECMLocX.size() ;  i++) {
//	cout<< nodeECMLocX[i]<<", "<<nodeECMLocY[i]<<", "<<peripORexcm[i] << endl; 
//}

PrintECM(0.0) ;
std::string cSVFileName = "./ECMFolder/EnergyExport_" + uniqueSymbolOutput + ".CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str());

			EnergyExport <<"Time,"<<"TotalMorseEnergyECM," << "TotalAdhEnergyECM,"<<"TotalLinSpringEnergy,"<<"TotalEnergy, " <<"TotalEnergyDerivative"<< std::endl;

} //initilaization function finished






void SceECM:: ApplyECMConstrain(int currentActiveCellCount, int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_Coef, bool cellPolar, bool subCellPolar, bool isInitPhase){   
#ifdef debugModeECM 
	cudaEvent_t start1, start2, start3, start4, start5, start6, start7, start8, stop;
	float elapsedTime1, elapsedTime2, elapsedTime3, elapsedTime4, elapsedTime5, elapsedTime6,  elapsedTime7 , elapsedTime8 ; 
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&start3);
	cudaEventCreate(&start4);
	cudaEventCreate(&start5);
	cudaEventCreate(&start6);
	cudaEventCreate(&start7);
	cudaEventCreate(&start8);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start1, 0);
#endif
thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
//isBasalMemNode.resize(totalNodeCountForActiveCellsECM,false) ;
adhPairECM_Cell.resize(totalNodeCountForActiveCellsECM,-1) ;
morseEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
adhEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
thrust::copy(nodesPointerECM->getInfoVecs().nodeLocX.begin(),nodesPointerECM->getInfoVecs().nodeLocX.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocX.begin()) ; 
thrust::copy(nodesPointerECM->getInfoVecs().nodeLocY.begin(),nodesPointerECM->getInfoVecs().nodeLocY.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocY.begin()) ; 
//cout << " max all node per cell in ECM module is " << maxAllNodePerCell << endl ; 
//cout<< "max membrane node per cell in ECM module is " << maxMembrNodePerCell<< endl ; 
//cout<< "I am in ECM module and dt is: " << dt << endl ; 
#ifdef debugModeECM
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif

double eCMBendStiff=6.0 ; // need to be an input

//if (cellPolar) {eCMLinSpringStiff=100 ; }
//if (subCellPolar) {eCMLinSpringStiff=100 ; }

//cout << "test to make sure ECM class reads cells class variables "<< cellsPointerECM->getCellInfoVecs().basalLocX[0] << endl ; 
//cout << "test to make sure ECM class reads cells class variables "<< cellsPointerECM->getCellInfoVecs().basalLocY[0] << endl ; 




double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 

EType* peripORexcmAddr= thrust::raw_pointer_cast (
			&peripORexcm[0]) ; 

// move the nodes of epithelial cells 
//// find the closest ECM node to each each cell //

 int numCells = cellsPointerECM->getCellInfoVecs().basalLocX.size() ;

counter ++ ; 
if (counter>=1000000 || curTime<(100*dt)) {
	counter=0 ; 
	thrust:: transform (
		thrust::make_zip_iterator (
					thrust:: make_tuple (
						cellsPointerECM->getCellInfoVecs().basalLocX.begin(),
						cellsPointerECM->getCellInfoVecs().basalLocY.begin())),
		thrust::make_zip_iterator (
					thrust:: make_tuple (
					 	cellsPointerECM->getCellInfoVecs().basalLocX.begin(),
                     	cellsPointerECM->getCellInfoVecs().basalLocY.begin()))+numCells, 
	    cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
		FindECMNeighborPerCell(nodeECMLocXAddr,nodeECMLocYAddr,numNodesECM ));

 	double * basalCellLocXAddr= thrust::raw_pointer_cast ( & ( cellsPointerECM->getCellInfoVecs().basalLocX[0]) ) ; 
 	double * basalCellLocYAddr= thrust::raw_pointer_cast ( & ( cellsPointerECM->getCellInfoVecs().basalLocY[0]) ) ;
	// int numCells = cellsPointerECM->getCellInfoVecs().basalLocX.size() ;
 	//cout << " Number of cells in ECM class is equal to " << numCells << endl; 
	thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 nodeECMLocX.begin(),
                     nodeECMLocY.begin()))+numNodesECM,
	    cellNeighborId.begin(),
		FindCellNeighborPerECMNode(basalCellLocXAddr,basalCellLocYAddr, numCells));
}




thrust::counting_iterator<int> iBegin2(0) ; 
//////////////////////////////////////////
 thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_permutation_iterator(
						cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor2(
												maxAllNodePerCell))),
					make_transform_iterator (iBegin,
							DivideFunctor2(maxAllNodePerCell)),
					make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					nodesPointerECM->getInfoVecs().nodeLocX.begin(),
					nodesPointerECM->getInfoVecs().nodeLocY.begin(), 
					nodesPointerECM->getInfoVecs().nodeIsActive.begin(),
					nodesPointerECM->getInfoVecs().memNodeType1.begin()
					)), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_permutation_iterator(
						cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor2(
												maxAllNodePerCell))),
					make_transform_iterator (iBegin,
							DivideFunctor2(maxAllNodePerCell)),
					 make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					 nodesPointerECM->getInfoVecs().nodeLocX.begin(),
                     nodesPointerECM->getInfoVecs().nodeLocY.begin(),
					 nodesPointerECM->getInfoVecs().nodeIsActive.begin(),
					 nodesPointerECM->getInfoVecs().memNodeType1.begin()
					 ))+totalNodeCountForActiveCellsECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodesPointerECM->getInfoVecs().nodeLocX.begin(),
					nodesPointerECM->getInfoVecs().nodeLocY.begin(),
					adhPairECM_Cell.begin(),
					morseEnergyCell.begin(),
					adhEnergyCell.begin())),
				MoveNodes2_Cell(nodeECMLocXAddr,nodeECMLocYAddr,maxMembrNodePerCell,numNodesECM,dt,Damp_Coef,isInitPhase,peripORexcmAddr,curTime,currentActiveCellCount));

#ifdef debugModeECM
	cudaEventRecord(start3, 0);
	cudaEventSynchronize(start3);
	cudaEventElapsedTime(&elapsedTime2, start2, start3);
#endif
/* To reduce computational cost
energyECM.totalMorseEnergyCellECM = thrust::reduce( morseEnergyCell.begin(),morseEnergyCell.begin()+totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() ); 
energyECM.totalAdhEnergyCellECM   = thrust::reduce( adhEnergyCell.begin()  ,adhEnergyCell.begin()  +totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() );
*/
bool* nodeIsActiveAddr= thrust::raw_pointer_cast (
			& (nodesPointerECM->getInfoVecs().nodeIsActive[0])) ; 

int * adhPairECM_CellAddr= thrust::raw_pointer_cast (
			&adhPairECM_Cell[0]) ; 


thrust:: transform (peripORexcm.begin(), peripORexcm.begin()+numNodesECM,
                    thrust::make_zip_iterator (thrust::make_tuple (
											                        stiffLevel.begin(),
																	sponLen.begin()    ))
					                           ,MechProp(isInitPhase));


double* stiffLevelAddr=thrust::raw_pointer_cast (
			&stiffLevel[0]) ; 

double*  sponLenAddr =thrust::raw_pointer_cast (
			&sponLen[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					linSpringAvgTension.begin(),
					linSpringEnergy.begin())),
				LinSpringForceECM(numNodesECM,nodeECMLocXAddr,nodeECMLocYAddr,stiffLevelAddr,sponLenAddr));

//////////////////////////////////// find the closest Cell to each ECM node ///////////



///////////////////////////////////

//cout << " I am after FindCellNeighbor functor" << endl ; 







#ifdef debugModeECM
	cudaEventRecord(start4, 0);
	cudaEventSynchronize(start4);
	cudaEventElapsedTime(&elapsedTime3, start3, start4);
#endif

//energyECM.totalLinSpringEnergyECM = 0.5 * ( thrust::reduce( linSpringEnergy.begin(),linSpringEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() )); 
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					fBendCenterX.begin(),
					fBendCenterY.begin(),
					fBendLeftX.begin(),
					fBendLeftY.begin(),
					fBendRightX.begin(),
					fBendRightY.begin())),
				CalBendECM(nodeECMLocXAddr,nodeECMLocYAddr,numNodesECM,eCMBendStiff));

#ifdef debugModeECM
	cudaEventRecord(start5, 0);
	cudaEventSynchronize(start5);
	cudaEventElapsedTime(&elapsedTime4, start4, start5);
#endif
double* fBendLeftXAddr= thrust::raw_pointer_cast (
			&fBendLeftX[0]) ; 
double* fBendLeftYAddr= thrust::raw_pointer_cast (
			&fBendLeftY[0]) ; 
double* fBendRightXAddr= thrust::raw_pointer_cast (
			&fBendRightX[0]) ; 
double* fBendRightYAddr= thrust::raw_pointer_cast (
			&fBendRightY[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					fBendCenterX.begin(),
					fBendCenterY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 fBendCenterX.begin(),
                                         fBendCenterY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin())),
				SumBendForce(fBendLeftXAddr,fBendLeftYAddr,fBendRightXAddr,fBendRightYAddr,numNodesECM));

#ifdef debugModeECM
	cudaEventRecord(start6, 0);
	cudaEventSynchronize(start6);
	cudaEventElapsedTime(&elapsedTime5, start5, start6);
#endif
//to make sure it is based on the distance used for action force calculation.
double* nodeCellLocXAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocX[0]) ; 
double* nodeCellLocYAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocY[0]) ;
 

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					cellNeighborId.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                     nodeECMLocY.begin(),
					 cellNeighborId.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin(),
					morseEnergy.begin(),
					adhEnergy.begin())),
				MorseAndAdhForceECM(numCells,maxAllNodePerCell,maxMembrNodePerCell,nodeCellLocXAddr,nodeCellLocYAddr,nodeIsActiveAddr,adhPairECM_CellAddr));

//cout << " I am after MorseAndAdhForceECM functor" << endl ; 

#ifdef debugModeECM
	cudaEventRecord(start7, 0);
	cudaEventSynchronize(start7);
	cudaEventElapsedTime(&elapsedTime6, start6, start7);
#endif

/* To reduce computational cost 
energyECM.totalMorseEnergyECMCell = thrust::reduce( morseEnergy.begin(),morseEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() ); 
energyECM.totalAdhEnergyECMCell   = thrust::reduce( adhEnergy.begin()  ,adhEnergy.begin()  +numNodesECM,(double) 0.0, thrust::plus<double>() );
*/

double dummy=0.0 ;

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					totalForceECMX.begin(),
					totalForceECMY.begin())),
				TotalECMForceCompute(dummy));

#ifdef debugModeECM
	cudaEventRecord(start8, 0);
	cudaEventSynchronize(start8);
	cudaEventElapsedTime(&elapsedTime7, start7, start8);
#endif



double Damp_CoefECM=Damp_Coef ; 
double Damp_CoefPerip=Damp_Coef ; 
// move the nodes of ECM 
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin(),
					peripORexcm.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin(),
					peripORexcm.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
				MoveNodeECM(dt,Damp_CoefECM,Damp_CoefPerip,numNodesECM,curTime));

/* To reduce computational cost
cout << "total Morse energy for cell-ECM is= "<< energyECM.totalMorseEnergyCellECM << endl ; 
cout << "total Morse energy for ECM-cell  is= "<< energyECM.totalMorseEnergyECMCell << endl ;
cout << "total adhesion energy for cell-ECM is= "<<  energyECM.totalAdhEnergyCellECM << endl ; 
cout << "total adhesion energy for ECM-cell  is= "<< energyECM.totalAdhEnergyECMCell << endl ; 
//assert (abs (energyECM.totalMorseEnergyCellECM-energyECM.totalMorseEnergyECMCell)<1.0) ;
//assert (abs (energyECM.totalAdhEnergyCellECM-  energyECM.totalAdhEnergyECMCell)  <1.0) ;


if (  (abs (energyECM.totalMorseEnergyCellECM-energyECM.totalMorseEnergyECMCell)>1.0) || 
	  (abs (energyECM.totalAdhEnergyCellECM-  energyECM.totalAdhEnergyECMCell)  >1.0)
   ) {

	cout << "Warning: Action and reaction forces in the ECM do not match each other" << endl ; 
}
*/ 

# ifdef debugModeECM 
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime8, start8, stop);
	std::cout << "time 1 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime1 << endl ; 
	std::cout << "time 2 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime2 << endl ; 
	std::cout << "time 3 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime3 << endl ; 
	std::cout << "time 4 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime4 << endl ; 
	std::cout << "time 5 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime5 << endl ; 
	std::cout << "time 6 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime6 << endl ; 
	std::cout << "time 7 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime7 << endl ; 
	std::cout << "time 8 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime8 << endl ; 
#endif

PrintECM(curTime); 

}

void  SceECM:: PrintECM(double curTime) {
		lastPrintECM=lastPrintECM+1 ; 
        if (lastPrintECM>=freqPlotData) { 
			outputFrameECM++ ; 
			lastPrintECM=0 ;

			// First ECM output file for paraview //
			std::string vtkFileName = "./ECMFolder/ECM_" + uniqueSymbolOutput +patch::to_string(outputFrameECM-1) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << nodeECMLocX.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
				<< 0.0 << std::endl;
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()<< " " << 3 *nodeECMLocX.size()<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << indexECM[i] << " "
				<< indexECM[i+1] << std::endl;
			}
			ECMOut << 2 << " " << indexECM[nodeECMLocX.size()-1] << " "<< indexECM[0] << std::endl; //last point to the first point


			ECMOut << "CELL_TYPES " << nodeECMLocX.size()<< endl;
			for (uint i = 0; i < nodeECMLocX.size() ; i++) {
				ECMOut << "3" << endl;
			}
			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
			ECMOut << "SCALARS Avg_Tension " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<linSpringAvgTension[i] <<endl ; 
			}
			
			ECMOut << "SCALARS Node_Type " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<peripORexcm[i] <<endl ; 
			}

			ECMOut.close();



			// second output file for curvature estimation //
			std::string txtFileName = "./ECMFolder/ECMLocationExport_" + uniqueSymbolOutput+ patch::to_string(outputFrameECM-1) + ".txt";
			ofstream ECMLocationExport ;
			ECMLocationExport.open(txtFileName.c_str());
			//ECMExport << "ECM pouch coordinates" << std::endl;

			for (uint i = 0; i < nodeECMLocX.size(); i++) {
		//		if (peripORexcm[i]==excm) {
					ECMLocationExport<< nodeECMLocX[i] << " " << nodeECMLocY[i] << " " << 0.0 << " "<< peripORexcm[i]<<std::endl;
		//		}
			}

			//ECMExport << "ECM lumen side coordinates" << std::endl;
		//	for (uint i = 0; i < nodeECMLocX.size(); i++) {
		//		if (peripORexcm[i]==perip) {
		//			ECMLocationExport << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
		//			<< 0.0 << std::endl;
		//		}
		//	}

			ECMLocationExport.close();
		

			//Third write file for ECM
			txtFileName = "./ECMFolder/ECMTensionExport_" + uniqueSymbolOutput+ patch::to_string(outputFrameECM-1) + ".txt";
			ofstream ECMTensionExport ;
			ECMTensionExport.open(txtFileName.c_str());

			for (uint i = 0; i < nodeECMLocX.size(); i++) {
					ECMTensionExport<< linSpringAvgTension[i]<< " " << peripORexcm[i]<< std::endl;
				}

			ECMTensionExport.close();
			///
			//Fourth write file for ECM
			energyECM.totalEnergyECMOld=energyECM.totalEnergyECM ; 
			energyECM.totalEnergyECM=       energyECM.totalMorseEnergyECMCell
			                         +      energyECM.totalAdhEnergyECMCell
								     +      energyECM.totalLinSpringEnergyECM ;
									 

			std::string cSVFileName = "./ECMFolder/EnergyExport_" + uniqueSymbolOutput+ ".CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str(),ofstream::app);
			
			//EnergyExport <<"totalMorseEnergyCell " << "totalAdhEnergyCell "<<  "totalMorseEnergy "<<"totalAdhEnergy "<< "totalLinSpringEnergy " << std::endl;
			EnergyExport <<curTime<<","<<energyECM.totalMorseEnergyECMCell << "," << energyECM.totalAdhEnergyECMCell<< "," << energyECM.totalLinSpringEnergyECM <<"," << energyECM.totalEnergyECM <<","<<energyECM.totalEnergyPrimeECM <<std::endl;


		}

}




