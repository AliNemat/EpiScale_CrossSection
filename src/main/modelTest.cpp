//============================================================================
// Name        : MeshGen.cpp
// Author      : Wenzhao Sun
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "MeshGen.h"
#include "commonData.h"
#include "CellInitHelper.h"
#include <vector>
#include "SimulationDomainGPU.h"

using namespace std;

GlobalConfigVars globalConfigVars;

void generateStringInputs(std::string &loadMeshInput,
		std::string &animationInput, std::string &animationFolder,
		std::vector<std::string> &boundaryMeshFileNames) {
	std::string meshLocation =
			globalConfigVars.getConfigValue("MeshLocation").toString();
	std::string meshName =
			globalConfigVars.getConfigValue("MeshName").toString();
	std::string meshExtention =
			globalConfigVars.getConfigValue("MeshExtention").toString();
	loadMeshInput = meshLocation + meshName + meshExtention;

	animationFolder =
			globalConfigVars.getConfigValue("AnimationFolder").toString();
	animationInput = animationFolder
			+ globalConfigVars.getConfigValue("AnimationName").toString();

	std::string boundaryMeshLocation = globalConfigVars.getConfigValue(
			"BoundaryMeshLocation").toString();
	std::string boundaryMeshName = globalConfigVars.getConfigValue(
			"BoundaryMeshName").toString();
	std::string boundaryMeshExtention = globalConfigVars.getConfigValue(
			"BoundaryMeshExtention").toString();
	std::string boundaryMeshInput = boundaryMeshLocation + boundaryMeshName
			+ boundaryMeshExtention;
	boundaryMeshFileNames.push_back(boundaryMeshInput);
}

int main() {
	srand(time(NULL));
	ConfigParser parser;
	std::string configFileName = "./resources/modelTest.cfg";
	globalConfigVars = parser.parseConfigFile(configFileName);

	std::string loadMeshInput;
	std::string animationInput;
	std::vector<std::string> boundaryMeshFileNames;
	std::string animationFolder;
	generateStringInputs(loadMeshInput, animationInput, animationFolder,
			boundaryMeshFileNames);

	double SimulationTotalTime = globalConfigVars.getConfigValue(
			"SimulationTotalTime").toDouble();
	double SimulationTimeStep = globalConfigVars.getConfigValue(
			"SimulationTimeStep").toDouble();
	int TotalNumOfOutputFrames = globalConfigVars.getConfigValue(
			"TotalNumOfOutputFrames").toInt();

	const double simulationTime = SimulationTotalTime;
	const double dt = SimulationTimeStep;
	const int numOfTimeSteps = simulationTime / dt;
	const int totalNumOfOutputFrame = TotalNumOfOutputFrames;
	const int outputAnimationAuxVarible = numOfTimeSteps
			/ totalNumOfOutputFrame;

	/*

	 CellInitHelper initHelper;

	 SimulationDomainGPU simuDomain;
	 SimulationInitData initData = initHelper.generateDiskInput(loadMeshInput);
	 simuDomain.initialize_V2(initData);

	 simuDomain.checkIfAllDataFieldsValid();

	 */

	//GEOMETRY::MeshGen meshGen;
	//std::vector<GEOMETRY::Point2D> points = meshGen.createBdryPointsOnCircle(7,
	//		8);
	//Criteria criteria(0.125, 2.0);
	//GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2D(points,
	//		GEOMETRY::MeshGen::default_list_of_seeds, criteria);
	//std::string testString = "./resources/BdryData_unit_test.txt";
	//std::string testString = "./resources/beakBdryInput.txt";
	//GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DFromFile(testString);
	//mesh.outputVtkFile("modelTest.vtk");
	//double SimulationTimeStep = globalConfigVars.getConfigValue(
	//		"SimulationTimeStep").toDouble();
	CellInitHelper initHelper;

	RawDataInput rawInput = initHelper.generateRawInput_stab(loadMeshInput);

	SimulationInitData simuData = initHelper.initInputsV2(rawInput);

	SimulationDomainGPU simuDomain;

	std::vector<CVector> stabilizedCenters = simuDomain.stablizeCellCenters(
			simuData);

	RawDataInput rawInput2 = initHelper.generateRawInputWithProfile(
			stabilizedCenters);

	SimulationInitData simuData2 = initHelper.initInputsV2(rawInput2);

	simuDomain.initialize(simuData2);

	AnimationCriteria aniCri;
	aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();
	aniCri.isStressMap = false;

	for (int i = 0; i <= numOfTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % outputAnimationAuxVarible == 0) {
			simuDomain.outputVtkFilesWithColor(animationInput, i, aniCri);
		}
		simuDomain.runAllLogic(SimulationTimeStep);
	}

	return 0;
}