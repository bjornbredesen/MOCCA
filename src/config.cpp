////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "./lib/libsvm-3.17/svm.h"
#include "vaux.hpp"
#include "sequences.hpp"
#include "sequencelist.hpp"
#include "config.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Configuration

char*SVMKernelName[]={
	(char*)"Invalid",
	(char*)"Linear",
	(char*)"Quadratic",
	(char*)"Cubic",
	(char*)"RBF",
};

char*getKernelName(int degree){
	return degree < 0 || degree > 4 ? SVMKernelName[0] : SVMKernelName[degree];
}

config _config = config {
	0,
	true,
	true,
	cSVMMOCCA,
	500,100,250,
	true, true,
	kLinear,
	dmCenters,
	true,
	0.0,1.0,0.5,0.995,1.0,
	0.0,
	500, 0, 8,
	//bool MOCCA_nOcc, MOCCA_GC, MOCCA_DNT;
	false, false, false,
	EPSILON_SVR,
	wmPREdictor,
	0.00000001,
	"",
	"","","",
	"",
	"","","",
	-1.,
	4, // Background model order
	cpmNone,
	false
};

config*getConfiguration(){
	return &_config;
}

void config::printInfo(){
	cmdSection("Settings");
	cout << t_indent << "Random seed: " << randSeed << cmdNewline;
	cout << t_indent << "IUPAC parsing: " << (useFSM?(char*)"Finite State Machine":(char*)"Naive") << cmdNewline;
	cout << t_indent << "Standard threshold: " << threshold << cmdNewline;
	string cv="None";
	cout << t_indent << "Window size: " << windowSize << cmdNewline
		<< t_indent << "Window stepping: " << windowStep << cmdNewline
		<< t_indent << "Window stepping (training): " << windowStepTrain << cmdNewline;
	cout << t_indent << "Distance mode: ";
	switch(distanceMode){
		case dmCenters:cout << "Centers";break;
		case dmBetween:cout << "Between";break;
	}
	if(motifPairsCanOverlap)cout << " (can overlap)\n";
	else cout << " (can not overlap)\n";
	cout << t_indent << "Core prediction: ";
	switch(corePredictionMode){
		case cpmNone:cout << "Disabled";break;
		case cpmContinuous:cout << "Continuous";break;
		case cpmMotifs:cout << "Motifs";break;
		case cpmMotifsStrong:cout << "Motifs (strong)";break;
	}
	if(corePredictionMax) cout << " - maximum";
	cout << cmdNewline;
	cout << t_indent << "Log-odds mode: ";
	switch(wmMode){
		case wmPREdictor:cout << "PREdictor";break;
		case wmZero:cout << "Zero for missing frequencies";break;
		case wmConstant:cout << "Constant pseudocounts (beta=" << loBeta << ")";break;
		case wmPPV:cout << "Positive Predictive Value weights";break;
		case wmBiPPV:cout << "Bi-directional Positive Predictive Value weights";break;
		default:{}
	}
	cout << cmdNewline;
	if(classifier == cSVMMOCCA || classifier == cSEQSVM)
		cout << t_indent << "SVM kernel: " << getKernelName(kernel) << cmdNewline;
	if(genomeFASTAPath.length() > 0) cout << t_indent << "Genome: " << genomeFASTAPath << cmdNewline;
}

