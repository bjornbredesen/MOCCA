////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Configuration

char*getKernelName(int degree);

enum weightMode{
	wmPREdictor,
	wmZero,
	wmConstant,
	wmZZPF, // Zero frequency=>Zero weight, Pseudocounts Finally
	wmPPV, // Weighting by Positive Predictive Value (no need for pseudocounts)
	wmBiPPV, // Bi-directional weighting by Positive Predictive Value (no need for pseudocounts)
};

enum pairDistanceMode{
	dmCenters,
	dmBetween,
};

enum SVMKernel{
	kInvalid,
	kLinear,
	kQuadratic,
	kCubic,
	kRBF,
};

enum classifierT{
	cInvalid,
	cCPREdictor,
	cSVMMOCCA,
	cDummyPREdictor,
	cSEQSVM,
};

class config{
public:
	int randSeed;
	bool useFSM;
	bool validate;
	classifierT classifier;
	int windowSize,windowStep,windowStepTrain;
	bool allowHomoPairing, allowHeteroPairing;
	SVMKernel kernel;
	pairDistanceMode distanceMode;
	bool motifPairsCanOverlap;
	double SVM_c0,SVM_C,SVM_nu,SVM_p,SVM_gamma;
	double threshold;
	char*inFASTA,*outWig,*outCoreSequence;
	bool MOCCA_nOcc, MOCCA_GC, MOCCA_DNT;
	int svmtype;
	weightMode wmMode;
	double loBeta;
	char*CAnalysisExportPath;
	char*outSCVal;
	/*
	printInfo
		Prints out information
	*/
	void printInfo();
};

config*getConfiguration();

