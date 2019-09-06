////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

#include "sequencelist.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Validation

typedef struct{
	double score;
	seqClass*cls;
}validationPair;

/*
getROCAUCxVP
	Gets ROC AUC(x) for a validation pair set.
*/
double getROCAUCxVP(validationPair*vp,int nvp,double x);

/*
getPRCAUCxVP
	Gets PRC AUC(x) for a validation pair set.
*/
double getPRCAUCxVP(validationPair*vp,int nvp,double x);

bool printValidationMeasures(validationPair*vp,int nvp,double threshold);

bool saveVPairTable(char*outpath,validationPair*vp,int nvp);

