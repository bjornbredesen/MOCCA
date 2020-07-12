////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Base classifiers
//	These classifiers only classify vectors of values.

class baseClassifierSmp{
public:
	autofree<double> vec;	// Training vector.
	int vecl;				// Vector length.
	seqClass*cls;			// Class.
	double clsE;			// Double value class, such as for regression.
	baseClassifierSmp(double*v,seqClass*c,double e);
	~baseClassifierSmp(){  }
	static baseClassifierSmp*create(double*v,int vl,seqClass*c);
};

class baseClassifier{
private:
public:
	deletevector<baseClassifierSmp> trainingExamples;
	int nFeatures;
	double threshold;
	bool trained;
	baseClassifier(int nf);
	virtual ~baseClassifier(){  }
	baseClassifierSmp*addTrain(double*v,int vl,seqClass*c);
	baseClassifierSmp*addTrainV(double*v,int vl,seqClass*c,double val);
	double apply(double*v,int vl);
	bool train();
	virtual bool exportAnalysisData(FILE*f, char*title, char*indent);
	virtual bool do_train() = 0;
	virtual double do_apply(double*v) = 0;
	virtual void printInfo(char*header) = 0;
};


////////////////////////////////////////////////////////////////////////////////////
// Log-odds classifier
//	Works similarly to the base classifier in the PREdictor.

class logoddsClassifier:public baseClassifier{
private:
	int nP,nN;
	autofree<double> weights, cP, cN;
	logoddsClassifier(int nf);
public:
	std::vector<std::string> featureNames;
	virtual ~logoddsClassifier(){ };
	static logoddsClassifier*create(int nf);
	double getWeight(int i);
	bool do_train();
	double do_apply(double*vec);
	void printInfo(char*header);
	bool exportAnalysisData(FILE*f, char*title, char*indent);
};


////////////////////////////////////////////////////////////////////////////////////
// Fast SVM
//	Trains SVM with LibSVM, but uses optimized classification procedures.
//	Implements feature value scaling internally.

class fastSVMClassifier:public baseClassifier{
private:
	svm_problem svmprob;		// LibSVM problem
	svm_parameter svmparam;	// LibSVM parameters
	svm_model*svmmdl;			// LibSVM model
	svm_node*svmvector;			// LibSVM hypervector
	int nFeatures;
	autofree<double> SVcoef;
	float**mSVcoef;
	autofree<double> vMin, vMax;
	autofree<float> fVec;
	fastSVMClassifier(int type,int nf,std::string _name);
	bool IaddTrain(double*fv,double cls);
public:
	std::string name;
	std::vector<std::string> featureNames;
	~fastSVMClassifier();
	static fastSVMClassifier*create(int type,int nf,std::string _name);
	void scaleVector(svm_node*n);
	void scaleVectorDoubleFloat(double*in,float*out);
	bool do_train();
	double do_apply(double*fv);
	void printInfo(char*header);
	bool exportAnalysisData(FILE*f, char*title, char*indent);
};


////////////////////////////////////////////////////////////////////////////////////
// Multi-class SVM

typedef struct{
	seqClass*cls;
	int votes;
}MultiClassClass;

typedef struct{
	MultiClassClass*clsP,*clsN;
	baseClassifier*classifier;
}MultiClassClassPair;

class MultiClassSVM{
private:
	std::vector<MultiClassClass>classes;
	std::vector<MultiClassClassPair>borders;
	deletevector<fastSVMClassifier>classifiers;
	int nFeatures;
	deletevector<baseClassifierSmp> trainingExamples;
	bool trained;
	int svmtype;
	//
	MultiClassSVM(int nf,int _svmt,std::string _name);
public:
	std::string name;
	std::vector<std::string> featureNames;
	//
	static MultiClassSVM*create(int _svmtype,int nfeatures,std::string _name);
	~MultiClassSVM(){  }
	//
	bool addTrain(double*v,int vl,seqClass*c,double cE=0);
	bool train();
	seqClass*apply(double*vec,int vecl);
	void printInfo(char*header);
	bool exportAnalysisData(FILE*f, char*title, char*indent);
};

////////////////////////////////////////////////////////////////////////////////////
// Random Forest
//	Trains RF with Ranger.

#include "../lib/ranger/src/globals.h"
#include "../lib/ranger/src/ForestClassification.h"
#include "../lib/ranger/src/ForestRegression.h"
#include "../lib/ranger/src/ForestProbability.h"
#include "../lib/ranger/src/utility.h"
#include "../lib/ranger/src/Data.h"
#include "../lib/ranger/src/DataDouble.h"

class RangerData: public ranger::DataDouble {
public:
	void setDataT(int nFeatures, std::vector<baseClassifierSmp*> dat) {
		// Set header
		for(int i=0;i<nFeatures;i++){
			variable_names.push_back("f" + std::to_string(i+1));
		}
		num_cols = variable_names.size();
		num_cols_no_snp = num_cols;
		for(baseClassifierSmp*t: dat)
			num_rows++;
		reserveMemory(1);
		bool error = false;
		// Set rows
		int row = 0;
		for(baseClassifierSmp*t: dat){
			for(int i=0;i<nFeatures;i++)
				set_x(i, row, t->vec[i], error);
			set_y(0, row, t->cls->flag ? 1. : -1., error);
			row++;
		}
	}
	void setVector(double*vec, int ncol){
		num_rows = 1;
		num_cols = ncol;
		num_cols_no_snp = num_cols;
		reserveMemory(0);
		bool error = false;
		for(int i=0;i<ncol;i++)
			set_x(i, 0, vec[i], error);
	}
};

class RangerRandomForest: public ranger::ForestProbability {
public:
	void train(){
		grow();
		computePredictionError();
		dependent_variable_names.push_back(std::string("target"));
	}
	double predictVec(double*vec){
		RangerData*rd = (RangerData*)data.get();
		int ncol = rd->getNumCols();
		rd->setVector(vec, ncol);
		num_samples = 1;
		ranger::equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
		predict();
		double ppos = 0.;
		double pneg = 0.;
		for(int icls = 0; icls < class_values.size(); icls++){
			if (class_values[icls] == 1.0)
				ppos = predictions[0][0][icls];
			else if (class_values[icls] == -1.0)
				pneg = predictions[0][0][icls];
		}
		return ppos - pneg;
	}
};

class RFClassifier:public baseClassifier{
private:
	RFClassifier(int nf);
	autodelete<RangerRandomForest> rf;
	std::string name;
public:
	std::vector<std::string> featureNames;
	virtual ~RFClassifier(){ };
	static RFClassifier*create(int nf, std::string name);
	bool do_train();
	double do_apply(double*vec);
	void printInfo(char*header);
};

////////////////////////////////////////////////////////////////////////////////////
// Linear Discriminant Analysis

#ifdef USE_SHOGUN
#include <shogun/base/init.h>

#include <shogun/lib/config.h>
#include <shogun/labels/MulticlassLabels.h>
#include <shogun/classifier/LDA.h>
#include <shogun/features/DenseFeatures.h>
#include <shogun/io/SGIO.h>
#include <shogun/lib/common.h>
#include <shogun/features/DataGenerator.h>

class LDAClassifier:public baseClassifier{
private:
	LDAClassifier(int nf);
	autodelete<shogun::CLDA> LDA;
	std::string name;
public:
	std::vector<std::string> featureNames;
	virtual ~LDAClassifier(){ };
	static LDAClassifier*create(int nf, std::string name);
	bool do_train();
	double do_apply(double*vec);
	void printInfo(char*header);
};

#endif

////////////////////////////////////////////////////////////////////////////////////
// Quadratic Discriminant Analysis

#ifdef USE_SHOGUN
#include <shogun/classifier/AveragedPerceptron.h>

class PerceptronClassifier:public baseClassifier{
private:
	PerceptronClassifier(int nf);
	autodelete<shogun::CAveragedPerceptron> perceptron;
	std::string name;
public:
	std::vector<std::string> featureNames;
	virtual ~PerceptronClassifier(){ };
	static PerceptronClassifier*create(int nf, std::string name);
	bool do_train();
	double do_apply(double*vec);
	void printInfo(char*header);
};

#endif

