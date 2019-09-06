////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
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

