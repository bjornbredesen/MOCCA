////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Motif classifier

enum MotifClassifier_feature{
	MF_Invalid,
	MF_nOcc,
	MF_GC,
	MF_DNT,
};

typedef struct{
	MotifClassifier_feature f;
	int ia,ib,ic;
	double da,db;
}MotifClassifier_featureI;

class MotifClassifier_featureSet{
private:
	MotifClassifier_featureSet(motifList*_motifs);
public:
	~MotifClassifier_featureSet();
	static MotifClassifier_featureSet*create(motifList*_motifs);
	motifList*motifs;
	MotifClassifier_featureI*features;
	std::vector<std::string> featureNames;
	int nfeatures;
	MotifClassifier_featureI*addFeature(MotifClassifier_feature f,int ia,int ib,int ic,double da,double db);
	bool addFeatures(MotifClassifier_feature f,int ia,int ib,int ic,double da,double db);
	bool getFeatures(double*out,motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs);
};

class MotifOccClassifier{
private:
	config*cfg;
	autodelete<MultiClassSVM> classifier;
	int motifInd;
	int nfeatures;
	autofree<double> features;
	motifList*motifs;
	autodelete<MotifClassifier_featureSet> featureSet;
	MotifOccClassifier(int mi,motifList*ml);
public:
	~MotifOccClassifier(){  }
	static MotifOccClassifier*create(int svmtype,int mi,motifList*ml);
	bool trainOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs,seqClass*_cls);
	bool trainFinish();
	double applyOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs);
	void printInfo();
	bool exportAnalysisData(FILE*f);
};

////////////////////////////////////////////////////////////////////////////////////
// SVM-MOCCA

class SVMMOCCA:public sequenceClassifier{
private:
	deletevector<MotifOccClassifier> subcls;
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	int svmtype;
	SVMMOCCA(int nf,int _svmtype);
	//
	autodelete<logoddsClassifier> classifier;
	autofree<double>fvec;
	deletevector<trainingSequence> trainSeq;
public:
	static SVMMOCCA*create(motifList*motifs,int _svmtype=C_SVC);
	virtual ~SVMMOCCA(){  }
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool flush();
	bool printInfo();
	bool exportAnalysisData(string path);
	virtual vector<prediction> predictWindow(char*buf,long long pos,int bufs, corePredictionModeT cpm);
};

