////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// SEQSVM

class SEQSVM:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	featureSet*features;
	autodelete<featureWindow> fwin;
	autodelete<fastSVMClassifier> classifier;
	int svmtype;
	SEQSVM(int nf);
public:
	static SEQSVM*create(motifList*motifs,featureSet*fs,int _svmtype);
	virtual ~SEQSVM(){  }
	int getNPair(int ia,int ib);
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

