////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Motif classifier

class RFMotifOccClassifier{
private:
	config*cfg;
	autodelete<RFClassifier> classifier;
	int motifInd;
	int nfeatures;
	autofree<double> features;
	motifList*motifs;
	autodelete<MotifClassifier_featureSet> featureSet;
	RFMotifOccClassifier(int mi,motifList*ml);
public:
	~RFMotifOccClassifier(){  }
	static RFMotifOccClassifier*create(int mi,motifList*ml);
	bool trainOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs,seqClass*_cls);
	bool trainFinish();
	double applyOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs);
	void printInfo();
	bool exportAnalysisData(FILE*f);
};

////////////////////////////////////////////////////////////////////////////////////
// RF-MOCCA

class RFMOCCA:public sequenceClassifier{
private:
	deletevector<RFMotifOccClassifier> subcls;
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	RFMOCCA(int nf);
	//
	autodelete<logoddsClassifier> classifier;
	autofree<double>fvec;
	deletevector<trainingSequence> trainSeq;
public:
	static RFMOCCA*create(motifList*motifs);
	virtual ~RFMOCCA(){  }
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool flush();
	bool printInfo();
	bool exportAnalysisData(string path);
	virtual vector<prediction> predictWindow(char*buf,long long pos,int bufs, corePredictionModeT cpm);
};

