////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// CPREdictor
//	Reimplementation of the PREdictor algorithm.

class CPREdictor:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	autofree<double> fvec;
	autodelete<logoddsClassifier> classifier;
	CPREdictor(int nf);
public:
	static CPREdictor*create(motifList*motifs);
	virtual ~CPREdictor(){  }
	int getNPair(int ia,int ib);
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

