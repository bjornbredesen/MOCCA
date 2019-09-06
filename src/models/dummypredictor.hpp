////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// DummyPREdictor
//	Dummy of the CPREdictor (no weights are used).

class DummyPREdictor:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	autofree<double> fvec;
	DummyPREdictor(int nf);
public:
	static DummyPREdictor*create(motifList*motifs);
	virtual ~DummyPREdictor(){  }
	int getNPair(int ia,int ib);
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

