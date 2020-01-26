////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// SEQRF

class SEQRF:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	featureSet*features;
	autodelete<featureWindow> fwin;
	autodelete<RFClassifier> classifier;
	SEQRF(int nf);
public:
	static SEQRF*create(motifList*motifs,featureSet*fs);
	virtual ~SEQRF(){  }
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

