////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// SEQLDA

class SEQLDA:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	featureSet*features;
	autodelete<featureWindow> fwin;
	autodelete<LDAClassifier> classifier;
	SEQLDA(int nf);
public:
	static SEQLDA*create(motifList*motifs,featureSet*fs);
	virtual ~SEQLDA(){  }
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

