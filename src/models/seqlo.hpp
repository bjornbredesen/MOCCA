////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// SEQLO

class SEQLO:public sequenceClassifier{
private:
	autodelete<motifWindow> mwin;
	motifOccContainer*moc;
	motifList*motifs;
	featureSet*features;
	autodelete<featureWindow> fwin;
	autodelete<logoddsClassifier> classifier;
	SEQLO(int nf);
public:
	static SEQLO*create(motifList*motifs,featureSet*fs);
	virtual ~SEQLO(){  }
	bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls);
	bool trainFinish();
	bool flush();
	double do_applyWindow(char*buf,long long pos,int bufs);
	bool printInfo();
	bool exportAnalysisData(string path);
};

