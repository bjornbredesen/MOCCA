////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Sequence classifiers

class sequenceClassifier{
private:
public:
	config*cfg;
	int nFeatures;
	double threshold;
	bool trained;
	sequenceClassifier(int nf);
	virtual ~sequenceClassifier(){ }
	bool train(seqList*sl);
	double getSequenceScore(seqListSeq*sls);
	bool getValidationTable(seqList*sl,validationPair*&vp,int&nvp);
	double applyWindow(char*buf,long long pos,int bufs);
	//
	/*
	applyFASTA
		Applies the classifier to a FASTA file and writes scores to a Wig file.
	*/
	bool applyFASTA(char*inpath,char*outpath);
	bool predictCoreSequence(char*inpath,char*outpath);
	//
	virtual bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls) = 0;
	virtual bool trainFinish() = 0;
	virtual double do_applyWindow(char*buf,long long pos,int bufs) = 0;
	virtual bool flush() = 0;
	virtual bool printInfo() = 0;
	virtual bool exportAnalysisData(string path) = 0;
};

////////////////////////////////////////////////////////////////////////////////////
// Training sequences

class trainingSequence{
public:
	int length;
	char*seq;
	seqClass*cls;
	bool own;
	trainingSequence(char*s,int l,seqClass*c,bool o=false);
	~trainingSequence();
	static trainingSequence*createSequenceClone(char*s,int l,seqClass*c);
};

