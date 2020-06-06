////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

struct prediction{
	int start, end, center;
	double score;
	int mstart, mend;
	prediction(int _start, int _end, double _score, int _mstart = -1, int _mend = -1){
		start = _start;
		end = _end;
		center = (start + end) / 2;
		score = _score;
		mstart = _mstart;
		mend = _mend;
	}
};

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
	bool applyFASTA(std::string inpath,std::string outpath);
	bool predictCoreSequence(std::string inpath, std::string outpath);
	/*
	calibrateThresholdGenomewidePrecision
		Calibrates the classifier threshold for a desired genome-wide precision.
	*/
	bool calibrateThresholdGenomewidePrecision(seqList*calpos,double wantPrecision);
	/*
	predictGenomewideFASTA
		Applies the classifier to a genome FASTA file and writes scores to output GFF- and Wig-files.
	*/
	bool predictGenomewideFASTA(std::string inFASTAPath,std::string outGFFPath,std::string outWigPath);
	//
	virtual bool trainWindow(char*buf,long long pos,int bufs,seqClass*cls) = 0;
	virtual bool trainFinish() = 0;
	virtual double do_applyWindow(char*buf,long long pos,int bufs) = 0;
	virtual bool flush() = 0;
	virtual bool printInfo() = 0;
	virtual bool exportAnalysisData(string path) = 0;
	virtual vector<prediction> predictWindow(char*buf,long long pos,int bufs, corePredictionModeT cpm);
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

