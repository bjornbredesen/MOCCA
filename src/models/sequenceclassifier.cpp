////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "../common.hpp"
#include "../lib/libsvm-3.17/svm.h"
#include "../aux.hpp"
#include "../config.hpp"
#include "../validation.hpp"
#include "../motifs.hpp"
#include "../sequences.hpp"
#include "../sequencelist.hpp"
#include "baseclassifier.hpp"
#include "sequenceclassifier.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Sequence classifiers
//	These classifiers classify sequences.

sequenceClassifier::sequenceClassifier(int nf){
	cfg=getConfiguration();
	threshold=cfg->threshold;
	trained=false;
	nFeatures=nf;
}

bool sequenceClassifier::train(seqList*sl){
	if(!sl){
		cmdError("Training sequence list is a null-pointer.");
		return false;
	}
	if(trained){
		cmdError("Already trained.");
		return false;
	}
	cmdTask task((char*)"Training");
	cmdTask::refresh();
	seqListSeq*sls=sl->seq;
	for(int l=0;l<sl->nseq;l++,sls++){
		if(!flush()){
			return false;
		}
		if(sls->trainMode==train_Full){
			// Full, normalized
			if(!trainWindow(sls->buf,0,sls->bufs,sls->cls)){
				return false;
			}
		}else{
			// Windows
			autodelete<seqStreamBuffer> ssb((seqStreamBuffer*)0);
			autodelete<seqStreamWindow> ssw((seqStreamWindow*)0);
			if(!createBufferWindow(sls->buf,sls->bufs,ssb.ptr,ssw.ptr,cfg->windowSize,cfg->windowStepTrain)){
				return false;
			}
			char*rb;
			int rbn;
			for(long i=0;(rbn=ssw.ptr->get(rb));i+=cfg->windowStepTrain){
				if(!trainWindow(rb,i,rbn,sls->cls)){
					return false;
				}
			}
		}
	}
	if(!trainFinish()){
		return false;
	}
	trained=true;
	return true;
}

double sequenceClassifier::getSequenceScore(seqListSeq*sls){
	if(!flush()){
		return -1.0;
	}
	// Windows
	autodelete<seqStreamBuffer> ssb((seqStreamBuffer*)0);
	autodelete<seqStreamWindow> ssw((seqStreamWindow*)0);
	if(!createBufferWindow(sls->buf,sls->bufs,ssb.ptr,ssw.ptr,cfg->windowSize,cfg->windowStep)){
		return -1.0;
	}
	char*rb;
	int rbn;
	double maxScore=0;
	for(long i=0;(rbn=ssw.ptr->get(rb));i+=cfg->windowStep){
		double v=applyWindow(rb,i,rbn);
		if(!i||v>maxScore)maxScore=v;
	}
	return maxScore;
}

bool sequenceClassifier::getValidationTable(seqList*sl,validationPair*&vp,int&nvp){
	if(!sl){
		cmdError("Validation sequence list is a null-pointer.");
		return false;
	}
	if(!trained){
		cmdError("Tried to validate untrained classifier.");
		return false;
	}
	seqListSeq*sls=sl->seq;
	autofree<validationPair>vPairs(sl->nseq);
	if(!vPairs.ptr){
		outOfMemory();
		return false;
	}
	for(int l=0;l<sl->nseq;l++,sls++){
		double maxScore=getSequenceScore(sls);
		vPairs.ptr[l].score=maxScore+threshold; // Store max. score without threshold
		vPairs.ptr[l].cls=sls->cls;
	}
	vp=vPairs.disown();
	nvp=sl->nseq;
	return true;
}

double sequenceClassifier::applyWindow(char*buf,long long pos,int bufs){
	if(!trained){
		cmdWarning("Tried to apply untrained classifier.");
		return -1;
	}
	double r=do_applyWindow(buf,pos,bufs);
	return r-threshold;
}

bool sequenceClassifier::applyFASTA(std::string inpath, std::string outpath){
	FILE*fout=0;
	if(!inpath.length() || !outpath.length())return false;
	timer mainTimer((char*)"Sequence scoring");
	cmdTask task((char*)"Sequence scoring");
	seqStreamFastaBatch*ssfb=seqStreamFastaBatch::load((char*)inpath.c_str());
	if(!ssfb){
		return false;
	}
	long long bptotal=0;
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb->getBlock());){
		bptotal+=ssfbblk->getLength();
	}
	delete ssfb;
	fout=fopen((char*)outpath.c_str(),"wb");
	if(!fout){
		cmdError("Could not open file for writing.");
		return false;
	}
	ssfb=seqStreamFastaBatch::load((char*)inpath.c_str());
	if(!ssfb){
		fclose(fout);
		return false;
	}
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb->getBlock());){
		seqStreamWindow*ssw=seqStreamWindow::create(ssfbblk,cfg->windowSize,cfg->windowStep);
		if(!ssw){
			fclose(fout);
			delete ssfb;
			return false;
		}
		fprintf(fout,"fixedStep start=1 step=%d span=%d chrom=%s\n",cfg->windowStep,cfg->windowSize,ssfbblk->getName());
		// Apply in windows.
		char*rb;
		int rbn;
		long nextSi=0;
		double cvalue;
		flush();
		for(long i=0;(rbn=ssw->get(rb));i+=cfg->windowStep){
			if(i>=nextSi){
				nextSi+=50000;
				task.setPercent((double(i)/double(bptotal))*100.0);
			}
			cvalue=applyWindow(rb,i,rbn);
			fprintf(fout,"%lf\n",cvalue+threshold);
		}
		delete ssw;
	}
	delete ssfb;
	fclose(fout);
	return true;
}

bool sequenceClassifier::predictCoreSequence(std::string inpath, std::string outpath){
	if(!inpath.length() || !outpath.length())return false;
	FILE*fout=0;
	fout=fopen((char*)outpath.c_str(),"wb");
	if(!fout){
		cmdError("Could not open file for writing.");
		return false;
	}
	autodelete<seqStreamFastaBatch> ssfb(seqStreamFastaBatch::load((char*)inpath.c_str()));
	if(!ssfb.ptr){
		fclose(fout);
		return false;
	}
	cout << "Predicted sequence core:\n";
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb.ptr->getBlock());){
		char*sname=ssfbblk->getName();
		if(!sname)sname=(char*)"Unnamed";
		autofree<char> buf((char*)0);
		int bufs = (int)ssfbblk->buffer(buf.ptr);
		// Get maximum window score with current window size
		int mwA = 0;
		int mwB = 0;
		double mwScore = 0.0;
		const int N_WSIZES = 8;
		int wSizes[N_WSIZES] = {500, 600, 750, 1000, 1500, 2000, 2500, 3000};
		for(int wsi = 0; wsi < N_WSIZES; wsi++){
			int cwSize = wSizes[wsi];
			//int cwStep = 5;
			int cwStep = 50;
			// Apply in windows
			autodelete<seqStreamBuffer> ssb((seqStreamBuffer*)0);
			autodelete<seqStreamWindow> ssw((seqStreamWindow*)0);
			if(!createBufferWindow(buf, bufs, ssb.ptr, ssw.ptr, cwSize, cwStep)){
				fclose(fout);
				return false;
			}
			char*rb;
			int rbn;
			double cvalue;
			flush();
			for(int i = 0; (rbn = ssw.ptr->get(rb)); i += cwStep){
				cvalue = applyWindow(rb,(long int)i,rbn);
				if(cvalue > mwScore){
					mwScore = cvalue;
					mwA = i;
					mwB = i + rbn;
				}
			}
		}
		fprintf(fout,"%s\t%d\t%d\t%f\n",sname,mwA,mwB,mwScore);
		cout << " - " << sname << ": " << mwA << ".." << mwB << " (" << (mwB-mwA) << " / " << bufs << " bp) - Score: " << mwScore << "\n";
	}
	fclose(fout);
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
// Training sequences

trainingSequence::trainingSequence(char*s,int l,seqClass*c,bool o){
	seq = s;
	length = l;
	cls = c;
	own = o;
}

trainingSequence::~trainingSequence(){
	if(own&&seq){
		free(seq);
	}
}

trainingSequence*trainingSequence::createSequenceClone(char*s,int l,seqClass*c){
	if(!s||l<=0||!c){
		cmdError("Invalid training sequence specification");
		return 0;
	}
	char*sc = (char*)malloc(sizeof(char)*l);
	if(!sc){
		outOfMemory();
		return 0;
	}
	memcpy(sc,s,sizeof(char)*l);
	trainingSequence*r = new trainingSequence(sc,l,c,true);
	if(!r){
		outOfMemory();
		free(sc);
		return 0;
	}
	r->seq = sc;
	return r;
}

