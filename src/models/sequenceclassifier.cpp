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
	int bptotal=0;
	for(int l=0;l<sl->nseq;l++,sls++)
		bptotal+=sls->bufs;
	sls = sl->seq;
	cmdTask task((char*)"Scoring sequences");
	int cbp = 0;
	int nextbp = 0;
	for(int l=0;l<sl->nseq;l++,sls++){
		double maxScore=getSequenceScore(sls);
		cbp+=sls->bufs;
		if(cbp>=nextbp){
			task.setPercent((double(cbp)/double(bptotal))*100.0);
			nextbp+=10000;
		}
		vPairs.ptr[l].score=maxScore+threshold; // Store max. score without threshold
		vPairs.ptr[l].cls=sls->cls;
	}
	vp=vPairs.disown();
	nvp=sl->nseq;
	cmdTask::wipe();
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

bool sequenceClassifier::predictGenomewideFASTA(std::string inFASTAPath, std::string outGFFPath, std::string outWigPath){
	if(!inFASTAPath.length())return false;
	timer mainTimer((char*)"Genome-wide prediction");
	cmdTask task((char*)"Genome-wide prediction");
	seqStreamFastaBatch*ssfb=seqStreamFastaBatch::load((char*)inFASTAPath.c_str());
	if(!ssfb){
		return false;
	}
	long long bptotal=0;
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb->getBlock());){
		bptotal+=ssfbblk->getLength();
	}
	delete ssfb;
	ssfb=seqStreamFastaBatch::load((char*)inFASTAPath.c_str());
	if(!ssfb){
		return false;
	}
	ofstream ofGFF;
	if(outGFFPath.length())
		ofGFF.open(outGFFPath);
	ofstream ofWig;
	if(outWigPath.length())
		ofWig.open(outWigPath);
	int nPredictions = 0;
	int cit = 0;
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb->getBlock());){
		std::string streamName = std::string(ssfbblk->getName());
		size_t ti = streamName.find(" ");
		std::string chromName = ti == std::string::npos ? streamName : streamName.substr(0, ti);
		cmdTask task((char*)chromName.c_str());
		if(ofWig.is_open())
			ofWig << "fixedStep chrom=" << chromName << " start=1 step=" << cfg->windowStep << " span=" << cfg->windowSize << "\n";
		//
		int pStart = -1;
		int pEnd = -1;
		double pScore = 0.;
		//
		seqStreamWindow*ssw=seqStreamWindow::create(ssfbblk,cfg->windowSize,cfg->windowStep);
		if(!ssw){
			delete ssfb;
			return false;
		}
		// Apply in windows.
		char*rb;
		int rbn;
		long nextSi=0;
		double cvalue;
		flush();
		for(long i=0;(rbn=ssw->get(rb));i+=cfg->windowStep){
			cit += rbn - (cfg->windowSize-cfg->windowStep);
			if(cit>=nextSi){
				nextSi+=50000;
				task.setPercent((double(cit)/double(bptotal))*100.0);
			}
			cvalue=applyWindow(rb,i,rbn);
			if(ofWig.is_open())
				ofWig << (cvalue+threshold) << "\n";
			if(cvalue>=0.){
				if(pEnd==-1){
					pStart=i;
					pScore=cvalue+threshold;
				}else{
					pScore=max(cvalue+threshold,pScore);
				}
				pEnd=i+rbn;
			}else if(pEnd != -1 && i > pEnd){
				cmdTask::wipe();
				cout << t_indent << "Predicted: " << chromName << ":" << pStart << ".." << pEnd << " - score: " << pScore << "\n";
				cmdTask::refresh();
				if(ofGFF.is_open())
					ofGFF << chromName << "\tMOCCA\tPrediction\t" << pStart << "\t" << pEnd << "\t" << pScore << "\t.\t.\t1\n";
				pEnd = -1;
				nPredictions++;
			}
		}
		if(pEnd != -1){
			cmdTask::wipe();
			cout << t_indent << "Predicted: " << chromName << ":" << pStart << ".." << pEnd << " - score: " << pScore << "\n";
			cmdTask::refresh();
			if(ofGFF.is_open())
				ofGFF << chromName << "\tMOCCA\tPrediction\t" << pStart << "\t" << pEnd << "\t" << pScore << "\t.\t.\t1\n";
			nPredictions++;
		}
		delete ssw;
	}
	delete ssfb;
	cmdTask::wipe();
	cout << t_indent << "Made " << nPredictions << " predictions genome-wide\n";
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

bool sequenceClassifier::calibrateThresholdGenomewidePrecision(seqList*calpos,double wantPrecision){
	threshold=0.;
	int nvp=0;
	autofree<validationPair> vp((validationPair*)0);
	timer mainTimer((char*)"Threshold calibration");
	cmdTask task((char*)"Calibrating threshold");
	{
		if(!sequenceClassifier::getValidationTable(calpos,vp.ptr,nvp))
			return false;
	}
	cmdTask taskp((char*)"Handling scores");
	vector<validationPair> svp;
	int nP = 0, nN = 0;
	for(int x=0;x<nvp;x++){
		svp.push_back(vp[x]);
		if(vp[x].cls->flag) nP++;
		else nN++;
	}
	sort(svp.begin(),svp.end(),
	[](const validationPair a,const validationPair b){
		return a.score < b.score;
	});
	int TP = nP, FP = nN;
	double aPrec = -1., bPrec = -1.;
	double aThr = 0., bThr = 0.;
	for(auto& cvp: svp){
		aPrec = bPrec;
		bPrec = (double)(TP) / (double)(TP+FP);
		aThr = bThr;
		bThr = cvp.score;
		if(bPrec>=wantPrecision)
			break;
		if(cvp.cls->flag){
			TP--;
		}else{
			FP--;
		}
	}
	threshold=aThr;
	if(aPrec!=bPrec)
		threshold = aThr + (wantPrecision-aPrec)*((bThr-aThr)/(bPrec-aPrec));
	cmdTask::wipe();
	cout << "Precision: " << wantPrecision << "\n";
	cout << "Calibrated threshold: " << threshold << "\n";
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

