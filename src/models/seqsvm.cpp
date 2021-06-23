////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "../common.hpp"
#include "../lib/libsvm-3.17/svm.h"
#include "../vaux.hpp"
#include "../config.hpp"
#include "../validation.hpp"
#include "../motifs.hpp"
#include "../sequences.hpp"
#include "../sequencelist.hpp"
#include "baseclassifier.hpp"
#include "sequenceclassifier.hpp"
#include "features.hpp"
#include "seqsvm.hpp"

////////////////////////////////////////////////////////////////////////////////////
// SEQSVM

SEQSVM::SEQSVM(int nf):sequenceClassifier(nf){  }

SEQSVM*SEQSVM::create(motifList*motifs,featureSet*fs,int _svmtype){
	if(!motifs){cmdError("Null-pointer argument.");return 0;}
	autodelete<motifWindow> _mwin;
	autodelete<featureWindow> _fwin;
	_mwin.ptr = motifWindow::create(motifs);
	if(!_mwin.ptr){
		return 0;
	}
	_fwin.ptr = featureWindow::create(_mwin.ptr, fs);
	if(!_fwin.ptr){
		return 0;
	}
	int nfeatures = _fwin.ptr->getNFeatures();
	if(!nfeatures){
		cmdError("No features for SEQSVM classifier");
		return 0;
	}
	autodelete<SEQSVM> r(new SEQSVM(nfeatures));
	if(!r.ptr){
		outOfMemory();
		return 0;
	}
	r.ptr->mwin.ptr = _mwin.disown();
	r.ptr->fwin.ptr = _fwin.disown();
	r.ptr->svmtype=_svmtype;
	r.ptr->moc=r.ptr->mwin.ptr->occContainer;
	r.ptr->motifs=r.ptr->mwin.ptr->motifs;
	r.ptr->features=fs;
	cout << r.ptr->nFeatures;
	r.ptr->classifier.ptr = fastSVMClassifier::create(_svmtype,nfeatures,std::string("SEQSVM"));
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	std::vector<std::string> featureNames = fs->getInstFeatureNames(motifs);
	for(int i = 0; i < nfeatures; i++){
		r.ptr->classifier.ptr->featureNames.push_back(featureNames[i]);
	}
	return r.disown();
}

bool SEQSVM::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return false;
	classifier.ptr->addTrain(fvec,nFeatures,cls);
	return true;
}

bool SEQSVM::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

bool SEQSVM::flush(){
	return mwin.ptr->flush();
}

double SEQSVM::do_applyWindow(char*buf,long long pos,int bufs){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return -1;
	return classifier.ptr->apply(fvec,nFeatures);
}

bool SEQSVM::printInfo(){
	cout << t_indent << "SEQSVM classifier" << cmdNewline;
	classifier.ptr->printInfo((char*)"SVM");
	return true;
}

bool SEQSVM::exportAnalysisData(string path){
	/*
	FILE*f=fopen(path.c_str(),"wb");
	int i=0;
	fprintf(f,"Motif pair\tWeight\n");
	for(int ia=0;ia<motifs->nmotifs;ia++){
		for(int ib=ia;ib<motifs->nmotifs;ib++){
			double w = classifier.ptr->getWeight(i);
			//fvec[i]=double(getNPair(ia,ib)*1000)/double(bufs);
			fprintf(f,"%s:%s\t\%f\n",motifs->motifs[ia].name,motifs->motifs[ib].name,w);
			i++;
		}
	}
	fclose(f);
	cout << "Saved classifier analysis data to \"" << path << "\"\n";*/
	cout << "Model analysis export not yet supported for SEQSVM" << cmdNewline;
	return false;
}

