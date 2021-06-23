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
#include "seqdummy.hpp"

////////////////////////////////////////////////////////////////////////////////////
// SEQDummy

SEQDummy::SEQDummy(int nf):sequenceClassifier(nf){  }

SEQDummy*SEQDummy::create(motifList*motifs,featureSet*fs){
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
		cmdError("No features for SEQDummy classifier");
		return 0;
	}
	autodelete<SEQDummy> r(new SEQDummy(nfeatures));
	if(!r.ptr){
		outOfMemory();
		return 0;
	}
	r.ptr->mwin.ptr = _mwin.disown();
	r.ptr->fwin.ptr = _fwin.disown();
	r.ptr->moc=r.ptr->mwin.ptr->occContainer;
	r.ptr->motifs=r.ptr->mwin.ptr->motifs;
	r.ptr->features=fs;
	cout << r.ptr->nFeatures;
	return r.disown();
}

bool SEQDummy::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	return true;
}

bool SEQDummy::trainFinish(){
	return true;
}

bool SEQDummy::flush(){
	return mwin.ptr->flush();
}

double SEQDummy::do_applyWindow(char*buf,long long pos,int bufs){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return -1;
	double sum = 0.;
	for(int l=0; l<nFeatures; l++)
		sum += fvec[l];
	return sum;
}

bool SEQDummy::printInfo(){
	cout << t_indent << "SEQDummy classifier" << cmdNewline;
	return true;
}

bool SEQDummy::exportAnalysisData(string path){
	cout << "No model analysis for dummy model" << cmdNewline;
	return false;
}

