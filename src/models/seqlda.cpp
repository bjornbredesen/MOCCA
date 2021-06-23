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
#include "seqlda.hpp"

////////////////////////////////////////////////////////////////////////////////////
// SEQLDA

SEQLDA::SEQLDA(int nf):sequenceClassifier(nf){  }

SEQLDA*SEQLDA::create(motifList*motifs,featureSet*fs){
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
		cmdError("No features for SEQLDA classifier");
		return 0;
	}
	autodelete<SEQLDA> r(new SEQLDA(nfeatures));
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
	r.ptr->classifier.ptr = LDAClassifier::create(nfeatures, std::string("cls"));
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	std::vector<std::string> featureNames = fs->getInstFeatureNames(motifs);
	for(int i = 0; i < nfeatures; i++){
		r.ptr->classifier.ptr->featureNames.push_back(featureNames[i]);
	}
	return r.disown();
}

bool SEQLDA::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return false;
	classifier.ptr->addTrain(fvec,nFeatures,cls);
	return true;
}

bool SEQLDA::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

bool SEQLDA::flush(){
	return mwin.ptr->flush();
}

double SEQLDA::do_applyWindow(char*buf,long long pos,int bufs){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return -1;
	return classifier.ptr->apply(fvec,nFeatures);
}

bool SEQLDA::printInfo(){
	cout << t_indent << "SEQLDA classifier" << cmdNewline;
	classifier.ptr->printInfo((char*)"LDA");
	return true;
}

bool SEQLDA::exportAnalysisData(string path){
	cout << "Model analysis export not yet supported for SEQLDA" << cmdNewline;
	return false;
}

