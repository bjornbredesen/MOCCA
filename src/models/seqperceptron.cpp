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
#include "seqperceptron.hpp"

////////////////////////////////////////////////////////////////////////////////////
// SEQPerceptron

SEQPerceptron::SEQPerceptron(int nf):sequenceClassifier(nf){  }

SEQPerceptron*SEQPerceptron::create(motifList*motifs,featureSet*fs){
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
		cmdError("No features for SEQPerceptron classifier");
		return 0;
	}
	autodelete<SEQPerceptron> r(new SEQPerceptron(nfeatures));
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
	r.ptr->classifier.ptr = PerceptronClassifier::create(nfeatures, std::string("cls"));
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	std::vector<std::string> featureNames = fs->getInstFeatureNames(motifs);
	for(int i = 0; i < nfeatures; i++){
		r.ptr->classifier.ptr->featureNames.push_back(featureNames[i]);
	}
	return r.disown();
}

bool SEQPerceptron::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return false;
	classifier.ptr->addTrain(fvec,nFeatures,cls);
	return true;
}

bool SEQPerceptron::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

bool SEQPerceptron::flush(){
	return mwin.ptr->flush();
}

double SEQPerceptron::do_applyWindow(char*buf,long long pos,int bufs){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return -1;
	return classifier.ptr->apply(fvec,nFeatures);
}

bool SEQPerceptron::printInfo(){
	cout << t_indent << "SEQPerceptron classifier" << cmdNewline;
	classifier.ptr->printInfo((char*)"Perceptron");
	return true;
}

bool SEQPerceptron::exportAnalysisData(string path){
	cout << "Model analysis export not yet supported for SEQPerceptron" << cmdNewline;
	return false;
}

