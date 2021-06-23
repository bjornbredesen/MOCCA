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
#include "seqlo.hpp"

////////////////////////////////////////////////////////////////////////////////////
// SEQLO

SEQLO::SEQLO(int nf):sequenceClassifier(nf){  }

SEQLO*SEQLO::create(motifList*motifs,featureSet*fs){
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
		cmdError("No features for SEQLO classifier");
		return 0;
	}
	autodelete<SEQLO> r(new SEQLO(nfeatures));
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
	r.ptr->classifier.ptr = logoddsClassifier::create(nfeatures);
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	std::vector<std::string> featureNames = fs->getInstFeatureNames(motifs);
	for(int i = 0; i < nfeatures; i++){
		r.ptr->classifier.ptr->featureNames.push_back(featureNames[i]);
	}
	return r.disown();
}

bool SEQLO::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return false;
	classifier.ptr->addTrain(fvec,nFeatures,cls);
	return true;
}

bool SEQLO::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

bool SEQLO::flush(){
	return mwin.ptr->flush();
}

double SEQLO::do_applyWindow(char*buf,long long pos,int bufs){
	double*fvec=fwin.ptr->extractFeatures(buf,pos,bufs,true);
	if(!fvec)return -1;
	return classifier.ptr->apply(fvec,nFeatures);
}

bool SEQLO::printInfo(){
	cout << t_indent << "SEQLO classifier" << cmdNewline;
	classifier.ptr->printInfo((char*)"Log-odds");
	return true;
}

bool SEQLO::exportAnalysisData(string path){
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
	cout << "Model analysis export not yet supported for SEQLO" << cmdNewline;
	return false;
}

