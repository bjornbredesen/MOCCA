////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
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
#include "dummypredictor.hpp"

////////////////////////////////////////////////////////////////////////////////////
// DummyPREdictor
//	Dummy of the CPREdictor (no weights are used).

DummyPREdictor::DummyPREdictor(int nf):sequenceClassifier(nf){  }

DummyPREdictor*DummyPREdictor::create(motifList*motifs){
	if(!motifs){cmdError("Null-pointer argument.");return 0;}
	int nfeatures=0;
	for(int ia=0;ia<motifs->nmotifs;ia++){
		for(int ib=ia;ib<motifs->nmotifs;ib++){
			nfeatures++;
		}
	}
	if(!nfeatures){
		cmdError("No features for PREdictor classifier");
		return 0;
	}
	autodelete<DummyPREdictor> r(new DummyPREdictor(nfeatures));
	if(!r.ptr){
		outOfMemory();
		return 0;
	}
	r.ptr->mwin.ptr = motifWindow::create(motifs);
	if(!r.ptr->mwin.ptr){
		return 0;
	}
	r.ptr->moc=r.ptr->mwin.ptr->occContainer;
	r.ptr->motifs=r.ptr->mwin.ptr->motifs;
	if(!r.ptr->fvec.resize((size_t)nfeatures)){
		return 0;
	}
	return r.disown();
}

int DummyPREdictor::getNPair(int ia,int ib){
	int nPairCut=219;
	int iv=0;
	if(ia==ib){
		if(!cfg->allowHomoPairing)return 0;
		motifOcc*o=moc->getFirst(ia);
		while(o){
			motifOcc*o2=o;
			while(o2){
				if(isMotifPair(o,o2,0,nPairCut))iv++;
				o2=moc->getNextSame(o2);
			}
			o=moc->getNextSame(o);
		}
	}else{
		if(!cfg->allowHeteroPairing)return 0;
		motifOcc*o=moc->getFirst(ia);
		while(o){
			motifOcc*o2=moc->getFirst(ib);
			while(o2){
				if(isMotifPair(o,o2,0,nPairCut))iv++;
				o2=moc->getNextSame(o2);
			}
			o=moc->getNextSame(o);
		}
	}
	return iv;
}

bool DummyPREdictor::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	return true;
}

bool DummyPREdictor::trainFinish(){
	return true;
}

bool DummyPREdictor::flush(){
	return mwin.ptr->flush();
}

double DummyPREdictor::do_applyWindow(char*buf,long long pos,int bufs){
	if(!mwin.ptr->readWindow(buf,pos,bufs)){
		return false;
	}
	int nm=motifs->nmotifs;
	int i=0;
	double ret = 0.0;
	for(int ia=0;ia<nm;ia++){
		for(int ib=ia;ib<nm;ib++){
			ret += double(getNPair(ia,ib)*1000)/double(bufs);
			i++;
		}
	}
	return ret;
}
bool DummyPREdictor::printInfo(){
	cout << t_indent << "DummyPREdictor classifier\n";
	cout << t_indent << "# features = " << nFeatures << "\n";
	return true;
}
bool DummyPREdictor::exportAnalysisData(string path){
	cmdError("Export for classifier analysis not implemented for DummyPREdictor\n");
	return false;
}

