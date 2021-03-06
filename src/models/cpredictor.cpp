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
#include "cpredictor.hpp"

////////////////////////////////////////////////////////////////////////////////////
// CPREdictor
//	Reimplementation of the PREdictor algorithm.

CPREdictor::CPREdictor(int nf):sequenceClassifier(nf){  }

CPREdictor*CPREdictor::create(motifList*motifs){
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
	autodelete<CPREdictor> r(new CPREdictor(nfeatures));
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
	r.ptr->classifier.ptr = logoddsClassifier::create(nfeatures);
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	if(!r.ptr->fvec.resize((size_t)nfeatures)){
		return 0;
	}
	for(int ia=0;ia<motifs->nmotifs;ia++){
		for(int ib=ia;ib<motifs->nmotifs;ib++){
			char fName[256];
			snprintf(fName, sizeof(fName), "%s vs %s", motifs->motifs[ia].name, motifs->motifs[ib].name);
			r.ptr->classifier.ptr->featureNames.push_back(std::string(fName));
		}
	}
	return r.disown();
}

int CPREdictor::getNPair(int ia,int ib){
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

bool CPREdictor::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	if(!mwin.ptr->readWindow(buf,pos,bufs)){
		return false;
	}
	int i=0;
	for(int ia=0;ia<motifs->nmotifs;ia++){
		for(int ib=ia;ib<motifs->nmotifs;ib++){
			fvec[i]=double(getNPair(ia,ib)*1000)/double(bufs);
			i++;
		}
	}
	classifier.ptr->addTrain(fvec.ptr,i,cls);
	return true;
}

bool CPREdictor::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

bool CPREdictor::flush(){
	return mwin.ptr->flush();
}

double CPREdictor::do_applyWindow(char*buf,long long pos,int bufs){
	if(!mwin.ptr->readWindow(buf,pos,bufs)){
		return false;
	}
	int nm=motifs->nmotifs;
	int i=0;
	for(int ia=0;ia<nm;ia++){
		for(int ib=ia;ib<nm;ib++){
			fvec[i]=double(getNPair(ia,ib)*1000)/double(bufs);
			i++;
		}
	}
	return classifier.ptr->apply(fvec.ptr,i);
}

bool CPREdictor::printInfo(){
	cout << t_indent << "CPREdictor classifier\n";
	classifier.ptr->printInfo((char*)"Log-odds");
	return true;
}

bool CPREdictor::exportAnalysisData(string path){
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
	cout << "Saved classifier analysis data to \"" << path << "\"\n";
	return false;
}

