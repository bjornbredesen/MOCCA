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
#include "svmmocca.hpp"
#include "rfmocca.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Motif classifier

RFMotifOccClassifier::RFMotifOccClassifier(int mi,motifList*ml){
	cfg=getConfiguration();
	motifInd=mi;
	motifs=ml;
	nfeatures=0;
}

RFMotifOccClassifier*RFMotifOccClassifier::create(int mi,motifList*ml){
	autodelete<RFMotifOccClassifier> r(new RFMotifOccClassifier(mi,ml));
	if(!r.ptr){
		outOfMemory();
		return 0;
	}
	r.ptr->featureSet.ptr=MotifClassifier_featureSet::create(ml);
	if(!r.ptr->featureSet.ptr){
		return 0;
	}
	config*cfg=r.ptr->cfg;
	// Add features
	if(cfg->MOCCA_nOcc)
		for(int l=0;l<ml->nmotifs;l++)
			r.ptr->featureSet.ptr->addFeature(MF_nOcc,l,250,0,0,0);
	if(cfg->MOCCA_GC)
		r.ptr->featureSet.ptr->addFeature(MF_GC,0,250,0,0,0);
	if(cfg->MOCCA_DNT)
		r.ptr->featureSet.ptr->addFeatures(MF_DNT,0,250,0,0,0);
	//
	r.ptr->nfeatures=r.ptr->featureSet.ptr->nfeatures;
	if(!r.ptr->nfeatures){
		cmdError("No features for motif classifier");
		return 0;
	}
	r.ptr->features.resize((size_t)r.ptr->nfeatures);
	if(!r.ptr->features.ptr){
		return 0;
	}
	r.ptr->features.fill((size_t)r.ptr->nfeatures,0);
	char cName[256];
	snprintf(cName, sizeof(cName), "Motif occurrence classifier - Motif: %s", ml->motifs[mi].name);
	r.ptr->classifier.ptr=RFClassifier::create(r.ptr->nfeatures,std::string(cName));
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	for(int l = 0; l < r.ptr->featureSet.ptr->nfeatures; l++){
		r.ptr->classifier.ptr->featureNames.push_back(std::string(r.ptr->featureSet.ptr->featureNames[l]));
	}
	return r.disown();
}

bool RFMotifOccClassifier::trainOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs,seqClass*_cls){
	if(!featureSet.ptr->getFeatures(features.ptr,o,moc,wpos,buf,bufs)){
		return false;
	}
	if(!classifier.ptr->addTrain(features.ptr,nfeatures,_cls)){
		return false;
	}
	return true;
}

bool RFMotifOccClassifier::trainFinish(){
	if(!classifier.ptr->train())return false;
	return true;
}

double RFMotifOccClassifier::applyOcc(motifOcc*o,motifOccContainer*moc,long long wpos,char*buf,int bufs){
	if(!featureSet.ptr->getFeatures(features,o,moc,wpos,buf,bufs)){
		return false;
	}
	return classifier.ptr->apply(features.ptr,nfeatures);
	//seqClass*r=classifier.ptr->apply(features.ptr,nfeatures);
	//if(!r)return 0;
	//return r->flag?1.0:-1.0;
}

void RFMotifOccClassifier::printInfo(){
	classifier.ptr->printInfo((char*)"Motif occurrence classifier");
}

bool RFMotifOccClassifier::exportAnalysisData(FILE*f){
	return classifier.ptr->exportAnalysisData(f,(char*)"Motif occurrence classifier",(char*)" - ");
}

////////////////////////////////////////////////////////////////////////////////////
// RF-MOCCA

RFMOCCA::RFMOCCA(int nf):sequenceClassifier(nf){
	// TODO Initialize
}

RFMOCCA*RFMOCCA::create(motifList*motifs){
	if(!motifs){cmdError("Null-pointer argument.");return 0;}
	if(!motifs->nmotifs){cmdError("No motifs.");return 0;}
	autodelete<RFMOCCA> r(new RFMOCCA(motifs->nmotifs));
	if(!r.ptr){
		outOfMemory();
		return 0;
	}
	r.ptr->mwin.ptr = motifWindow::create(motifs);
	if(!r.ptr->mwin.ptr){
		return 0;
	}
	r.ptr->moc = r.ptr->mwin.ptr->occContainer;
	r.ptr->motifs = r.ptr->mwin.ptr->motifs;
	r.ptr->classifier.ptr=logoddsClassifier::create(motifs->nmotifs);
	if(!r.ptr->classifier.ptr){
		return 0;
	}
	if(!r.ptr->fvec.resize((size_t)motifs->nmotifs)){
		return 0;
	}
	for(int l=0;l<motifs->nmotifs;l++){
		r.ptr->classifier.ptr->featureNames.push_back(std::string(motifs->motifs[l].name));
		RFMotifOccClassifier*sc = RFMotifOccClassifier::create(l,motifs);
		if(!sc){
			return 0;
		}
		r.ptr->subcls.push_back(sc);
	}
	return r.disown();
}

bool RFMOCCA::trainWindow(char*buf,long long pos,int bufs,seqClass*cls){
	trainingSequence*ts = trainingSequence::createSequenceClone(buf,bufs,cls);
	if(!ts)return false;
	trainSeq.push_back(ts);
	return true;
}

bool RFMOCCA::trainFinish(){
	for(trainingSequence*ts: trainSeq.v){
		mwin.ptr->flush();
		if(!mwin.ptr->readWindow(ts->seq,0,ts->length)){
			return false;
		}
		for(int l=0;l<motifs->nmotifs;l++){
			RFMotifOccClassifier*c=subcls[l];
			motifOcc*o=moc->getFirst(l);
			while(o){
				if(!c->trainOcc(o,moc,0,ts->seq,ts->length,ts->cls)){
					return false;
				}
				o=moc->getNextSame(o);
			}
		}
	}
	for(RFMotifOccClassifier*sc: subcls.v){
		if(!sc->trainFinish())return false;
	}
	for(trainingSequence*ts: trainSeq.v){
		mwin.ptr->flush();
		if(!mwin.ptr->readWindow(ts->seq,0,ts->length)){
			return false;
		}
		for(int x=0;x<motifs->nmotifs;x++){
			RFMotifOccClassifier*c=subcls[x];
			motifOcc*o=moc->getFirst(x);
			int nc=0;
			while(o){
				if( c->applyOcc(o,moc,0,ts->seq,ts->length)>0 ){
					nc++;
				}
				o=moc->getNextSame(o);
			}
			fvec[x]=double(nc*1000)/double(ts->length);
		}
		if(!classifier.ptr->addTrain(fvec.ptr,nFeatures,ts->cls)){
			return false;
		}
	}
	if(!classifier.ptr->train()){
		return false;
	}
	return true;
}

double RFMOCCA::do_applyWindow(char*buf,long long pos,int bufs){
	if(!mwin.ptr->readWindow(buf,pos,bufs)){
		return false;
	}
	for(int x=0;x<motifs->nmotifs;x++){
		RFMotifOccClassifier*c=subcls[x];
		motifOcc*o=moc->getFirst(x);
		int nc=0;
		while(o){
			if( c->applyOcc(o,moc,pos,buf,bufs)>0 ){
				nc++;
			}
			o=moc->getNextSame(o);
		}
		fvec[x]=double(nc*1000)/double(bufs);
	}
	return classifier.ptr->apply(fvec.ptr,nFeatures);
}

bool RFMOCCA::flush(){
	return mwin.ptr->flush();
}

bool RFMOCCA::printInfo(){
	cout << t_indent << "RF-MOCCA\n";
	cout << t_indent << "Trained?: " << (trained?"Yes":"No") << "\n";
	// Uncomment to always display model analysis
	//for(int l=0;l<motifs->nmotifs;l++){
	//	subcls[l]->printInfo();
	//}
	//classifier.ptr->printInfo((char*)"Log-odds");
	return true;
}

bool RFMOCCA::exportAnalysisData(string path){
	FILE*f=fopen(path.c_str(),"wb");
	fprintf(f, "RF-MOCCA model analysis\n");
	for(int l=0;l<motifs->nmotifs;l++){
		if(!subcls[l]->exportAnalysisData(f))
			return false;
	}
	if(!classifier.ptr->exportAnalysisData(f, (char*)"Log-odds", (char*)" - "))
		return false;
	fclose(f);
	cout << "Saved classifier analysis data to \"" << path << "\"\n";
	return true;
}

