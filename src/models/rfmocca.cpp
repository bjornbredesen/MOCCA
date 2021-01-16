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

vector<prediction> RFMOCCA::predictWindow(char*buf,long long pos,int bufs, corePredictionModeT cpm){
	vector<prediction> ret = vector<prediction>();
	if(!mwin.ptr->readWindow(buf,pos,bufs)){
		return ret;
	}
	vector<prediction> motpos = vector<prediction>();
	for(int x=0;x<motifs->nmotifs;x++){
		RFMotifOccClassifier*c=subcls[x];
		motifOcc*o=moc->getFirst(x);
		int nc=0;
		while(o){
			if( c->applyOcc(o,moc,pos,buf,bufs)>0 ){
				nc++;
				if(cpm != cpmNone){
					int center = (o->start + o->start + motifs->motifs[x].len) / 2;
					double mscore = classifier.ptr->getWeight(x);
					motpos.push_back(prediction(center - 250, center + 250, mscore,
						o->start,
						o->start + motifs->motifs[x].len));
				}
			}
			o=moc->getNextSame(o);
		}
		fvec[x]=double(nc*1000)/double(bufs);
	}
	double cvalue = classifier.ptr->apply(fvec.ptr,nFeatures);
	if(cpm != cpmNone){
		if(cvalue >= threshold){
			// Motif prediction--based core prediction algorithm
			if(cpm == cpmMotifs || cpm == cpmMotifsStrong){
				// Try to limit to central occurrences with cumulative scores above threshold
				vector<prediction> cmotpos = vector<prediction>();
				sort(motpos.begin(),motpos.end(),
				[](const prediction a,const prediction b){
					return a.mstart < b.mstart;
				});
				for(auto&oA: motpos){
					double cumscore = 0.;
					for(auto&oB: motpos){
						if(oB.mstart < oA.start) continue;
						if(oB.mend > oA.end) break;
						cumscore += oB.score;
					}
					if(cumscore >= threshold)
						cmotpos.push_back(prediction(
							max(oA.start, 0),
							min(oA.end, (int)pos + bufs),
							cumscore));
				}
				if(cpm == cpmMotifsStrong) return cmotpos;
				if(cmotpos.size() > 0) motpos = cmotpos;
				// Flatten
				for(auto&occ: motpos){
					//if(ret.size() > 0 && occ.start < ret.back().end+500){
					if(ret.size() > 0 && occ.start < ret.back().end){
						ret.back().end = max(ret.back().end, occ.end);
						ret.back().score = max(ret.back().score, occ.score);
					}else{
						ret.push_back(prediction(occ.start, occ.end, occ.score));
					}
				}
				// If there are predictions scoring above the threshold, predict those.
				// Otherwise, there must be many low-scoring predictions, so predict
				// all of them.
				// Thus: Get predictions with scores above threshold.
				vector<prediction> pred = vector<prediction>();
				for(auto&p: ret)
					if(p.score >= threshold)
						pred.push_back(p);
				// If any, return that list
				if(pred.size() > 0)
					return pred;
				// If none, return the regular list of predictions
			// Continuous core prediction algorithm
			}else if(cpm == cpmContinuous){
				// For continuous predictions, we want to predict a core delimited
				// by motif occurrences that spans potentially a larger region, and
				// has a high score per basepair.
				sort(motpos.begin(),motpos.end(),
				[](const prediction a,const prediction b){
					return a.center < b.center;
				});
				vector<prediction> cmotpos = vector<prediction>();
				for(auto&m: motpos) cmotpos.push_back(m);
				sort(cmotpos.begin(),cmotpos.end(),
				[](const prediction a,const prediction b){
					return a.mstart < b.mstart;
				});
				int maxWinStart = -1;
				int maxWinEnd = -1;
				double maxWinScore = -1.;
				double maxWinCumScore = -1.;
				int iCA = 0;
				for(int iA = 0; iA < motpos.size() - 1; iA++){
					auto&oA = motpos[iA];
					int iCB = iCA;
					double cumscoreA = 0.;
					for(int iB = iA; iB < motpos.size(); iB++){
						auto&oB = motpos[iB];
						if(oB.end > oA.start + cfg->windowSize)
							break;
						int wA = max(oA.start, 0);
						int wB = min(oB.end, (int)pos + bufs);
						if(wB-wA < oA.end-oA.start) continue;
						bool block = false;
						double cumscore = cumscoreA;
						for(int iC = iCB; iC < cmotpos.size(); iC++){
							auto&oC = cmotpos[iC];
							// Occurrences are sorted by motif start, so if they
							// start before the window, we can safely skip.
							if(oC.mstart < wA){
								iCA++;
								iCB++;
								continue;
							}
							//if(oC.mstart > wB) continue;
							// If the occurrence is past the window end, we can
							// safely break.
							if(oC.mstart > wB) break;
							// If only the end is past, later shorter occurrences
							// may still be inside, so just block updating of the
							// bookmark position.
							if(oC.mend > wB) {
								block = true;
								continue;
							}
							// The occurrence is inside the window, so add.
							cumscore += oC.score;
							// If an occurrence was past the end of the window,
							// do not update bookmarks.
							if(!block) {
								iCB = iC + 1;
								cumscoreA = cumscore;
							}
						}
						double winScore = cumscore / max(double(wB - wA), 1.);
						if(maxWinStart == -1 || winScore > maxWinScore){
							maxWinStart = wA;
							maxWinEnd = wB;
							maxWinScore = winScore;
							maxWinCumScore = cumscore;
						}
					}
				}
				if(maxWinStart == -1) return ret;
				ret.push_back(prediction(maxWinStart, maxWinEnd, maxWinCumScore));
				/*
				// Unoptimized base algorithm
				sort(motpos.begin(),motpos.end(),
				[](const prediction a,const prediction b){
					return a.center < b.center;
				});
				vector<prediction> mwnd = vector<prediction>();
				for(int iA = 0; iA < motpos.size() - 1; iA++){
					auto&oA = motpos[iA];
					for(int iB = iA; iB < motpos.size(); iB++){
						auto&oB = motpos[iB];
						if(oB.end > oA.start + cfg->windowSize)
							break;
						int wA = max(oA.start, 0);
						int wB = min(oB.end, (int)pos + bufs);
						if(wB-wA < oA.end-oA.start) continue;
						double cumscore = 0.;
						for(auto&oC: motpos){
							if(oC.mstart < wA || oC.mend > wB) continue;
							cumscore += oC.score;
						}
						mwnd.push_back(prediction(wA, wB, cumscore));
					}
				}
				if(mwnd.size() == 0) return mwnd;
				sort(mwnd.begin(),mwnd.end(),
				[](const prediction a,const prediction b){
					double sA = a.score / max(double(a.end - a.start), 1.);
					double sB = b.score / max(double(b.end - b.start), 1.);
					return sA > sB;
				});
				ret.push_back(prediction(mwnd[0].start, mwnd[0].end, mwnd[0].score));
				*/
			}
		}
	}else if(cvalue >= threshold){
		ret.push_back(prediction(pos, pos+bufs, cvalue));
	}
	return ret;
}

