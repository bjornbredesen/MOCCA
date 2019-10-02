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
#include "features.hpp"


////////////////////////////////////////////////////////////////////////////////////
// Feature info

/*
featureCardinality
	Enumeration of supported feature cardinalities.
*/
enum featureCardinality{
	featureCardinality_Invalid,
	featureCardinality_Singleton,
	featureCardinality_Linear,
	featureCardinality_Quadratic,
	featureCardinality_QuadraticReflexive,
};

/*
featureCardinalityL
	Indicates how many motif indices each cardinality used.
*/
int featureCardinalityL[]={
	0,	// Invalid
	0,	// Singleton
	1,	// Linear
	2,	// Quadratic
	2	// Quadratic reflexive
};

/*
	Feature specifications:
*/
typedef struct{
	char*name;
	featureCardinality cardinality;
}s_featureInfo;

s_featureInfo featureInfo[N_FEATURES]={
	{	(char*)"Invalid",		featureCardinality_Invalid,	},
	{	(char*)"nOcc", 			featureCardinality_Linear,	},
	{	(char*)"nPair", 		featureCardinality_QuadraticReflexive,	},
	{	(char*)"MDPA", 			featureCardinality_Linear,	},
	{	(char*)"MDP", 			featureCardinality_Quadratic,	},
	{	(char*)"MDM", 			featureCardinality_Quadratic,	},
	{	(char*)"MDDA", 			featureCardinality_Linear,	},
	{	(char*)"MDD", 			featureCardinality_Quadratic,	},
	{	(char*)"PEDI", 			featureCardinality_QuadraticReflexive,	},
	{	(char*)"nPairDH",		featureCardinality_QuadraticReflexive,	},
	{	(char*)"nPairCosI",		featureCardinality_QuadraticReflexive,	},
	{	(char*)"nOccPair",		featureCardinality_Quadratic,	},
	{	(char*)"GC content", 	featureCardinality_Singleton,	},
	{	(char*)"nPair2D",		featureCardinality_QuadraticReflexive,	},
	{	(char*)"nPairRgn",		featureCardinality_QuadraticReflexive,	},
};


////////////////////////////////////////////////////////////////////////////////////
// Feature set

featureSet::featureSet(){
	nfeatures=0;
	features.ptr=0;
	nifeatures=0;
	ifeatures.ptr=0;
	instantiated=false;
}

featureSet::~featureSet(){
	uninstantiate();
}

featureSet*featureSet::create(){
	featureSet*r=new featureSet();
	if(!r){
		outOfMemory();
		return 0;
	}
	return r;
}

featureSetFeature*featureSet::addFeature(featureType _f,int ia,int ib,int ic,double da,double db,int cf,int icf){
	nfeatures++;
	if(!features.resize(nfeatures)){
		return 0;
	}
	featureSetFeature*f=&features.ptr[nfeatures-1];
	memset(f,0,sizeof(featureSetFeature));
	f->f=_f;
	f->ia=ia;
	f->ib=ib;
	f->ic=ic;
	f->da=da;
	f->db=db;
	f->cf=cf;
	f->icf=icf;
	return f;
}

void featureSet::uninstantiate(){
	instantiated=false;
	nifeatures=0;
}

featureSetInstFeature*featureSet::addInstFeature(featureSetFeature*bf,int ia,int ib,int ic,int icf){
	if(!bf)return 0;
	nifeatures++;
	if(!ifeatures.resize(nifeatures)){
		return 0;
	}
	featureSetInstFeature*f=&ifeatures.ptr[nifeatures-1];
	memset(f,0,sizeof(featureSetInstFeature));
	f->fsf=bf;
	f->ia=ia;
	f->ib=ib;
	f->ic=ic;
	f->icf=icf;
	return f;
}

bool featureSet::instantiateFeature(featureSetFeature*f,motifList*ml){
	if(!f)return false;
	/*
		Instantiate features depending on cardinality and motifs
	*/
	int icfa=0,icfb=f->cf-1;
	if(f->icf!=featureMotif_All)icfa=icfb=f->icf;
	for(int icf=icfa;icf<=icfb;icf++)switch(featureInfo[(int)f->f].cardinality){
		case featureCardinality_Singleton:{
			if(!addInstFeature(f,featureMotif_All,featureMotif_All,featureMotif_All,icf))
				return false;
			break;}
		case featureCardinality_Linear:{
			if(f->ia==featureMotif_All){
				for(int l=0;l<ml->nmotifs;l++){
					if(!addInstFeature(f,l,featureMotif_All,featureMotif_All,icf))
						return false;
				}
			}else{
				if(!addInstFeature(f,f->ia,featureMotif_All,featureMotif_All,icf))
					return false;
			}
			break;}
		case featureCardinality_Quadratic:{
			if(f->ia==featureMotif_All&&f->ib==featureMotif_All){
				for(int x=0;x<ml->nmotifs;x++){
					for(int y=0;y<ml->nmotifs;y++){
						if(!addInstFeature(f,x,y,featureMotif_All,icf))
							return false;
					}
				}
			}else if(f->ia!=featureMotif_All&&f->ib!=featureMotif_All){
				if(!addInstFeature(f,f->ia,f->ib,featureMotif_All,icf))
					return false;
			}else{
				cmdError("Invalid feature.");
				return false;
			}
			break;}
		case featureCardinality_QuadraticReflexive:{
			if(f->ia==featureMotif_All&&f->ib==featureMotif_All){
				for(int x=0;x<ml->nmotifs;x++){
					for(int y=0;y<=x;y++){
						if(!addInstFeature(f,x,y,featureMotif_All,icf))
							return false;
					}
				}
			}else if(f->ia!=featureMotif_All&&f->ib!=featureMotif_All){
				if(!addInstFeature(f,f->ia,f->ib,featureMotif_All,icf))
					return false;
			}else{
				cmdError("Invalid feature.");
				return false;
			}
			break;}
		default:;
	}
	return true;
}

bool featureSet::instantiate(motifList*ml){
	if(instantiated){
		cmdError("Feature set already instantiated.");
		return false;
	}
	if(!ml||!ml->nmotifs)return false;
	featureSetFeature*f=features;
	for(int l=0;l<nfeatures;l++,f++){
		if(!instantiateFeature(f,ml))
			return false;
	}
	instantiated=true;
	return true;
}

bool featureSet::skipUnusedMotifs(motifList*ml){
	if(!instantiated){
		cmdError("Feature set not instantiated.");
		return false;
	}
	motifListMotif*m=ml->motifs;
	for(int l=0;l<ml->nmotifs;l++,m++)
		m->skip=true;
	featureSetInstFeature*f=ifeatures;
	for(int l=0;l<nifeatures;l++,f++){
		int cardl=featureCardinalityL[(int)featureInfo[(int)f->fsf->f].cardinality];
		if(cardl>=1)
			ml->motifs[f->ia].skip=false;
		if(cardl>=2)
			ml->motifs[f->ib].skip=false;
		if(cardl>=3)
			ml->motifs[f->ic].skip=false;
	}
	int nskip=0;
	m=ml->motifs;
	for(int l=0;l<ml->nmotifs;l++,m++)
		if(m->skip)nskip++;
	cout << t_indent << "Skipped motifs: " << nskip << "/" << ml->nmotifs << "\n";
	return true;
}

void featureSet::printInfo(){
	cmdSection("Features");
	featureSetFeature*f=features;
	for(int l=0;l<nfeatures;l++,f++){
		cout << t_indent << featureInfo[(int)f->f].name;
		switch(f->f){
			case featureType_Invalid:break;
			case featureType_nOcc:break;
			case featureType_nPair:cout << "(" << f->da << ")";break;
			case featureType_nPairDH:cout << "(" << f->da << ")";break;
			case featureType_PEDI:cout << "(" << f->da << ", " << f->db << ")";break;
			case featureType_MDPA:break;
			case featureType_MDP:break;
			case featureType_MDDA:break;
			case featureType_MDD:break;
			case featureType_MDM:break;
			case featureType_nOccPair:cout << "(" << f->da << ")";break;
			case featureType_GC:break;
			case featureType_nPair2D:cout << "(" << f->da << ", " << f->db << ", " << (f->icf?"Y-axis":"X-axis") << ")";break;
			case featureType_nPairRgn:cout << "(" << f->da << ")";break;
		}
		cout << "\n";
	}
}

void featureSet::printInfoI(motifList*ml){
	if(!nifeatures){
		return;
	}
	cmdSection("Instantiated features");
	cout << t_indent << "# instantiated features: " << nifeatures << "\n";
	featureSetInstFeature*f=ifeatures;
	for(int l=0;l<nifeatures;l++,f++){
		cout << t_indent;
		printInstFeatureName(f,ml);
		cout << "\n";
	}
}

void featureSet::printInstFeatureName(featureSetInstFeature*f,motifList*ml){
	int cardl=featureCardinalityL[(int)featureInfo[(int)f->fsf->f].cardinality];
	cout << featureInfo[(int)f->fsf->f].name;
	switch(f->fsf->f){
		case featureType_Invalid:break;
		case featureType_nOcc:break;
		case featureType_nPair:cout << "(" << f->fsf->da << ")";break;
		case featureType_nPairDH:cout << "(" << f->fsf->da << ")";break;
		case featureType_PEDI:cout << "(" << f->fsf->da << ", " << f->fsf->db << ")";break;
		case featureType_MDPA:break;
		case featureType_MDP:break;
		case featureType_MDDA:break;
		case featureType_MDD:break;
		case featureType_MDM:break;

		case featureType_nOccPair:cout << "(" << f->fsf->da << ")";break;
		case featureType_GC:break;
		case featureType_nPair2D:{
			cout << "(" << f->fsf->da << ", " << f->fsf->db << ")" << t_indent << (f->icf?"Y-axis":"X-axis");
			break;}
		case featureType_nPairRgn:{
			bool strand=f->icf&1;
			int rgn=f->icf>>1;
			char*rgnn[3]={ (char*)"Upstream", (char*)"Proximal", (char*)"Downstream" };
			char*strn[2]={ (char*)"Same strand", (char*)"Opposite strand" };
			cout << "(" << f->fsf->da << ")" << t_indent << rgnn[rgn] << t_indent << strn[strand];
			break;}
	}
	if(cardl>=1)
		cout << t_indent << ml->motifs[f->ia].name;
	if(cardl>=2)
		cout << ", " << ml->motifs[f->ib].name;
}

std::vector<std::string> featureSet::getInstFeatureNames(motifList*ml){
	std::vector<std::string> names;
	if(!nifeatures){
		return names;
	}
	featureSetInstFeature*f=ifeatures;
	for(int l=0;l<nifeatures;l++,f++){
		/*cout << t_indent;
		printInstFeatureName(f,ml);
		cout << "\n";*/
		std::stringstream cname;
		int cardl=featureCardinalityL[(int)featureInfo[(int)f->fsf->f].cardinality];
		cname << featureInfo[(int)f->fsf->f].name;
		switch(f->fsf->f){
			case featureType_Invalid:break;
			case featureType_nOcc:break;
			case featureType_nPair:cname << "(" << f->fsf->da << ")";break;
			case featureType_nPairDH:cname << "(" << f->fsf->da << ")";break;
			case featureType_PEDI:cname << "(" << f->fsf->da << ", " << f->fsf->db << ")";break;
			case featureType_MDPA:break;
			case featureType_MDP:break;
			case featureType_MDDA:break;
			case featureType_MDD:break;
			case featureType_MDM:break;
			case featureType_nOccPair:cname << "(" << f->fsf->da << ")";break;
			case featureType_GC:break;
			case featureType_nPair2D:{
				cname << "(" << f->fsf->da << ", " << f->fsf->db << ")" << t_indent << (f->icf?"Y-axis":"X-axis");
				break;}
			case featureType_nPairRgn:{
				bool strand=f->icf&1;
				int rgn=f->icf>>1;
				char*rgnn[3]={ (char*)"Upstream", (char*)"Proximal", (char*)"Downstream" };
				char*strn[2]={ (char*)"Same strand", (char*)"Opposite strand" };
				cname << "(" << f->fsf->da << ")" << t_indent << rgnn[rgn] << t_indent << strn[strand];
				break;}
		}
		if(cardl>=1)
			cname << t_indent << ml->motifs[f->ia].name;
		if(cardl>=2)
			cname << ", " << ml->motifs[f->ib].name;
		names.push_back(cname.str());
	}
	return names;
}


////////////////////////////////////////////////////////////////////////////////////
// Feature window

featureWindow::featureWindow(motifWindow*mw,featureSet*fs){
	cfg=getConfiguration();
	fvec=0;
	features=fs;
	mwin=mw;
	motifs=mw->motifs;
	occContainer=mw->occContainer;
}

featureWindow::~featureWindow(){
	if(fvec)free(fvec);
}

featureWindow*featureWindow::create(motifWindow*mw,featureSet*fs){
	if(!mw||!fs){
		cmdError("featureWindow::create(): Null-pointer arguments.");
		return 0;
	}
	featureWindow*r=new featureWindow(mw,fs);
	if(!r){
		outOfMemory();
		return 0;
	}
	if(!fs->instantiated){
		// If the feature set has not already been instantiated, do it now.
		if(!fs->instantiate(mw->motifs)){
			delete r;
			return 0;
		}
	}
	if(!fs->skipUnusedMotifs(mw->motifs)){
		delete r;
		return 0;
	}
	int nf=fs->nifeatures;
	if(!nf){
		cmdError("Empty feature set.");
		delete r;
		return 0;
	}
	r->fvec=(double*)malloc(sizeof(double)*nf);
	if(!r->fvec){
		outOfMemory();
		delete r;
		return 0;
	}
	return r;
}

double featureWindow::getMDP(int a,int b){
	if(!occContainer->nOccT[a]||!occContainer->nOccT[b]||(a==b&&occContainer->nOccT[a]<2)){
		return double(cfg->windowSize);
	}
	double mdp=0;
	long nmdp=0;
	motifOcc*o=occContainer->getFirst(a);
	while(o){
		if(o->skip){o=occContainer->getNextSame(o);continue;}
		motifOcc*o2=occContainer->getFirst(b);
		double dp=cfg->windowSize;
		while(o2){
			if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
			if(o!=o2){
				double d=(double(o->start)+double(o->mot->len)*0.5)-(double(o2->start)+double(o2->mot->len)*0.5);
				if(d<0)d=-d;
				if(d<dp)dp=d;
			}
			o2=occContainer->getNextSame(o2);
		}
		if(dp<cfg->windowSize){
			mdp+=dp;
			nmdp++;
		}
		o=occContainer->getNextSame(o);
	}
	return nmdp?mdp/double(nmdp):double(cfg->windowSize);
}

double featureWindow::getMDPA(int a){
	if(!occContainer->nOccT[a]||occContainer->nOcc<2){
		return double(cfg->windowSize);
	}
	double mdp=0;
	long nmdp=0;
	motifOcc*o=occContainer->getFirst(a);
	while(o){
		if(o->skip){o=occContainer->getNextSame(o);continue;}
		motifOcc*o2=occContainer->getFirst();
		double dp=cfg->windowSize;
		while(o2){
			if(o2->skip){o2=occContainer->getNext(o2);continue;}
			if(o!=o2){
				double d=(double(o->start)+double(o->mot->len)*0.5)-(double(o2->start)+double(o2->mot->len)*0.5);
				if(d<0)d=-d;
				if(d<dp)dp=d;
			}
			o2=occContainer->getNext(o2);
		}
		if(dp<cfg->windowSize){
			mdp+=dp;
			nmdp++;
		}
		o=occContainer->getNextSame(o);
	}
	return nmdp?mdp/double(nmdp):double(cfg->windowSize);
}

double featureWindow::getMDM(int a,int b){
	if(!occContainer->nOccT[a]||!occContainer->nOccT[b]||(a==b&&occContainer->nOccT[a]<2)){
		return double(cfg->windowSize);
	}
	double mdm=0;
	long nmdm=0;
	motifOcc*o=occContainer->getFirst(a);
	while(o){
		if(o->skip){o=occContainer->getNextSame(o);continue;}
		motifOcc*o2=occContainer->getFirst(b);
		int dm=0;
		int ndm=0;
		while(o2){
			if(o!=o2){
				if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
				int d1=int(o2->start-(o->start+o->mot->len));
				int d2=int(o->start-(o2->start+o2->mot->len));
				int d_alpha=max(max(d1,d2),0);
				if(d_alpha<=cfg->windowSize){
					dm+=d_alpha;
					ndm++;
				}
			}
			o2=occContainer->getNextSame(o2);
		}
		if(ndm){
			mdm+=double(dm)/double(ndm);
			nmdm++;
		}
		o=occContainer->getNextSame(o);
	}
	return nmdm?double(mdm)/double(nmdm):double(cfg->windowSize);
}

double featureWindow::getMDDA(int a){
	if(!occContainer->nOccT[a]||occContainer->nOcc<2){
		return double(cfg->windowSize);
	}
	long mdd=0;
	long nmdd=0;
	motifOcc*o=occContainer->getFirst(a);
	while(o){
		if(o->skip){o=occContainer->getNextSame(o);continue;}
		motifOcc*o2=occContainer->getFirst();
		int dd=-1;
		while(o2){
			if(o2->skip){o2=occContainer->getNext(o2);continue;}
			if(o!=o2){
				int d1=int(o2->start-(o->start+o->mot->len));
				int d2=int(o->start-(o2->start+o2->mot->len));
				int d_alpha=max(max(d1,d2),0);
				if(d_alpha<=cfg->windowSize){
					if(dd==-1){
						dd=d_alpha;
					}else{
						dd=max(dd,d_alpha);
					}
				}
			}
			o2=occContainer->getNext(o2);
		}
		if(dd!=-1){
			mdd+=dd;
			nmdd++;
		}
		o=occContainer->getNextSame(o);
	}
	return nmdd?double(mdd)/double(nmdd):double(cfg->windowSize);
}

double featureWindow::getMDD(int a,int b){
	if(!occContainer->nOccT[a]||!occContainer->nOccT[b]||(a==b&&occContainer->nOccT[a]<2)){
		return double(cfg->windowSize);
	}
	long mdd=0;
	long nmdd=0;
	motifOcc*o=occContainer->getFirst(a);
	while(o){
		if(o->skip){o=occContainer->getNextSame(o);continue;}
		motifOcc*o2=occContainer->getFirst(b);
		int dd=-1;
		while(o2){
			if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
			if(o!=o2){
				int d1=int(o2->start-(o->start+o->mot->len));
				int d2=int(o->start-(o2->start+o2->mot->len));
				int d_alpha=max(max(d1,d2),0);
				if(d_alpha<=cfg->windowSize){
					if(dd==-1){
						dd=d_alpha;
					}else{
						dd=max(dd,d_alpha);
					}
				}
			}
			o2=occContainer->getNextSame(o2);
		}
		if(dd!=-1){
			mdd+=dd;
			nmdd++;
		}
		o=occContainer->getNextSame(o);
	}
	return nmdd?double(mdd)/double(nmdd):double(cfg->windowSize);
}

int featureWindow::getNFeatures(){
	return features->nifeatures;
}

double*featureWindow::extractFeatures(char*wseq,long long wpos,int wlen,bool doRead){
	if(doRead){
		if(!mwin->readWindow(wseq,wpos,wlen)){
			return 0;
		}
	}
	featureSetFeature*fsf;
	featureSetInstFeature*fsif=features->ifeatures;
	double v;
	int iv;
	double normv=1.0/double(wlen);
	for(int l=0;l<features->nifeatures;l++,fsif++){
		fsf=fsif->fsf;
		switch(fsf->f){
			case featureType_nOcc:{
				int n=0;
				motifOcc*o=occContainer->getFirst(fsif->ia);
				while(o){
					if(!o->skip)n++;
					o=occContainer->getNextSame(o);
				}
				v=double(n)*normv;
				break;}
			case featureType_nOccPair:{
				iv=0;
				int nPairCut=int(fsf->da);
				motifOcc*o=occContainer->getFirst(fsif->ia);
				while(o){
					if(o->skip){o=occContainer->getNextSame(o);continue;}
					motifOcc*o2=occContainer->getFirst(fsif->ib);
					while(o2){
						if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
						if(o!=o2){
							double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
												-(double(o->start)+double(o->mot->len)/2.0);
							if(d_gamma<0)d_gamma=-d_gamma;
							if(d_gamma<=nPairCut){
								iv++;
								// This only needs to know if the first motif occurrence is paired
								// with one of the other type and then count it, so break when paired.
								break;
							}
						}
						o2=occContainer->getNextSame(o2);
					}
					o=occContainer->getNextSame(o);
				}
				v=double(iv)*normv;
				break;}
			case featureType_nPair:{
				iv=0;
				int nPairCut=int(fsf->da);
				// New
				if(fsif->ia==fsif->ib){
					motifOcc*o=occContainer->getFirst(fsif->ia),*o2,*no;
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						no=o2=occContainer->getNextSame(o);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							int d1=int(o2->start-(o->start+o->mot->len));
							int d2=int(o->start-(o2->start+o2->mot->len));
							int d_alpha=max(max(d1,d2),0);
							if(d_alpha<=nPairCut){
								iv++;
							}
							o2=occContainer->getNextSame(o2);
						}
						o=no;
					}
				}else{
					motifOcc*o=occContainer->getFirst(fsif->ia);
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						motifOcc*o2=occContainer->getFirst(fsif->ib);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							int d1=int(o2->start-(o->start+o->mot->len));
							int d2=int(o->start-(o2->start+o2->mot->len));
							int d_alpha=max(max(d1,d2),0);
							if(d_alpha<=nPairCut){
								iv++;
							}
							o2=occContainer->getNextSame(o2);
						}
						o=occContainer->getNextSame(o);
					}
				}
				v=double(iv)*normv;
				break;}
			case featureType_nPair2D:{
				double s=0;
				bool axis=fsif->icf&2;
				double freq=fsf->db;
				int nPairCut=int(fsf->da);
				if(fsif->ia==fsif->ib){
					motifOcc*o=occContainer->getFirst(fsif->ia),*o2,*no;
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						no=o2=occContainer->getNextSame(o);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
												-(double(o->start)+double(o->mot->len)/2.0);
							if(d_gamma<0)d_gamma=-d_gamma;
							if(d_gamma<=nPairCut){
								double ph=d_gamma/freq;
								if(ph<0)ph=-ph;
								s+=axis?sin(ph*3.141592654*2.0):cos(ph*3.141592654*2.0);
							}
							o2=occContainer->getNextSame(o2);
						}
						o=no;
					}
				}else{
					motifOcc*o=occContainer->getFirst(fsif->ia);
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						motifOcc*o2=occContainer->getFirst(fsif->ib);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
												-(double(o->start)+double(o->mot->len)/2.0);
							if(d_gamma<=nPairCut&&d_gamma>=-nPairCut){
								double ph=d_gamma/freq;
								if(o->strand)ph=-ph;
								s+=axis?sin(ph*3.141592654*2.0):cos(ph*3.141592654*2.0);
							}
							o2=occContainer->getNextSame(o2);
						}
						o=occContainer->getNextSame(o);
					}
				}
				v=s*normv;
				break;}
			case featureType_nPairRgn:{
				// Cardinality x 6 (regions and strands)
				double s=0;
				
				int nPairCut=int(fsf->da);
				bool strand=fsif->icf&1;
				int rgn=fsif->icf>>1;
				
				double dMin=0,dMax=0;
				
				double left=-double(nPairCut)/2.0;
				double dd=double(nPairCut)/3.0;
				
				switch(rgn){
				case 0: dMin=left; dMax=left+dd; break;
				case 1: dMin=left+dd; dMax=left+dd*2.0; break;
				case 2: dMin=left+dd*2.0; dMax=left+dd*3.0; break;
				default:cmdWarning("Invalid region");
				}

				
				motifOcc*o=occContainer->getFirst(fsif->ia);
				while(o){
					if(o->skip){o=occContainer->getNextSame(o);continue;}
					motifOcc*o2=occContainer->getFirst(fsif->ib);
					while(o2){
						if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
						
						// Only consider pairs with the right strand combination.
						if( (!strand&&o->strand!=o2->strand)
						||	(strand&&o->strand==o2->strand) ){
							o2=occContainer->getNextSame(o2);
							continue;
						}

						double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
											-(double(o->start)+double(o->mot->len)/2.0);

						// If main occurrence is on the other strand, reverse the distance
						if(o->strand)d_gamma=-d_gamma;
						
						if(d_gamma>=dMin&&d_gamma<=dMax){
							s+=1.0;
						}

						o2=occContainer->getNextSame(o2);
					}
					o=occContainer->getNextSame(o);
				}

				v=s*normv;
				break;}
			case featureType_nPairDH:{
				v=0;
				int nPairCut=int(fsf->da);
				if(fsif->ia==fsif->ib){
					motifOcc*o=occContainer->getFirst(fsif->ia),*o2,*no;
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						no=o2=occContainer->getNextSame(o);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
												-(double(o->start)+double(o->mot->len)/2.0);
							if(d_gamma<=nPairCut&&d_gamma>=-nPairCut){
								double phaseshift=0;
								if(o->strand!=o2->strand)phaseshift=5.25;
								v+=cos(((double(d_gamma)+phaseshift)/10.5)*3.141592654*2.0)+1.0;
							}
							o2=occContainer->getNextSame(o2);
						}
						o=no;
					}
				}else{
					motifOcc*o=occContainer->getFirst(fsif->ia);
					while(o){
						if(o->skip){o=occContainer->getNextSame(o);continue;}
						motifOcc*o2=occContainer->getFirst(fsif->ib);
						while(o2){
							if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
							double d_gamma=(double(o2->start)+double(o2->mot->len)/2.0)
												-(double(o->start)+double(o->mot->len)/2.0);
							if(d_gamma<=nPairCut&&d_gamma>=-nPairCut){
								double phaseshift=0;
								if(o->strand!=o2->strand)phaseshift=5.25;
								v+=cos(((double(d_gamma)+phaseshift)/10.5)*3.141592654*2.0)+1.0;
							}
							o2=occContainer->getNextSame(o2);
						}
						o=occContainer->getNextSame(o);
					}
				}
				v*=normv/2.0;
				break;}
			case featureType_PEDI:{
				v=0;
				int nPairCut=int(fsf->da);
				double delta=fsf->db;
				motifOcc*o=occContainer->getFirst(fsif->ia);
				while(o){
					if(o->skip){o=occContainer->getNextSame(o);continue;}
					motifOcc*o2=occContainer->getFirst(fsif->ib);
					while(o2){
						if(o2->skip){o2=occContainer->getNextSame(o2);continue;}
						int d1=int(o2->start-(o->start+o->mot->len));
						int d2=int(o->start-(o2->start+o2->mot->len));
						int d_alpha=max(max(d1,d2),0);
						if(d_alpha<=nPairCut){
							double C=cos((double(d_alpha)/delta)*3.141592654*2.0);
							v+=(C+1.0)/2.0;
						}
						o2=occContainer->getNextSame(o2);
					}
					o=occContainer->getNextSame(o);
				}
				v*=normv;
				break;}
			case featureType_MDPA:{
				v=getMDPA(fsif->ia);
				break;}
			case featureType_MDP:{
				v=getMDP(fsif->ia,fsif->ib);
				break;}
			case featureType_MDDA:{
				v=getMDDA(fsif->ia);
				break;}
			case featureType_MDD:{
				v=getMDD(fsif->ia,fsif->ib);
				break;}
			case featureType_MDM:{
				v=getMDM(fsif->ia,fsif->ib);
				break;}
			case featureType_GC:{
				char*c=wseq;
				int gc=0;
				for(int i=0;i<wlen;i++,c++){
					if(*c=='G'||*c=='C')gc++;
				}
				v=double(gc)/double(wlen);
				break;}
			default:cmdError("Unsupported feature."); return 0;
		}
		fvec[l]=v;
	}
	return fvec;
}

