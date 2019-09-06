////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "aux.hpp"
#include "validation.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Validation

double getROCAUCxVP(validationPair*vp,int nvp,double x){
	if(!vp||!nvp){
		cmdError("Invalid validation set.");
		return false;
	}
	int TP = 0, FP = 0, nP = 0, nN = 0;
	vector<validationPair> svp;
	for(int x=0;x<nvp;x++){
		svp.push_back(vp[x]);
		if(vp[x].cls->flag) nP++;
		else nN++;
	}
	sort(svp.begin(),svp.end(),
	[](const validationPair a,const validationPair b){
		return a.score > b.score;
	});
	double aTPR = 0, bTPR = 0;
	double aFPR = 0, bFPR = 0;
	double AUC = 0.0;
	for(auto& cvp: svp){
		if(cvp.cls->flag){
			TP++;
		}else{
			FP++;
		}
		aFPR = bFPR;
		bFPR = double(FP)/double(nN);
		aTPR = bTPR;
		bTPR = double(TP)/double(nP);
		if(bFPR>=x){
			bTPR = aTPR + (bTPR-aTPR)*(x-aFPR)/(bFPR-aFPR);
			AUC += ( x-aFPR ) * ( bTPR+aTPR ) * 0.5;
			break;
		}
		AUC += ( bFPR-aFPR ) * ( bTPR+aTPR ) * 0.5;
	}
	return AUC;
}

double getPRCAUCxVP(validationPair*vp,int nvp,double x){
	if(!vp||!nvp){
		cmdError("Invalid validation set.");
		return false;
	}
	int TP = 0, FP = 0, nP = 0, nN = 0;
	vector<validationPair> svp;
	for(int x=0;x<nvp;x++){
		svp.push_back(vp[x]);
		if(vp[x].cls->flag) nP++;
		else nN++;
	}
	sort(svp.begin(),svp.end(),
	[](const validationPair a,const validationPair b){
		return a.score > b.score;
	});
	double aPrecision = 0, bPrecision = 0;
	double aRecall = 0, bRecall = 0;
	double AUC = 0.0;
	for(auto& cvp: svp){
		if(cvp.cls->flag){
			TP++;
		}else{
			FP++;
		}
		aRecall = bRecall;
		bRecall = double(TP)/double(nP);
		aPrecision = bPrecision;
		bPrecision = double(TP)/double(TP+FP);
		if(bRecall>=x){
			bPrecision = aPrecision + (bPrecision-aPrecision)*(x-aRecall)/(bRecall-aRecall);
			AUC += ( x-aRecall ) * ( bPrecision+aPrecision ) * 0.5;
			break;
		}
		AUC += ( bRecall-aRecall ) * ( bPrecision+aPrecision ) * 0.5;
	}
	return AUC;
}

bool printValidationMeasures(validationPair*vp,int nvp,double threshold){
	if(!vp||!nvp){
		cmdError("Invalid validation set.");
		return false;
	}
	cout << t_indent << t_indent << "Validation measures:\n";
	int TP=0,FP=0,TN=0,FN=0;
	for(int l=0;l<nvp;l++){
		bool pcls=vp[l].score>=threshold;
		bool tcls=vp[l].cls->flag;
		if(pcls&&tcls)TP++;
		else if(pcls&&!tcls)FP++;
		else if(!pcls&&tcls)FN++;
		else TN++;
	}
	double ACC=double(TP+TN)/double(TP+TN+FP+FN);
	double PPV=double(TP)/double(TP+FP);
	cout << t_indent << t_indent << t_indent << "TP = " << TP << t_indent << "TN = " << TN << "\n";
	cout << t_indent << t_indent << t_indent << "FP = " << FP << t_indent << "FN = " << FN << "\n";
	cout << t_indent << t_indent << t_indent << "Accuracy = " << (ACC*100.0) << "%\n";
	cout << t_indent << t_indent << t_indent << "PPV = " << (PPV*100.0) << "%\n";
	cout << t_indent << t_indent << t_indent << "ROC AUC(0.1) = " << getROCAUCxVP(vp,nvp,0.1) << "\n";
	cout << t_indent << t_indent << t_indent << "ROC AUC(0.5) = " << getROCAUCxVP(vp,nvp,0.5) << "\n";
	cout << t_indent << t_indent << t_indent << "ROC AUC(1.0) = " << getROCAUCxVP(vp,nvp,1.0) << "\n";
	cout << t_indent << t_indent << t_indent << "PRC AUC(0.1) = " << getPRCAUCxVP(vp,nvp,0.1) << "\n";
	cout << t_indent << t_indent << t_indent << "PRC AUC(0.5) = " << getPRCAUCxVP(vp,nvp,0.5) << "\n";
	cout << t_indent << t_indent << t_indent << "PRC AUC(1.0) = " << getPRCAUCxVP(vp,nvp,1.0) << "\n";
	return true;
}

bool saveVPairTable(char*outpath,validationPair*vp,int nvp){
	FILE*fout=fopen(outpath,"wb");
	if(!fout){
		cout << m_error << "Could not open file \"" << outpath << "\" for writing.\n";
		return false;
	}
	fprintf(fout,"score\tclass\n");
	for(int l=0;l<nvp;l++){
		fprintf(fout,"%.10f\t%s\n",vp[l].score,vp[l].cls->flag?(char*)"+":(char*)"-");
	}
	fclose(fout);
	return true;
}

