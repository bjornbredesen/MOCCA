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
#include "../sequencelist.hpp"
#include "baseclassifier.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Base classifiers
//	These classifiers only classify vectors of values.

baseClassifierSmp::baseClassifierSmp(double*v,seqClass*c,double e){
	cls = c;
	vecl = 0;
	clsE = e;
}

baseClassifierSmp*baseClassifierSmp::create(double*v,int vl,seqClass*c){
	if(!v||!c){
		cmdError("Invalid training example.");
		return 0;
	}
	double e = c->flag?1.0:-1.0;
	baseClassifierSmp*r = new baseClassifierSmp(v,c,e);
	if(!r){
		outOfMemory();
		return 0;
	}
	r->vec.resize((size_t)vl);
	if(!r->vec.ptr){
		delete r;
		return 0;
	}
	memcpy(r->vec.ptr,v,sizeof(double)*vl);
	r->vecl = vl;
	return r;
}

baseClassifier::baseClassifier(int nf){
	nFeatures=nf;
	threshold=0.0;
	trained=false;
}

baseClassifierSmp*baseClassifier::addTrain(double*v,int vl,seqClass*c){
	if(vl!=nFeatures){
		cmdError("Training vector of incorrect size.");
		cout << t_indent << vl << " != " << nFeatures << "\n";
		return 0;
	}
	baseClassifierSmp*t = baseClassifierSmp::create(v,vl,c);
	if(!t)return 0;
	trainingExamples.push_back(t);
	return t;
}

baseClassifierSmp*baseClassifier::addTrainV(double*v,int vl,seqClass*c,double val){
	baseClassifierSmp*t=addTrain(v,vl,c);
	t->clsE=val;
	return t;
}

double baseClassifier::apply(double*v,int vl){
	if(!trained){
		cmdWarning("Tried to apply untrained classifier.");
		return -1;
	}
	if(vl!=nFeatures){
		cmdWarning("Vector to classify is of incorrect size.");
		return -1;
	}
	return do_apply(v)-threshold;
}

bool baseClassifier::train(){
	if(trained){
		cmdError("Classifier already trained.");
		return false;
	}
	if(!do_train()){
		return false;
	}
	trained=true;
	return true;
}

bool baseClassifier::exportAnalysisData(FILE*f, char*title, char*indent){
	cmdError("Model export has not yet been implemented for this classifier.");
	return false;
}


////////////////////////////////////////////////////////////////////////////////////
// Log-odds classifier
//	Works similarly to the base classifier in the PREdictor.

logoddsClassifier::logoddsClassifier(int nf):baseClassifier(nf){
	nP=nN=0;
}

logoddsClassifier*logoddsClassifier::create(int nf){
	if(!nf){
		cmdError("No features.");
		return 0;
	}
	logoddsClassifier*r=new logoddsClassifier(nf);
	if(!r){
		outOfMemory();
		return 0;
	}
	r->weights.resize(nf);
	r->cP.resize(nf);
	r->cN.resize(nf);
	if(!r->weights.ptr||!r->cP.ptr||!r->cN.ptr){
		outOfMemory();
		delete r;
		return 0;
	}
	r->weights.fill(nf, 0);
	r->cP.fill(nf, 0);
	r->cN.fill(nf, 0);
	return r;
}

double logoddsClassifier::getWeight(int i){
	if(i<0||i>=nFeatures){
		cmdError("Invalid feature index.");
		return 0;
	}
	return weights[i];
}

bool logoddsClassifier::do_train(){
	for(baseClassifierSmp*t: trainingExamples.v){
		for(int i=0;i<nFeatures;i++){
			if(t->cls->flag){
				nP++;
				cP[i]+=t->vec[i];
			}else{
				nN++;
				cN[i]+=t->vec[i];
			}
		}
	}
	config*cfg = getConfiguration();
	for(int i=0;i<nFeatures;i++){
		switch(cfg->wmMode){
			case wmPREdictor:
				// Identical to jPREdictor
				if(cP[i]==0||cN[i]==0)cP[i]+=1,cN[i]+=1;
				weights[i]=log(cP[i])-log(double(nP))-(log(cN[i])-log(double(nN)));
				break;
			case wmZero:
				// Zero weight instead of pseudocount
				if(cP[i]==0||cN[i]==0)weights[i]=0;
				else weights[i]=log(cP[i])-log(double(nP))-(log(cN[i])-log(double(nN)));
				break;
			case wmConstant:{
				// Constant pseudocount
				double beta=cfg->loBeta;
				weights[i]=log(cP[i]+beta)-log(double(nP))-(log(cN[i]+beta)-log(double(nN)));
				break;}
			case wmPPV:{
				double TP=cP[i]/double(nP); // Predictions in positives = TP
				double FP=cN[i]/double(nN); // Predictions in negatives = FP
				if(TP==0&&FP==0){
					weights[i]=0;
				}else{
					// Weight with Positive Predictive Value (PPV).
					double PPV=TP/(TP+FP);
					weights[i]=PPV;
				}
				break;}
			case wmBiPPV:{
				double TP=cP[i]/double(nP); // Predictions in positives = TP
				double FP=cN[i]/double(nN); // Predictions in negatives = FP
				if(TP==0&&FP==0){
					weights[i]=0;
				}else{
					// score = frq * PPVpos - frq * PPVneg
					// = frq * (PPVpos - PPVneg)
					// = frq * (TP/(TP+FP) - FP/(TP+FP))
					// = frq * (TP-FP)/(TP+FP)
					// This corresponds to adding counts for the ones we think
					// are positive, and subtracting (penalizing) those that
					// we think are negatives.
					weights[i]=(TP-FP)/(TP+FP);
				}
				break;}
			default:{}
		}
	}
	return true;
}

double logoddsClassifier::do_apply(double*vec){
	double r=0;
	for(int i=0;i<nFeatures;i++){
		r+=vec[i]*weights[i];
	}
	return r;
}

void logoddsClassifier::printInfo(char*header){
	cout << t_indent << header << "\n";
	cout << t_indent << t_indent << "Weights:\n";
	if(featureNames.size() == (unsigned int)nFeatures){
		for(int i=0;i<nFeatures;i++){
			cout << t_indent << t_indent << t_indent << featureNames[i] << ": " << weights[i] << "\n";
		}
	}else{
		for(int i=0;i<nFeatures;i++){
			cout << t_indent << t_indent << t_indent << weights[i] << "\n";
		}
	}
}

bool logoddsClassifier::exportAnalysisData(FILE*f, char*title, char*indent){
	fprintf(f,"%s%s\n",indent,title);
	fprintf(f,"%s - Weights:\n",indent);
	if(featureNames.size() == (unsigned int)nFeatures){
		for(int i=0;i<nFeatures;i++){
			fprintf(f,"%s -  - %s: %.14f\n",indent,featureNames[i].c_str(),weights[i]);
		}
	}else{
		for(int i=0;i<nFeatures;i++){
			fprintf(f,"%s -  - %.14f\n",indent,weights[i]);
		}
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////////////
// Fast SVM
//	Trains SVM with LibSVM, but uses optimized classification procedures.
//	Implements feature value scaling internally.

fastSVMClassifier::fastSVMClassifier(int type,int nf,std::string _name):baseClassifier(nf){
	memset(&svmprob,0,sizeof(svm_problem));
	config*cfg=getConfiguration();
	struct svm_parameter param={
		type, // svm type
		POLY,	// kernel type
		2,		// degree
		1.0/double(nf),// gamma
		cfg->SVM_c0,// coef0
		10000,	// cache size
		1e-3,	// eps
		cfg->SVM_C,	// C
		0,		// nr_weight
		NULL,	// weight label
		NULL,	// weight
		cfg->SVM_nu,	// nu
		cfg->SVM_p,	// p
		1,		// shrinking
		0,		// probability
	};
	svmparam=param;
	svmmdl=0;
	svmvector=0;
	nFeatures=nf;
	mSVcoef=0;
	name = _name;
}

bool fastSVMClassifier::IaddTrain(double*fv,double cls){
	svmprob.l++;
	svmprob.y=(double*)realloc(svmprob.y,sizeof(double)*svmprob.l);
	if(!svmprob.y){
		svmprob.l--;
		outOfMemory();
		return false;
	}
	svmprob.x=(svm_node**)realloc(svmprob.x,sizeof(svm_node*)*svmprob.l);
	if(!svmprob.x){
		svmprob.l--;
		outOfMemory();
		return false;
	}
	svm_node*nv=(svm_node*)malloc(sizeof(svm_node)*(nFeatures+1));
	if(!nv){
		svmprob.l--;
		outOfMemory();
		return false;
	}
	memset(nv,0,sizeof(svm_node)*(nFeatures+1));
	nv[nFeatures].index=-1;
	for(int i=0;i<nFeatures;i++){
		nv[i].index=i;
		nv[i].value=fv[i];
	}
	svmprob.x[svmprob.l-1]=nv;
	svmprob.y[svmprob.l-1]=cls;
	return true;
}

fastSVMClassifier::~fastSVMClassifier(){
	if(mSVcoef){
		for(int l=0;l<svmmdl->l;l++){
			if(mSVcoef[l])free(mSVcoef[l]);
		}
		free(mSVcoef);
	}
	if(svmmdl)svm_free_and_destroy_model(&svmmdl);
	if(svmvector)free(svmvector);
	if(svmprob.x){
		for(int l=0;l<svmprob.l;l++){
			if(svmprob.x[l])free(svmprob.x[l]);
		}
		free(svmprob.x);
	}
	if(svmprob.y)free(svmprob.y);
}

fastSVMClassifier*fastSVMClassifier::create(int type,int nf,std::string _name){
	fastSVMClassifier*r=new fastSVMClassifier(type,nf,_name);
	if(!r){
		outOfMemory();
		return 0;
	}
	switch(getConfiguration()->kernel){
		case kLinear:
			r->svmparam.kernel_type=LINEAR;
			r->svmparam.degree=0;
			r->svmparam.gamma=0;
			r->svmparam.coef0=0;
			break;
		case kQuadratic:
			r->svmparam.kernel_type=POLY;
			r->svmparam.degree=2;
			break;
		case kCubic:
			r->svmparam.kernel_type=POLY;
			r->svmparam.degree=3;
			break;
		case kRBF:
			r->svmparam.kernel_type=RBF;
			break;
		default:
			cmdError("Invalid kernel.");
			delete r;
			return 0;
	}
	r->svmvector=(svm_node*)malloc(sizeof(svm_node)*(nf+1));
	if(!r->svmvector){
		outOfMemory();
		delete r;
		return 0;
	}
	memset(r->svmvector,0,sizeof(svm_node)*(nf+1));
	for(int l=0;l<nf;l++)r->svmvector[l].index=l;
	r->svmvector[nf].index=-1;
	r->vMin.resize(nf);
	r->vMax.resize(nf);
	r->SVcoef.resize(nf);
	if(!r->vMin.ptr||!r->vMax.ptr||!r->SVcoef.ptr){
		outOfMemory();
		delete r;
		return 0;
	}
	return r;
}

void fastSVMClassifier::scaleVector(svm_node*n){
	if(!n)return;
	for(int x=0;x<nFeatures;x++){
		double range=vMax[x]-vMin[x];
		double s=0;
		if(range!=0)s=2.0/range;
		n[x].value=(n[x].value-vMin[x]-range*0.5)*s;
	}
}

void fastSVMClassifier::scaleVectorDoubleFloat(double*in,float*out){
	if(!in||!out)return;
	for(int x=0;x<nFeatures;x++,in++,out++){
		double range=vMax[x]-vMin[x];
		double s=0;
		if(range!=0)s=2.0/range;
		*out=float(((*in)-vMin[x]-range*0.5)*s);
	}
}

bool fastSVMClassifier::do_train(){
	// Make LibSVM vectors from stored double vectors
	for(baseClassifierSmp*t: trainingExamples.v){
		if(svmparam.svm_type==ONE_CLASS&&!t->cls->flag)continue;
		if(!IaddTrain(t->vec,t->clsE)){
			return false;
		}
	}
	// Scale and train
	for(int l=0;l<svmprob.l;l++){
		svm_node*v=svmprob.x[l];
		for(int x=0;x<nFeatures;x++){
			if(!l||v[x].value>vMax[x]){
				vMax[x]=v[x].value;
			}
			if(!l||v[x].value<vMin[x]){
				vMin[x]=v[x].value;
			}
		}
	}
	for(int l=0;l<svmprob.l;l++)scaleVector(svmprob.x[l]);
	svmmdl=svm_train(&svmprob,&svmparam);
	if(!svmmdl){
		cmdError("LibSVM did not return a classifier");
		return false;
	}
	if(!svmmdl->sv_coef){
		svmparam.kernel_type=-1;
		return true;
	}
	
	if(svmparam.kernel_type==LINEAR){
		double*sv_coef=svmmdl->sv_coef[0];
		for(int l=0;l<nFeatures;l++){
			double w=0;
			for(int i=0;i<svmmdl->l;i++){
				svm_node*sv=svmmdl->SV[i];
				w+=sv_coef[i]*sv[l].value;
			}
			SVcoef[l]=w;
		}
	}else if(svmparam.kernel_type==POLY){
		fVec.resize(nFeatures);
		if(!fVec.ptr){
			outOfMemory();
			return false;
		}
		mSVcoef=(float**)malloc(sizeof(float*)*svmmdl->l);
		if(!mSVcoef){
			outOfMemory();
			return false;
		}
		svm_node**sv=svmmdl->SV;
		for(int l=0;l<svmmdl->l;l++,sv++){
			mSVcoef[l]=(float*)malloc(sizeof(float)*nFeatures);
			if(!mSVcoef[l]){
				outOfMemory();
				return false;
			}
			svm_node*csv=*sv;
			float*v=mSVcoef[l];
			for(int x=0;x<nFeatures;x++,csv++,v++){
				*v=(float)csv->value;
			}
		}
	}
	
	return true;
}

double fastSVMClassifier::do_apply(double*fv){
	if(!svmmdl)return 0;
	if(svmparam.kernel_type==-1){
		return 0;
	}else if(svmparam.kernel_type==LINEAR){
		for(int i=0;i<nFeatures;i++){
			svmvector[i].value=fv[i];
		}
		scaleVector(svmvector);
		double r=0;
		for(int l=0;l<nFeatures;l++){
			r+=svmvector[l].value*SVcoef[l];
		}
		return r-svmmdl->rho[0];
		
	}else if(svmparam.kernel_type==POLY){
		scaleVectorDoubleFloat(fv,fVec.ptr);
		double r=0;
		double kv;
		float dot;
		float**sv=mSVcoef;
		double*coef=svmmdl->sv_coef[0];
		for(int l=0;l<svmmdl->l;l++,coef++,sv++){
			dot=0;
			int x=0;
			float*vc=fVec.ptr;
			float*svc=*sv;
			for(;x<(nFeatures&((-1)^7));x+=8,vc+=8,svc+=8){
				dot+=vc[0]*svc[0];
				dot+=vc[1]*svc[1];
				dot+=vc[2]*svc[2];
				dot+=vc[3]*svc[3];
				dot+=vc[4]*svc[4];
				dot+=vc[5]*svc[5];
				dot+=vc[6]*svc[6];
				dot+=vc[7]*svc[7];
			}
			for(;x<nFeatures;x++,vc++,svc++){
				dot+=vc[0]*svc[0];
			}

			kv=svmparam.gamma*dot+svmparam.coef0;
			if(svmparam.degree==2)r+=kv*kv*coef[0];
			else if(svmparam.degree==3)r+=kv*kv*kv*coef[0];
		}
		return r-svmmdl->rho[0];
	}else if(svmparam.kernel_type==RBF){
		for(int i=0;i<nFeatures;i++){
			svmvector[i].value=fv[i];
		}
		scaleVector(svmvector);
		double r=0;
		double*coef=svmmdl->sv_coef[0];
		svm_node**sv=svmmdl->SV,*csv,*mv;
		double qdist,tmp;
		for(int l=0;l<svmmdl->l;l++,coef++,sv++){
			csv=sv[0];
			mv=svmvector;
			qdist=0;
			for(int x=0;x<nFeatures;x++,csv++,mv++){
				tmp=mv->value-csv->value;
				qdist+=tmp*tmp;
			}
			r+=exp(-svmparam.gamma*qdist)*coef[0];
		}
		return r-svmmdl->rho[0];
	}
	cmdWarning("Unsupported kernel");
	return 0;
}

void fastSVMClassifier::printInfo(char*header){
	if(name.length() > 0){
		cout << t_indent << header << " - " << name << "\n";
	}else{
		cout << t_indent << header << "\n";
	}
	config*cfg=getConfiguration();
	cout << t_indent << t_indent << "Kernel: " << getKernelName(cfg->kernel) << "\n";
	cout << t_indent << t_indent << "Type: ";
	switch(svmparam.svm_type){
	case C_SVC: cout << "C_SVC\n"; break;
	case NU_SVC: cout << "NU_SVC\n"; break;
	case ONE_CLASS: cout << "ONE_CLASS\n"; break;
	case EPSILON_SVR: cout << "EPSILON_SVR\n"; break;
	case NU_SVR: cout << "NU_SVR\n"; break;
	default: cout << "Invalid\n"; break;
	}
	if(svmparam.svm_type==C_SVC||svmparam.svm_type==EPSILON_SVR){
		cout << t_indent << t_indent << "C: " << cfg->SVM_C << "\n";
	}
	if(svmparam.svm_type==NU_SVC||svmparam.svm_type==NU_SVR){
		cout << t_indent << t_indent << "nu: " << cfg->SVM_nu << "\n";
	}
	if(cfg->kernel!=kLinear){
		cout << t_indent << t_indent << "gamma: " << cfg->SVM_gamma << "\n";
		cout << t_indent << t_indent << "c0: " << cfg->SVM_c0 << "\n";
	}
	if(svmparam.svm_type==EPSILON_SVR){
		cout << t_indent << t_indent << "p: " << cfg->SVM_p << "\n";
	}
	cout << t_indent << t_indent << "# SV: " << svmmdl->l << "\n";
	
	if(svmparam.kernel_type == LINEAR){
		
		cout << t_indent << t_indent << "Model weights (linear SVM)\n";
		double modelBias = 0.0;
		for(int l=0;l<nFeatures;l++){
			double range = vMax[l] - vMin[l];
			double s = 0;
			if(range != 0.0) s = 2.0 / range;
			modelBias -= (vMin[l] * s + 1.0) * SVcoef[l];
		}
		modelBias -= svmmdl->rho[0];
		cout << t_indent << t_indent << t_indent << "Model bias: " << modelBias << "\n";
		for(int l=0;l<nFeatures;l++){
			double range = vMax[l] - vMin[l];
			double s = 0;
			if(range != 0.0) s = 2.0 / range;
			// Show coefficient
			double weight = s * SVcoef[l];
			cout << t_indent << t_indent << t_indent << "Feature weight - " << featureNames[l] << ": " << weight << "\n";
		}
		
	}else if(svmparam.kernel_type == POLY && svmparam.degree == 2 && svmparam.coef0 == 0.0){
		
		cout << t_indent << t_indent << "Model weights (quadratic SVM)\n";
		
		double gamma = svmparam.gamma;
		
		// Quadratic component
		for(int i = 0; i < nFeatures; i++){
			double alpha_i = vMin[i];
			double beta_i = vMax[i];
			if(alpha_i == beta_i) continue;
			for(int j = 0; j < nFeatures; j++){
				double alpha_j = vMin[j];
				double beta_j = vMax[j];
				if(alpha_j == beta_j) continue;
				double w_ij = 0.0;
				for(int x = 0; x < svmmdl->l; x++){
					double c = svmmdl->sv_coef[0][x];
					double v_i = svmmdl->SV[x][i].value;
					double v_j = svmmdl->SV[x][j].value;
					//
					w_ij += gamma * gamma * c * ( (2.0*v_i)/(beta_i - alpha_i) ) * ( (2.0*v_j)/(beta_j - alpha_j) );
					//
				}
				cout << t_indent << t_indent << t_indent << "Feature pair weight - " << featureNames[i] << " / " << featureNames[j] << ": " << w_ij << "\n";
			}
		}
		
		// Linear component
		for(int i = 0; i < nFeatures; i++){
			double w_i = 0.0;
			//
			double alpha_i = vMin[i];
			double beta_i = vMax[i];
			if(alpha_i == beta_i) continue;
			for(int j = 0; j < nFeatures; j++){
				double alpha_j = vMin[j];
				double beta_j = vMax[j];
				if(alpha_j == beta_j) continue;
				for(int x = 0; x < svmmdl->l; x++){
					double c = svmmdl->sv_coef[0][x];
					double v_i = svmmdl->SV[x][i].value;
					double v_j = svmmdl->SV[x][j].value;
					//
					w_i -= 2.0 * gamma * gamma * c * ( (2.0*v_i)/(beta_i - alpha_i) ) * ( (v_j*(alpha_j + beta_j))/(beta_j - alpha_j) );
				}
			}
			//
			cout << t_indent << t_indent << t_indent << "Single feature weight - " << featureNames[i] << ": " << w_i << "\n";
		}
		
		// Bias
		double modelBias = 0.0;
		for(int x = 0; x < svmmdl->l; x++){
			double c = svmmdl->sv_coef[0][x];
			
			double Q = 0.0;
			for(int i = 0; i < nFeatures; i++){
				double alpha_i = vMin[i];
				double beta_i = vMax[i];
				if(alpha_i == beta_i) continue;
				double v_i = svmmdl->SV[x][i].value;
				Q += (v_i*(alpha_i + beta_i))/(beta_i - alpha_i);
			}
			modelBias += gamma * gamma * c * Q * Q;
		}
		modelBias -= svmmdl->rho[0];
		cout << t_indent << t_indent << t_indent << "Model bias: " << modelBias << "\n";
		
	}

}

bool fastSVMClassifier::exportAnalysisData(FILE*f, char*title, char*indent){
	if(name.length() > 0){
		fprintf(f, "%s%s - %s\n", indent, title, name.c_str());
	}else{
		fprintf(f, "%s%s\n", indent, title);
	}
	config*cfg=getConfiguration();
	fprintf(f, "%s - Kernel: %s\n", indent, getKernelName(cfg->kernel));
	fprintf(f, "%s - Type: ", indent);
	switch(svmparam.svm_type){
	case C_SVC: fprintf(f, "C_SVC\n"); break;
	case NU_SVC: fprintf(f, "NU_SVC\n"); break;
	case ONE_CLASS: fprintf(f, "ONE_CLASS\n"); break;
	case EPSILON_SVR: fprintf(f, "EPSILON_SVR\n"); break;
	case NU_SVR: fprintf(f, "NU_SVR\n"); break;
	default: fprintf(f, "Invalid\n"); break;
	}
	if(svmparam.svm_type==C_SVC||svmparam.svm_type==EPSILON_SVR){
		fprintf(f, "%s - C: %.14f\n", indent, cfg->SVM_C);
	}
	if(svmparam.svm_type==NU_SVC||svmparam.svm_type==NU_SVR){
		fprintf(f, "%s - nu: %.14f\n", indent, cfg->SVM_nu);
	}
	if(cfg->kernel!=kLinear){
		fprintf(f, "%s - gamma: %.14f\n", indent, cfg->SVM_gamma);
		fprintf(f, "%s - c0: %.14f\n", indent, cfg->SVM_c0);
	}
	if(svmparam.svm_type==EPSILON_SVR){
		fprintf(f, "%s - p: %.14f\n", indent, cfg->SVM_p);
	}
	fprintf(f, "%s - # SV: %d\n", indent, svmmdl->l);
	if(svmparam.kernel_type == LINEAR){
		fprintf(f, "%s - Model weights (linear SVM, scaled)\n", indent);
		double modelBias = 0.0;
		for(int l=0;l<nFeatures;l++){
			double range = vMax[l] - vMin[l];
			double s = 0;
			if(range != 0.0) s = 2.0 / range;
			modelBias -= (vMin[l] * s + 1.0) * SVcoef[l];
		}
		modelBias -= svmmdl->rho[0];
		fprintf(f, "%s - Model bias: %.14f\n", indent, modelBias);
		for(int l=0;l<nFeatures;l++){
			double range = vMax[l] - vMin[l];
			double s = 0;
			if(range != 0.0) s = 2.0 / range;
			// Show coefficient
			double weight = s * SVcoef[l];
			fprintf(f, "%s -  - Feature weight - %s: %.14f\n", indent, featureNames[l].c_str(), weight);
		}
	}else if(svmparam.kernel_type == POLY && svmparam.degree == 2 && svmparam.coef0 == 0.0){
		fprintf(f, "%s - Model weights (quadratic SVM, scaled)\n", indent);
		
		double gamma = svmparam.gamma;
		
		// Quadratic component
		for(int i = 0; i < nFeatures; i++){
			double alpha_i = vMin[i];
			double beta_i = vMax[i];
			if(alpha_i == beta_i) continue;
			for(int j = 0; j < nFeatures; j++){
				double alpha_j = vMin[j];
				double beta_j = vMax[j];
				if(alpha_j == beta_j) continue;
				double w_ij = 0.0;
				for(int x = 0; x < svmmdl->l; x++){
					double c = svmmdl->sv_coef[0][x];
					double v_i = svmmdl->SV[x][i].value;
					double v_j = svmmdl->SV[x][j].value;
					//
					w_ij += gamma * gamma * c * ( (2.0*v_i)/(beta_i - alpha_i) ) * ( (2.0*v_j)/(beta_j - alpha_j) );
					//
				}
				fprintf(f, "%s -  - Feature pair weight - %s / %s: %.14f\n", indent, featureNames[i].c_str(), featureNames[j].c_str(), w_ij);
			}
		}
		
		// Linear component
		for(int i = 0; i < nFeatures; i++){
			double w_i = 0.0;
			//
			double alpha_i = vMin[i];
			double beta_i = vMax[i];
			if(alpha_i == beta_i) continue;
			for(int j = 0; j < nFeatures; j++){
				double alpha_j = vMin[j];
				double beta_j = vMax[j];
				if(alpha_j == beta_j) continue;
				for(int x = 0; x < svmmdl->l; x++){
					double c = svmmdl->sv_coef[0][x];
					double v_i = svmmdl->SV[x][i].value;
					double v_j = svmmdl->SV[x][j].value;
					//
					w_i -= 2.0 * gamma * gamma * c * ( (2.0*v_i)/(beta_i - alpha_i) ) * ( (v_j*(alpha_j + beta_j))/(beta_j - alpha_j) );
				}
			}
			//
			fprintf(f, "%s -  - Single feature weight - %s: %.14f\n", indent, featureNames[i].c_str(), w_i);
		}
		
		// Bias
		double modelBias = 0.0;
		for(int x = 0; x < svmmdl->l; x++){
			double c = svmmdl->sv_coef[0][x];
			
			double Q = 0.0;
			for(int i = 0; i < nFeatures; i++){
				double alpha_i = vMin[i];
				double beta_i = vMax[i];
				if(alpha_i == beta_i) continue;
				double v_i = svmmdl->SV[x][i].value;
				Q += (v_i*(alpha_i + beta_i))/(beta_i - alpha_i);
			}
			modelBias += gamma * gamma * c * Q * Q;
		}
		modelBias -= svmmdl->rho[0];
		fprintf(f, "%s -  - Model bias:: %.14f\n", indent, modelBias);
		
		{
		fprintf(f, "%s - Model weights (quadratic SVM, unscaled)\n", indent);
		
		double gamma = svmparam.gamma;
		
		// Quadratic component
		for(int i = 0; i < nFeatures; i++){
			for(int j = 0; j < nFeatures; j++){
				double w_ij = 0.0;
				for(int x = 0; x < svmmdl->l; x++){
					double c = svmmdl->sv_coef[0][x];
					double v_i = svmmdl->SV[x][i].value;
					double v_j = svmmdl->SV[x][j].value;
					//
					w_ij += gamma * gamma * c * v_i * v_j;
					//
				}
				fprintf(f, "%s -  - Feature pair weight - %s / %s: %.14f\n", indent, featureNames[i].c_str(), featureNames[j].c_str(), w_ij);
			}
		}
		
		// Bias
		modelBias = -svmmdl->rho[0];
		fprintf(f, "%s -  - Model bias:: %.14f\n", indent, modelBias);
		}
	}
	
	return true;
}


////////////////////////////////////////////////////////////////////////////////////
// Multi-class SVM

MultiClassSVM::MultiClassSVM(int nf,int _svmt,std::string _name){
	trained=false;
	nFeatures=nf;
	svmtype=_svmt;
	name = _name;
}

MultiClassSVM*MultiClassSVM::create(int _svmtype,int nfeatures,std::string _name){
	MultiClassSVM*r=new MultiClassSVM(nfeatures,_svmtype,_name);
	if(!r){
		outOfMemory();
		return 0;
	}
	return r;
}

bool MultiClassSVM::addTrain(double*v,int vl,seqClass*c,double cE){
	if(trained){
		cmdError("Classifier already trained.");
		return false;
	}
	if(vl!=nFeatures){
		cmdError("Training vector of incorrect size.");
		cout << t_indent << vl << " != " << nFeatures << "\n";
		return false;
	}
	baseClassifierSmp*t = baseClassifierSmp::create(v,vl,c);
	if(!t)return 0;
	t->clsE = cE;
	trainingExamples.push_back(t);
	return true;
}

bool MultiClassSVM::train(){
	if(trained){
		cmdError("Classifier already trained.");
		return false;
	}
	// Note classes.
	for(baseClassifierSmp*t: trainingExamples.v){
		seqClass*cls=t->cls;
		// Skip any classes already noted.
		bool found=false;
		for(int y=0;y<(int)classes.size();y++){
			if(classes[y].cls==cls){
				found=true;
				break;
			}
		}
		if(found)continue;
		// Note it.
		MultiClassClass c={0,0};
		c.cls=cls;
		classes.push_back(c);
	}
	// Note class pairs.
	for(int x=0;x<(int)classes.size()-1;x++){
		for(int y=x+1;y<(int)classes.size();y++){
			// Add pair and make classifier.
			MultiClassClassPair p;
			p.clsP=&classes[x];
			p.clsN=&classes[y];
			char cName[256];
			snprintf(cName, sizeof(cName), "%s - Class boundary model - %s (+) vs. %s (-)", name.c_str(), p.clsP->cls->name.c_str(), p.clsN->cls->name.c_str());
			fastSVMClassifier*cls=fastSVMClassifier::create(svmtype,nFeatures,std::string(cName));
			if(!cls)return false;
			if(featureNames.size() > 0){
				for(unsigned int z = 0; z < featureNames.size(); z++){
					cls->featureNames.push_back(std::string(featureNames[z]));
				}
			}
			classifiers.push_back(cls);
			p.classifier=cls;
			borders.push_back(p);
			// Train.
			for(baseClassifierSmp*t: trainingExamples.v){
				if(t->cls!=p.clsP->cls
				&&t->cls!=p.clsN->cls){
					continue;
				}
				if(!cls->addTrainV(t->vec,t->vecl,t->cls,t->cls==p.clsP->cls?1.0:-1.0)){
					return false;
				}
			}
			if(!cls->train())return false;
		}
	}
	trained=true;
	return true;
}

seqClass*MultiClassSVM::apply(double*vec,int vecl){
	if(!trained){
		cmdError("Tried to apply untrained classifier.");
		return 0;
	}
	if(vecl!=nFeatures){
		cmdWarning("Vector to classify is of incorrect size.");
		return 0;
	}
	// Have classifiers vote for the correct class.
	for(auto&c:classes){
		c.votes=0;
	}
	for(auto&p:borders){
		double v=p.classifier->apply(vec,vecl);
		if(v>0)p.clsP->votes++;
		else p.clsN->votes++;
	}
	// Find the (first) class with the most votes.
	int maxVotes=-1;
	seqClass*r=0;
	for(auto&c:classes){
		if(c.votes>maxVotes){
			maxVotes=c.votes;
			r=c.cls;
		}
	}
	return r;
}

void MultiClassSVM::printInfo(char*header){
	if(name.length() > 0){
		cout << t_indent << header << " - " << name << "\n";
	}else{
		cout << t_indent << header << "\n";
	}
	cout << t_indent << t_indent << "Training examples: " << (trainingExamples.v.size()) << "\n";
	cout << t_indent << t_indent << "Classes: " << classes.size() << "\n";
	cout << t_indent << t_indent << "Class borders: " << borders.size() << "\n";
	for(auto&p:borders){
		p.classifier->printInfo((char*)"Boundary");
	}
}

bool MultiClassSVM::exportAnalysisData(FILE*f, char*title, char*indent){
	if(name.length() > 0){
		fprintf(f,"%s%s - %s\n",indent,title,name.c_str());
	}else{
		fprintf(f,"%s%s\n",indent,title);
	}
	fprintf(f,"%s - Training examples: %d\n",indent,(int)trainingExamples.v.size());
	fprintf(f,"%s - Classes: %d\n",indent,(int)classes.size());
	fprintf(f,"%s - Class borders: %d\n",indent,(int)borders.size());
	for(auto&p:borders){
		if(!p.classifier->exportAnalysisData(f, (char*)"Boundary", (char*)(std::string(indent)+" - ").c_str() ))
			return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
// Random Forest
//	Trains RF with Ranger.

RFClassifier::RFClassifier(int nf):baseClassifier(nf){
}

RFClassifier*RFClassifier::create(int nf, std::string name){
	if(!nf){
		cmdError("No features.");
		return 0;
	}
	RFClassifier*r=new RFClassifier(nf);
	if(!r){
		outOfMemory();
		return 0;
	}
	r->rf.ptr = new RangerRandomForest();
	r->name = name;
	return r;
}

bool RFClassifier::do_train(){
	// Fill in data values
	std::unique_ptr<RangerData> data{};
	data = ranger::make_unique<RangerData>(); //new ranger::DataDouble();
	data->setDataT(nFeatures, trainingExamples.v);
	
	// Train model
	std::vector<std::string> catvars;
	std::vector<double> regfac;
	std::vector<double> sample_fraction_vector = { ranger::DEFAULT_SAMPLE_FRACTION_REPLACE };
	rf.ptr->init(
		ranger::MemoryMode::MEM_DOUBLE,
		std::move(data),
		(ranger::uint)0, // uint mtry,
		std::string(""), //std::string output_prefix,
		getConfiguration()->RF_nTrees, //uint num_trees,
		(ranger::uint)0, //uint seed,
		(ranger::uint)getConfiguration()->nThreads, //uint num_threads,
		//ranger::ImportanceMode::IMP_NONE, //ImportanceMode importance_mode,
		ranger::ImportanceMode::IMP_PERM_CASEWISE, //ImportanceMode importance_mode,
		(ranger::uint)0, //uint min_node_size,
		false, //bool prediction_mode,
		true, //bool sample_with_replacement,
		catvars, //const std::vector<std::string>& unordered_variable_names,
		false, //bool memory_saving_splitting,
		//ranger::SplitRule::LOGRANK,//SplitRule splitrule,
		ranger::SplitRule::HELLINGER,//SplitRule splitrule,
		false, //bool predict_all,
		sample_fraction_vector, //std::vector<double>& sample_fraction,
		ranger::DEFAULT_ALPHA, //double alpha,
		ranger::DEFAULT_MINPROP, //double minprop,
		false, //bool holdout,
		ranger::PredictionType::RESPONSE, //PredictionType prediction_type,
		ranger::DEFAULT_NUM_RANDOM_SPLITS, //uint num_random_splits,
		false, //bool order_snps,
		ranger::DEFAULT_MAXDEPTH, //uint max_depth,
		regfac,//const std::vector<double>& regularization_factor,
		false //bool regularization_usedepth
	);
	rf.ptr->train();
	
	
	return true;
}

double RFClassifier::do_apply(double*vec){
	return rf.ptr->predictVec(vec);
}

void RFClassifier::printInfo(char*header){
	cout << t_indent << header << "\n";
}

////////////////////////////////////////////////////////////////////////////////////
// Linear Discriminant Analysis
#ifdef USE_SHOGUN

LDAClassifier::LDAClassifier(int nf):baseClassifier(nf){
}

LDAClassifier*LDAClassifier::create(int nf, std::string name){
	if(!nf){
		cmdError("No features.");
		return 0;
	}
	LDAClassifier*r=new LDAClassifier(nf);
	if(!r){
		outOfMemory();
		return 0;
	}
	r->name = name;
	return r;
}

bool LDAClassifier::do_train(){
	shogun::init_shogun_with_defaults();
	int ntrain = trainingExamples.v.size();
	shogun::SGVector<double> lvec(ntrain);
	shogun::SGMatrix<double> fmat(nFeatures, ntrain);
	int row = 0;
	for(baseClassifierSmp*t: trainingExamples.v){
		for(int i=0; i<nFeatures; i++)
			fmat(i, row) = t->vec[i];
		lvec[row] = t->cls->flag ? 1 : -1;
		row++;
	}
	shogun::CBinaryLabels*lbl = new shogun::CBinaryLabels(lvec);
	shogun::CDenseFeatures<double>*fv = new shogun::CDenseFeatures<double>(fmat);
	LDA.ptr = new shogun::CLDA(0.0001, fv, lbl);
	LDA.ptr->train();
	return true;
}

double LDAClassifier::do_apply(double*vec){
	shogun::SGMatrix<double> fmat(nFeatures, 1);
	for(int i=0; i<nFeatures; i++)
		fmat(i, 0) = vec[i];
	shogun::CDenseFeatures<double>*fv = new shogun::CDenseFeatures<double>(fmat);
	autodelete<shogun::CBinaryLabels> pred(LDA.ptr->apply_binary(fv));
	return pred.ptr->get_value(0);
}

void LDAClassifier::printInfo(char*header){
	cout << t_indent << header << "\n";
}
#endif

