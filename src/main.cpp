////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
////////////////////////////////////////////////////////////////////////////////////
// Basics

#include "common.hpp"
#include "./lib/rapidxml-1.13/rapidxml.hpp"
#include "./lib/libsvm-3.17/svm.h"
using namespace rapidxml;

////////////////////////////////////////////////////////////////////////////////////
// Program parts

#include "aux.hpp"
#include "config.hpp"
#include "sequences.hpp"
#include "sequencelist.hpp"
#include "validation.hpp"
#include "motifs.hpp"
#include "models/features.hpp"
#include "models/baseclassifier.hpp"
#include "models/sequenceclassifier.hpp"
#include "models/svmmocca.hpp"
#include "models/cpredictor.hpp"
#include "models/dummypredictor.hpp"
#include "models/seqsvm.hpp"

////////////////////////////////////////////////////////////////////////////////////
// File list

/*
fileListFile
	Structure for file list entries.
*/
typedef struct{
	char*name;
	char*path;
}fileListFile;

std::vector<fileListFile> fileList;

bool registerFile(char*name,char*path){
	fileListFile f = fileListFile {
		name,
		path
	};
	fileList.push_back(f);
	return true;
}

void printRegisteredFiles(){
	cmdSection("Files");
	if(fileList.size() > 0){
		for(auto& f: fileList)
			cout << t_indent << f.name << t_indent << f.path << "\n";
	}else{
		cout << t_indent << "None\n";
	}
}

////////////////////////////////////////////////////////////////////////////////////
// Main

/*
constructClassifier
	Constructs a classifier and trains it.
*/
sequenceClassifier*constructClassifier(motifList*motifs,featureSet*features,seqList*trainseq){
	sequenceClassifier*cls=0;
	config*cfg=getConfiguration();
	switch(cfg->classifier){
		case cSVMMOCCA:cls=SVMMOCCA::create(motifs,cfg->svmtype);break;
		case cCPREdictor:cls=CPREdictor::create(motifs);break;
		case cDummyPREdictor:cls=DummyPREdictor::create(motifs);break;
		case cSEQSVM:cls=SEQSVM::create(motifs,features,cfg->svmtype);break;
		default:cmdError("Invalid classifier.");return 0;
	}
	if(!cls)return 0;
	if(!trainseq||!trainseq->nseq){
		cmdError("No training sequences specified.");
		delete cls;
		return 0;
	}
	if(!cls->train(trainseq)){
		delete cls;
		return 0;
	}
	return cls;
}

/*
runPipeline
	Runs the main application pipeline.
*/
bool runPipeline(motifList*&motifs,featureSet*&features,seqList*&trainseq,seqList*valseq){
	cout << sepline;
	
	autodelete<sequenceClassifier>cls((sequenceClassifier*)0);
	cls.ptr=constructClassifier(motifs,features,trainseq);
	config*cfg=getConfiguration();
	if(!cls.ptr){
		return false;
	}
	// Basic pipeline
	{
		if(cfg->classifier == cSEQSVM)
			features->printInfo();
		cmdSection("Classifier");
		cls.ptr->printInfo();
		cmdSection("Validation");
		if(cfg->validate){
			autofree<validationPair>vp((validationPair*)0);
			int nvp=0;
			if(!cls.ptr->getValidationTable(trainseq,vp.ptr,nvp))return false;
			cout << t_indent << "Training set\n";
			printValidationMeasures(vp,nvp,cls.ptr->threshold);
		}
		if(valseq->nseq){
			autofree<validationPair>vp((validationPair*)0);
			int nvp=0;
			if(!cls.ptr->getValidationTable(valseq,vp.ptr,nvp))return false;
			cout << t_indent << "Validation set\n";
			if(cfg->validate)printValidationMeasures(vp,nvp,cls.ptr->threshold);
			if(cfg->outSCVal)if(!saveVPairTable(cfg->outSCVal,vp.ptr,nvp))return false;
		}
		if(cfg->inFASTA&&cfg->outWig){
			cmdSection("FASTA scoring");
			cls.ptr->applyFASTA(cfg->inFASTA,cfg->outWig);
		}
		if(cfg->inFASTA&&cfg->outCoreSequence){
			cmdSection("FASTA scoring");
			cls.ptr->predictCoreSequence(cfg->inFASTA,cfg->outCoreSequence);
		}
		if(cfg->CAnalysisExportPath){
			cmdSection("Classifier analysis export");
			cls.ptr->exportAnalysisData(cfg->CAnalysisExportPath);
		}
	}
	cout << sepline;
	return true;
}

/*
print_help
	Outputs help message.
*/
void print_help(){
	cout << " Usage:\n";
	char argFmtStr[] = " %35s - %s\n";
	printf(argFmtStr, "-C:SVM-MOCCA", "Use SVM-MOCCA.");
	printf(argFmtStr, "-C:SVM-MOCCA:C-SVC", "Use SVM-MOCCA with C-SVC.");
	printf(argFmtStr, "-C:SVM-MOCCA:nu-SVC", "Use SVM-MOCCA with nu-SVC.");
	printf(argFmtStr, "-C:CPREdictor", "Use CPREdictor.");
	printf(argFmtStr, "-C:DummyPREdictor", "Use dummy PREdictor.");
	printf(argFmtStr, "-k:linear", "Use linear kernel.");
	printf(argFmtStr, "-k:quadratic", "Use quadratic kernel.");
	printf(argFmtStr, "-k:cubic", "Use cubic kernel.");
	printf(argFmtStr, "-k:RBF", "Use radial basis function kernel.");
	printf(argFmtStr, "-wSize SIZE", "Sets the window size to SIZE.");
	printf(argFmtStr, "-wStep SIZE", "Sets the window step size to SIZE.");
	printf(argFmtStr, "-motif:IUPAC NAME MOTIF MISMATCHES", "Adds an IUPAC-motif with the name NAME, and sequence MOTIF,");
	printf(argFmtStr, "", "with MISMATCHES mismatches allowed.");
	printf(argFmtStr, "-motif:XML PATH", "Adds motifs specified in an XML-file PATH.");
	printf(argFmtStr, "-f:MOCCA:nOcc", "Adds motif occurrence frequency features to SVM-MOCCA.");
	printf(argFmtStr, "-f:MOCCA:DNT", "Adds dinucleotide features to SVM-MOCCA.");
	printf(argFmtStr, "-train:FASTA PATH CLASS MODE", "Adds a training sequence file.");
	printf(argFmtStr, "", "PATH: Path to FASTA file.");
	printf(argFmtStr, "", "CLASS: A class ID, defined with `-class`, or one of the");
	printf(argFmtStr, "", "pre-specified binary classes: '+' for positive or '-'");
	printf(argFmtStr, "", "for negative.");
	printf(argFmtStr, "", "MODE: Can be \"win\", for training with all windows within");
	printf(argFmtStr, "", "each training sequence file, or \"full\", for training with");
	printf(argFmtStr, "", "the full sequences.");
	printf(argFmtStr, "-validate:FASTA PATH CLASS", "Adds a validation sequence file.");
	printf(argFmtStr, "", "PATH: Path to FASTA file.");
	printf(argFmtStr, "", "CLASS: A class ID, defined with \"-class\", or one of the");
	printf(argFmtStr, "", "pre-specified binary classes: \"+\" for positive or \"-\"");
	printf(argFmtStr, "", "for negative.");
	printf(argFmtStr, "-validate:outSCTable PATH", "Writes scores and classes for each validation sequence to");
	printf(argFmtStr, "", "a tab-separated table file, PATH.");
	printf(argFmtStr, "-class NAME VALUE FLAG", "Registers a sequence class with name NAME, value/ID VALUE,");
	printf(argFmtStr, "", "and flag FLAG (\"+\" for positive, or \"-\" for negative)");
	cout << sepline;
}

/*
parse_arg
	Parses application arguments.
	Also fills in motif list, feature set and training set.
*/
bool parse_arg(int _argc,char**argv,motifList*ml,featureSet*features,seqList*trainseq,seqList*valseq){
	if(!ml||!trainseq||!valseq)return false;
	config*cfg=getConfiguration();
	// Pass 1
	char**cargv=argv;
	int argc=_argc;
	for(int l=0;l<argc;l++,cargv++){
		char*a=*cargv;
		// Check them against the valid arguments.
		if(!strcmp(a,"-class")){
			if(l>=argc-3){
				argSyntaxError();
				return false;
			}
			bool cls;
			if(cargv[3][0]=='+'&&!cargv[3][1])cls=true;
			else if(cargv[3][0]=='-'&&!cargv[3][1])cls=false;
			else{
				argSyntaxError();
				return false;
			}
			if(!registerSeqClass(strtod(cargv[2],0),cargv[1],cls)){
				return false;
			}
			cargv[0]=cargv[1]=cargv[2]=cargv[3]=0;
			cargv+=3,argc-=3;
		}else if(!strcmp(a,"-seed")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->randSeed=(int)strtol(cargv[1],0,10);
			srand(cfg->randSeed);
			cargv[0]=cargv[1]=0;
			cargv++,argc--;
		}else if(!strcmp(a,"-h")||!strcmp(a,"-help")){
			print_help();
			return false;
		}
	}
	// Add default, binary classes if none were specified.
	if(!nSequenceClasses()){
		if(!registerSeqClass(1,(char*)"Positive",true)){
			return false;
		}
		if(!registerSeqClass(-1,(char*)"Negative",false)){
			return false;
		}
	}
	// Pass 2
	cargv=argv;
	argc=_argc;
	for(int l=0;l<argc;l++,cargv++){
		char*a=*cargv;
		// Skip arguments from the first pass.
		if(!a)continue;
		// Check them against the valid arguments.
		if(!strcmp(a,"-motif:IUPAC")){
			if(l>=argc-3){
				argSyntaxError();
				return false;
			}
			if(!ml->addIUPACMotif(cargv[1],cargv[2],(int)strtol(cargv[3],0,10))){
				return false;
			}
			cargv+=3,argc-=3;
		}else if(!strcmp(a,"-motif:kmer")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			if(!ml->addKMers((int)strtol(cargv[1],0,10))){
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-motif:Random")){
			if(l>=argc-2){
				argSyntaxError();
				return false;
			}
			if(!ml->addRandom((int)strtol(cargv[1],0,10),(int)strtol(cargv[2],0,10))){
				return false;
			}
			cargv+=2,argc-=2;
		}else if(!strcmp(a,"-motif:XML")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			if(!ml->addMotifsFromXML(cargv[1])){
				return false;
			}
			if(!registerFile((char*)"Motif XML",cargv[1])){
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-motif:FSM")){
			cfg->useFSM=true;
		}else if(!strcmp(a,"-motif:d:centers")){
			cfg->distanceMode=dmCenters;
		}else if(!strcmp(a,"-motif:d:between")){
			cfg->distanceMode=dmBetween;
		}else if(!strcmp(a,"-motif:d:noOverlap")){
			cfg->motifPairsCanOverlap=false;
		}else if(!strcmp(a,"-motif:d:overlap")){
			cfg->motifPairsCanOverlap=true;
			
		}else if(!strcmp(a,"-f:nOcc")){
			if(!features->addFeature(featureType_nOcc,featureMotif_All,0,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:nPair")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			if(!features->addFeature(featureType_nPair,featureMotif_All,featureMotif_All,featureMotif_All,strtod(cargv[1],0),0,1,0)){
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-f:nOccPair")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			if(!features->addFeature(featureType_nOccPair,featureMotif_All,featureMotif_All,featureMotif_All,strtod(cargv[1],0),0,1,0)){
				return false;
			}
			cargv++,argc--;
			
		}else if(!strcmp(a,"-f:MDPA")){
			if(!features->addFeature(featureType_MDPA,featureMotif_All,0,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:MDP")){
			if(!features->addFeature(featureType_MDP,featureMotif_All,featureMotif_All,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:MDM")){
			if(!features->addFeature(featureType_MDM,featureMotif_All,featureMotif_All,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:MDDA")){
			if(!features->addFeature(featureType_MDDA,featureMotif_All,0,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:MDD")){
			if(!features->addFeature(featureType_MDD,featureMotif_All,featureMotif_All,0,0,0,1,0)){
				return false;
			}
		}else if(!strcmp(a,"-f:GC")){
			if(!features->addFeature(featureType_GC,0,0,0,0,0,1,0)){
				return false;
			}
			
		}else if(!strcmp(a,"-wm:PREdictor")){
			cfg->wmMode=wmPREdictor;
		}else if(!strcmp(a,"-wm:Zero")){
			cfg->wmMode=wmZero;
		}else if(!strcmp(a,"-wm:Constant")){
			cfg->wmMode=wmConstant;
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->loBeta=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-wm:ZZPF")){
			cfg->wmMode=wmZZPF;
		}else if(!strcmp(a,"-wm:PPV")){
			cfg->wmMode=wmPPV;
		}else if(!strcmp(a,"-wm:BiPPV")){
			cfg->wmMode=wmBiPPV;
		}else if(!strcmp(a,"-train:FASTA")){
			if(l>=argc-3){
				argSyntaxError();
				return false;
			}
			if(!trainseq->loadFastaBatch(cargv[1],getSeqClassByName(cargv[2]),getTrainModeByName(cargv[3]))){
				return false;
			}
			if(!registerFile((char*)"Training sequences",cargv[1])){
				return false;
			}
			cargv+=3,argc-=3;
		}else if(!strcmp(a,"-validate:FASTA")){
			if(l>=argc-2){
				argSyntaxError();
				return false;
			}
			if(!valseq->loadFastaBatch(cargv[1],getSeqClassByName(cargv[2]),train_Full)){
				return false;
			}
			if(!registerFile((char*)"Validation sequences",cargv[1])){
				return false;
			}
			cargv+=2,argc-=2;
		}else if(!strcmp(a,"-validate:outSCTable")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->outSCVal=cargv[1];
			cargv++,argc--;
			
		}else if(!strcmp(a,"-validate:no")){
			cfg->validate = false;
			
		}else if(!strcmp(a,"-C:CPREdictor")){
			cfg->classifier=cCPREdictor;
		}else if(!strcmp(a,"-C:DummyPREdictor")){
			cfg->classifier=cDummyPREdictor;
		}else if(!strcmp(a,"-C:SVM")){
			cfg->classifier=cSEQSVM;
			cfg->svmtype=C_SVC;
			
		}else if(!strcmp(a,"-C:SVM-MOCCA")){
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=C_SVC;
		}else if(!strcmp(a,"-C:SVM-MOCCA:C-SVC")){
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=C_SVC;
		}else if(!strcmp(a,"-C:SVM-MOCCA:nu-SVC")){
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=NU_SVC;
		}else if(!strcmp(a,"-C:SVM-MOCCA:One-class")){
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=ONE_CLASS;
			
		}else if(!strcmp(a,"-C:analysis:export")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->CAnalysisExportPath=cargv[1];
			cargv++,argc--;

		}else if(!strcmp(a,"-k:linear")){
			cfg->kernel=kLinear;
		}else if(!strcmp(a,"-k:quadratic")){
			cfg->kernel=kQuadratic;
		}else if(!strcmp(a,"-k:cubic")){
			cfg->kernel=kCubic;
		}else if(!strcmp(a,"-k:RBF")){
			cfg->kernel=kRBF;
		}else if(!strcmp(a,"-in:FASTA")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->inFASTA=cargv[1];
			cargv++,argc--;
		}else if(!strcmp(a,"-out:Wig")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->outWig=cargv[1];
			cargv++,argc--;
		}else if(!strcmp(a,"-out:core-sequence")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->outCoreSequence=cargv[1];
			cargv++,argc--;
		}else if(!strcmp(a,"-SVM:C")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->SVM_C=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-SVM:gamma")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->SVM_gamma=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-SVM:c0")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->SVM_c0=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-SVM:p")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->SVM_p=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-SVM:nu")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->SVM_nu=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-threshold")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->threshold=strtod(cargv[1],0);
			cargv++,argc--;
		}else if(!strcmp(a,"-f:MOCCA:nOcc")){
			cfg->MOCCA_nOcc=true;
		}else if(!strcmp(a,"-f:MOCCA:GC")){
			cfg->MOCCA_GC=true;
		}else if(!strcmp(a,"-f:MOCCA:DNT")){
			cfg->MOCCA_DNT=true;
		}else if(!strcmp(a,"-wSize")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->windowSize=(int)strtol(cargv[1],0,10);
			if(cfg->windowSize<=0){
				argSyntaxError();
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-wStep")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->windowStep=(int)strtol(cargv[1],0,10);
			if(cfg->windowStep<=0){
				argSyntaxError();
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-wStepTrain")){
			if(l>=argc-1){
				argSyntaxError();
				return false;
			}
			cfg->windowStepTrain=(int)strtol(cargv[1],0,10);
			if(cfg->windowStep<=0){
				argSyntaxError();
				return false;
			}
			cargv++,argc--;
		}else if(!strcmp(a,"-no-homo-pairing")){
			cfg->allowHomoPairing = false;
		}else if(!strcmp(a,"-no-hetero-pairing")){
			cfg->allowHeteroPairing = false;
		}else{
			cout << m_error << "Invalid command-line argument \"" << a << "\". Aborting.\n";
			return false;
		}
	}
	return true;
}

/*
libsvm_print_null
	Empty printing function to disable LibSVM output.
*/
void libsvm_print_null(const char*s){}

/*
main
*/
int main(int argc,char**argv){
	timer mainTimer((char*)"Full run");
	cout << sepline << " \033[1;34mMOCCA\033[0m\n Copyright, Bjørn Bredesen, 2013-2019\n Bjorn.Bredesen@ii.uib.no\n" << sepline;
	
	svm_set_print_string_function(&libsvm_print_null);
	
	config*cfg=getConfiguration();
	{
		timeval tv;
		gettimeofday(&tv,0);
		cfg->randSeed=(int)tv.tv_usec;
		srand(cfg->randSeed);
	}
	initIUPACTbl();
	cout.precision(15);
	
	if(argc<=1){
		print_help();
		return -1;
	}
	bool err=false;
	
	// Parse settings
	motifList*motifs=motifList::create();
	seqList*trainseq=seqList::create();
	seqList*valseq=seqList::create();
	featureSet*features=featureSet::create();
	
	if(!parse_arg(argc-1,argv+1,motifs,features,trainseq,valseq)){
		err=true;
	}
	
	// Print settings
	if(!err){
		cfg->printInfo();
		motifs->printInfo();
	}
	if(!err){
		printSeqClasses();
		trainseq->printInfo((char*)"Training sequences");
		printRegisteredFiles();
	}
	
	// Run pipeline
	if(!err)if(!runPipeline(motifs,features,trainseq,valseq))err=true;
	
	// Free memory
	if(motifs)delete motifs;
	if(trainseq)delete trainseq;
	if(valseq)delete valseq;
	return err?-1:0;
}

