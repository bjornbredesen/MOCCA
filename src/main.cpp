////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
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
#include "models/seqlo.hpp"
#include "models/seqdummy.hpp"

////////////////////////////////////////////////////////////////////////////////////
// File list

/*
fileListFile
	Structure for file list entries.
*/
typedef struct{
	std::string name;
	std::string path;
}fileListFile;

std::vector<fileListFile> fileList;

bool registerFile(std::string name, std::string path){
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
// Arguments

#include <functional>
#include <unordered_map>
#define N_ARGUMENT_PASSES 2

void print_help();
void print_licenses();

struct cmdArg{
	std::string arg;
	int pass;
	int nParam;
	std::string example;
	std::vector<std::string> helpLines;
	std::function<bool (std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq)> apply;
};

cmdArg argumentTypes[] = {
	//---------------------------
	// Base
	{
		// Argument
		"-h",
		// Pass
		0,
		// Parameters
		0,
		// Documentation
		"-h",
		{ "Prints help text." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			print_help();
			return false;
		}
	},
	{
		// Argument
		"-license",
		// Pass
		0,
		// Parameters
		0,
		// Documentation
		"-license",
		{ "Prints license information." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			print_licenses();
			return false;
		}
	},
	{
		// Argument
		"-seed",
		// Pass
		0,
		// Parameters
		1,
		// Documentation
		"-seed VALUE",
		{ "Sets the random seed to VALUE." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->randSeed = (int) strtol(params[0].c_str(), 0, 10);
			srand(cfg->randSeed);
			return true;
		}
	},
	//---------------------------
	// Classes
	{
		// Argument
		"-class",
		// Pass
		0,
		// Parameters
		3,
		// Documentation
		"-class NAME VALUE FLAG",
		{ "Registers a sequence class with name NAME, value/ID VALUE,",
		  "and flag FLAG (\"+\" for positive, or \"-\" for negative)." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			bool cls;
			if(params[2][0] == '+' && !params[2][1]) cls=true;
			else if(params[2][0] == '-' && !params[2][1]) cls=false;
			else{
				argSyntaxError();
				return false;
			}
			if(!registerSeqClass(strtod(params[1].c_str(), 0), (char*)params[0].c_str(), cls)){
				return false;
			}
			return true;
		}
	},
	//---------------------------
	// Motifs
	{
		// Argument
		"-motif:IUPAC",
		// Pass
		1,
		// Parameters
		3,
		// Documentation
		"-motif:IUPAC NAME MOTIF MISMATCHES",
		{ "Adds an IUPAC-motif with the name NAME, and sequence MOTIF,",
		  "with MISMATCHES mismatches allowed." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!ml->addIUPACMotif((char*)params[0].c_str(), (char*)params[1].c_str(), (int)strtol(params[2].c_str(), 0, 10))){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-motif:XML",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-motif:XML PATH",
		{ "Adds motifs specified in an XML-file PATH." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!ml->addMotifsFromXML(params[0])){
				return false;
			}
			if(!registerFile((char*)"Motif XML", params[0])){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-motif:kmer",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-motif:kmer k",
		{ "Adds all k-mers." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!ml->addKMers((int)strtol(params[0].c_str(), 0, 10))){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-motif:Random",
		// Pass
		1,
		// Parameters
		2,
		// Documentation
		"-motif:Random N LEN",
		{ "Adds N random sequences of length LEN." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!ml->addRandom((int)strtol(params[0].c_str(), 0, 10), (int)strtol(params[1].c_str(), 0, 10))){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-motif:FSM",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:FSM",
		{ "Enable Finite State Machine for motif occurrence parsing." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->useFSM=true;
			return true;
		}
	},
	{
		// Argument
		"-motif:No-FSM",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:No-FSM",
		{ "Disable Finite State Machine for motif occurrence parsing." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->useFSM=false;
			return true;
		}
	},
	{
		// Argument
		"-motif:d:centers",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:d:centers",
		{ "Sets the distance mode to centered." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->distanceMode=dmCenters;
			return true;
		}
	},
	{
		// Argument
		"-motif:d:between",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:d:between",
		{ "Sets the distance mode to between." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->distanceMode=dmBetween;
			return true;
		}
	},
	{
		// Argument
		"-motif:d:noOverlap",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:d:noOverlap",
		{ "Disallows motif pairs to overlap when counting." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->motifPairsCanOverlap=false;
			return true;
		}
	},
	{
		// Argument
		"-motif:d:overlap",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-motif:d:overlap",
		{ "Allows motif pairs to overlap when counting." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->motifPairsCanOverlap=true;
			return true;
		}
	},
	{
		// Argument
		"-no-homo-pairing",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-no-homo-pairing",
		{ "Disables pairing of the same motifs for CPREdictor." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->allowHomoPairing = false;
			return true;
		}
	},
	{
		// Argument
		"-no-hetero-pairing",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-no-hetero-pairing",
		{ "Disables pairing of different motifs for CPREdictor." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->allowHeteroPairing = false;
			return true;
		}
	},
	//---------------------------
	// General features
	{
		// Argument
		"-f:nOcc",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:nOcc",
		{ "Adds motif occurrence frequency features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_nOcc, featureMotif_All, 0, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:nPair",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-f:nPair D",
		{ "Adds motif pair occurrence frequency features, with distance",
		  "cutoff D." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_nPair, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[0].c_str(), 0),0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:nOccPair",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-f:nOccPair D",
		{ "Adds motif occurrence pair occurrence frequency features, with",
		  "distance cutoff D." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_nOccPair, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[0].c_str(), 0), 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:MDPA",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MDPA",
		{ "Adds Mean Distance Proximal All features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_MDPA, featureMotif_All, 0, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:MDP",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MDP",
		{ "Adds Mean Distance Proximal features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_MDP, featureMotif_All, featureMotif_All, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:MDM",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MDM",
		{ "Adds Mean Distance Mean features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_MDM, featureMotif_All, featureMotif_All, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:MDDA",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MDDA",
		{ "Adds Mean Distance Distal All features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_MDDA, featureMotif_All, 0, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:MDD",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MDD",
		{ "Adds Mean Distance Distal features." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_MDD, featureMotif_All, 0, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:GC",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:GC",
		{ "Adds GC-content feature." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_GC, 0, 0, 0, 0, 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:nPair2D",
		// Pass
		1,
		// Parameters
		2,
		// Documentation
		"-f:nPair2D D F",
		{ "Adds 2D oscillatory motif pair feature.",
		  "Pair occurrence within D base pairs is weighted by the sine",
		  "and cosine of the distance separately, with the periodicity",
		  "scaled to F." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_nPair2D, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[0].c_str(), 0), strtod(params[1].c_str(), 0), 1, 0)){
				return false;
			}
			if(!features->addFeature(featureType_nPair2D, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[0].c_str(), 0), strtod(params[1].c_str(), 0), 1, 1)){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:PEDI",
		// Pass
		1,
		// Parameters
		3,
		// Documentation
		"-f:PEDI D F P",
		{ "Adds oscillatory motif pair feature.",
		  "Pair occurrence within D base pairs is weighted by the cosine",
		  "of the distance, with the periodicity scaled to F, and phase",
		  "shifted by P." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_PEDI, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[2].c_str(), 0), strtod(params[1].c_str(), 0), 1, strtod(params[0].c_str(), 0))){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-f:nPairDH",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-f:nPairDH D",
		{ "Adds B-DNA oscillatory motif pair feature.",
		  "Pairs within D base pairs are weighted by their distance with",
		  "a cosine curve with a frequency of 10.5 base pairs, with the",
		  "curve shifted by 5.25 base pairs for occurrences on opposite",
		  "strands." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!features->addFeature(featureType_nPairDH, featureMotif_All, featureMotif_All, featureMotif_All, strtod(params[0].c_str(), 0), 0, 1, 0)){
				return false;
			}
			return true;
		}
	},
	//---------------------------
	// SVM-MOCCA features
	{
		// Argument
		"-f:MOCCA:nOcc",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MOCCA:nOcc",
		{ "Adds motif occurrence frequency features to SVM-MOCCA." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->MOCCA_nOcc=true;
			return true;
		}
	},
	{
		// Argument
		"-f:MOCCA:DNT",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MOCCA:DNT",
		{ "Adds dinucleotide features to SVM-MOCCA." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->MOCCA_DNT=true;
			return true;
		}
	},
	{
		// Argument
		"-f:MOCCA:GC",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-f:MOCCA:GC",
		{ "Adds GC content feature to SVM-MOCCA." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->MOCCA_GC=true;
			return true;
		}
	},
	//---------------------------
	// Classifiers
	{
		// Argument
		"-C:CPREdictor",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:CPREdictor",
		{ "Sets the classifier to the CPREdictor." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cCPREdictor;
			return true;
		}
	},
	{
		// Argument
		"-C:DummyPREdictor",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:DummyPREdictor",
		{ "Sets the classifier to the DummyPREdictor." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cDummyPREdictor;
			return true;
		}
	},
	{
		// Argument
		"-C:LO",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:LO",
		{ "Sets the classifier to log-odds, using separately specified",
		  "feature spaces." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSEQLO;
			return true;
		}
	},
	{
		// Argument
		"-C:Dummy",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:Dummy",
		{ "Sets the classifier to dummy (un-weighted sum), using separately",
		  "specified, feature spaces." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSEQDummy;
			return true;
		}
	},
	{
		// Argument
		"-C:SVM",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:SVM",
		{ "Sets the classifier to SVM, using separately specified",
		  "feature spaces." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSEQSVM;
			cfg->svmtype=C_SVC;
			return true;
		}
	},
	{
		// Argument
		"-C:SVM-MOCCA",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:SVM-MOCCA",
		{ "Sets the classifier to SVM-MOCCA." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=C_SVC;
			return true;
		}
	},
	{
		// Argument
		"-C:SVM-MOCCA:C-SVC",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:SVM-MOCCA:C-SVC",
		{ "Sets the classifier to SVM-MOCCA, using the C-SVC formulation." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=C_SVC;
			return true;
		}
	},
	{
		// Argument
		"-C:SVM-MOCCA:nu-SVC",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-C:SVM-MOCCA:nu-SVC",
		{ "Sets the classifier to SVM-MOCCA, using the nu-SVC formulation." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->classifier=cSVMMOCCA;
			cfg->svmtype=NU_SVC;
			return true;
		}
	},
	//---------------------------
	// Classifier configuration
	{
		// Argument
		"-C:analysis:export",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-C:analysis:export PATH",
		{ "Exports model analysis to file." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->CAnalysisExportPath = params[0];
			return true;
		}
	},
	{
		// Argument
		"-wm:PREdictor",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-wm:PREdictor",
		{ "Sets the log-odds weight mode to the PREdictor formulation.",
		  "When either a positive or negative class summed feature is zero,",
		  "a pseudocount of 1 is added for both the positives and negatives." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->wmMode=wmPREdictor;
			return true;
		}
	},
	{
		// Argument
		"-wm:Zero",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-wm:Zero",
		{ "Sets the log-odds weight mode to zero.",
		  "When either a positive or negative class summed feature is zero,",
		  "the weight is set to zero." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->wmMode=wmZero;
			return true;
		}
	},
	{
		// Argument
		"-wm:Constant",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-wm:Constant VALUE",
		{ "Sets the log-odds weight mode to constant.",
		  "A constant pseudocount is added to summed positive and negative",
		  "features, so that they are never zero." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->loBeta=strtod(params[0].c_str(),0);
			cfg->wmMode=wmConstant;
			return true;
		}
	},
	{
		// Argument
		"-wm:PPV",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-wm:PPV",
		{ "Sets the log-odds weight mode to PPV.",
		  "Instead of log-odds, the Positive Predictive Value is calculated",
		  "based on motif occurrences in positives (TP) versus negatives (FP)." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->wmMode=wmPPV;
			return true;
		}
	},
	{
		// Argument
		"-wm:BiPPV",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-wm:BiPPV",
		{ "Sets the log-odds weight mode to BiPPV.",
		  "Instead of log-odds, the Bi-Positive Predictive Value is calculated",
		  "based on motif occurrences in positives (TP) versus negatives (FP).",
		  "BiPPV is the difference between PPV for the positives, and PPV",
		  "for switched labels (False Discovery Rate)." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->wmMode=wmBiPPV;
			return true;
		}
	},
	{
		// Argument
		"-k:linear",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-k:linear",
		{ "Sets the SVM kernel to linear." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->kernel=kLinear;
			return true;
		}
	},
	{
		// Argument
		"-k:quadratic",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-k:quadratic",
		{ "Sets the SVM kernel to quadratic (poly. degree 2)." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->kernel=kQuadratic;
			return true;
		}
	},
	{
		// Argument
		"-k:cubic",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-k:cubic",
		{ "Sets the SVM kernel to cubic (poly. degree 3)." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->kernel=kCubic;
			return true;
		}
	},
	{
		// Argument
		"-k:RBF",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-k:RBF",
		{ "Sets the SVM kernel to Radial Basis Function." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->kernel=kRBF;
			return true;
		}
	},
	{
		// Argument
		"-SVM:C",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-SVM:C VALUE",
		{ "Sets the SVM metaparameter C." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->SVM_C=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-SVM:nu",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-SVM:nu VALUE",
		{ "Sets the SVM metaparameter nu." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->SVM_nu=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-SVM:p",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-SVM:p VALUE",
		{ "Sets the SVM metaparameter p." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->SVM_p=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-SVM:gamma",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-SVM:gamma VALUE",
		{ "Sets the SVM kernel parameter gamma." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->SVM_gamma=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-SVM:c0",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-SVM:c0 VALUE",
		{ "Sets the SVM kernel parameter c0." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->SVM_c0=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-threshold",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-threshold VALUE",
		{ "Sets the classifier threshold." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->threshold=strtod(params[0].c_str(), 0);
			return true;
		}
	},
	{
		// Argument
		"-wSize",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-wSize VALUE",
		{ "Sets the classifier window size." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->windowSize = (int)strtol(params[0].c_str(), 0, 10);
			if(cfg->windowSize <= 0){
				argSyntaxError();
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-wStep",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-wStep VALUE",
		{ "Sets the classifier window step size." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->windowStep = (int)strtol(params[0].c_str(), 0, 10);
			if(cfg->windowStep <= 0){
				argSyntaxError();
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-wStepTrain",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-wStepTrain VALUE",
		{ "Sets the classifier window training step size.",
		  "When training with sequence windows, this can be set in order to",
		  "control model complexity growth." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->windowStepTrain = (int)strtol(params[0].c_str(), 0, 10);
			if(cfg->windowStepTrain <= 0){
				argSyntaxError();
				return false;
			}
			return true;
		}
	},
	//---------------------------
	// Sequences
	{
		// Argument
		"-train:FASTA",
		// Pass
		1,
		// Parameters
		3,
		// Documentation
		"-train:FASTA PATH CLASS MODE",
		{ "Adds a training sequence file.",
		  "PATH: Path to FASTA file.",
		  "CLASS: A class ID, defined with \"-class\", or one of the",
		  "pre-specified binary classes: \"+\" for positive or \"-\"",
		  "for negative.",
		  "MODE: Can be \"win\", for training with all windows within",
		  "each training sequence file, or \"full\", for training with",
		  "the full sequences." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!trainseq->loadFastaBatch((char*)params[0].c_str(), getSeqClassByName(params[1]), getTrainModeByName((char*)params[2].c_str()))){
				return false;
			}
			if(!registerFile("Training sequences",params[0])){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-validate:FASTA",
		// Pass
		1,
		// Parameters
		2,
		// Documentation
		"-validate:FASTA PATH CLASS",
		{ "Adds a validation sequence file.",
		  "PATH: Path to FASTA file.",
		  "CLASS: A class ID, defined with \"-class\", or one of the",
		  "pre-specified binary classes: \"+\" for positive or \"-\"",
		  "for negative." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			if(!valseq->loadFastaBatch((char*)params[0].c_str(), getSeqClassByName(params[1]), train_Full)){
				return false;
			}
			if(!registerFile("Validation sequences",params[0])){
				return false;
			}
			return true;
		}
	},
	{
		// Argument
		"-validate:no",
		// Pass
		1,
		// Parameters
		0,
		// Documentation
		"-validate:no",
		{ "Disables validation." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->validate = false;
			return true;
		}
	},
	{
		// Argument
		"-validate:outSCTable",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-validate:outSCTable PATH",
		{ "Outputs validation set score and class table to file." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->outSCVal = params[0];
			return true;
		}
	},
	{
		// Argument
		"-in:FASTA",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-in:FASTA PATH",
		{ "Sets an input FASTA file to be scored." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->inFASTA = params[0];
			return true;
		}
	},
	{
		// Argument
		"-out:Wig",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-out:Wig PATH",
		{ "Sets an output Wiggle file for scored sequences." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->outWig = params[0];
			return true;
		}
	},
	{
		// Argument
		"-out:core-sequence",
		// Pass
		1,
		// Parameters
		1,
		// Documentation
		"-out:core-sequence PATH",
		{ "Sets an output FASTA file for predicted core sequences",
		  "from input FASTA file." },
		// Code
		[](std::vector<std::string> params, config*cfg, motifList*ml, featureSet*features, seqList*trainseq, seqList*valseq) -> bool {
			cfg->outCoreSequence = params[0];
			return true;
		}
	},
};

/*
print_help
	Outputs help message.
*/
void print_help(){
	cout << " Usage:\n";
	char argFmtStr[] = " %35s - %s\n";
	char argFmtStrC[] = " %35s   %s\n";
	for(auto argType: argumentTypes) {
		printf(argFmtStr, argType.example.c_str(), argType.helpLines[0].c_str());
		for(unsigned int i = 1; i < argType.helpLines.size(); i++) {
			printf(argFmtStrC, "", argType.helpLines[i].c_str());
		}
	}
	cout << sepline;
}

/*
print_help
	Outputs license information.
*/
void print_licenses(){
	cout << " Licenses:\n";
	cout << sepline;
	cout << t_bold("SVM-MOCCA") << "\n\n";
	cout << "MIT License\n\
\n\
Copyright (c) 2019 Bjørn André Bredesen\n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy\n\
of this software and associated documentation files (the \"Software\"), to deal\n\
in the Software without restriction, including without limitation the rights\n\
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n\
copies of the Software, and to permit persons to whom the Software is\n\
furnished to do so, subject to the following conditions:\n\
\n\
The above copyright notice and this permission notice shall be included in all\n\
copies or substantial portions of the Software.\n\
\n\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n\
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n\
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n\
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n\
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n\
SOFTWARE.\n";
	cout << sepline;
	cout << t_bold("Dependency license: LibSVM") << "\n\
\n\
Copyright (c) 2000-2013 Chih-Chung Chang and Chih-Jen Lin\n\
All rights reserved.\n\
\n\
Redistribution and use in source and binary forms, with or without\n\
modification, are permitted provided that the following conditions\n\
are met:\n\
\n\
1. Redistributions of source code must retain the above copyright\n\
notice, this list of conditions and the following disclaimer.\n\
\n\
2. Redistributions in binary form must reproduce the above copyright\n\
notice, this list of conditions and the following disclaimer in the\n\
documentation and/or other materials provided with the distribution.\n\
\n\
3. Neither name of copyright holders nor the names of its contributors\n\
may be used to endorse or promote products derived from this software\n\
without specific prior written permission.\n\
\n\
\n\
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n\
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n\
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n\
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR\n\
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,\n\
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,\n\
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n\
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n\
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n\
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n\
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";
	cout << sepline;
	cout << t_bold("Dependency license: RapidXML") << "\n\
\n\
Copyright (c) 2006, 2007 Marcin Kalicinski\n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy \n\
of this software and associated documentation files (the \"Software\"), to deal \n\
in the Software without restriction, including without limitation the rights \n\
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies \n\
of the Software, and to permit persons to whom the Software is furnished to do so, \n\
subject to the following conditions:\n\
\n\
The above copyright notice and this permission notice shall be included in all \n\
copies or substantial portions of the Software.\n\
\n\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR \n\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, \n\
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL \n\
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER \n\
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, \n\
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS \n\
IN THE SOFTWARE.\n";
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
	std::unordered_map<std::string, cmdArg> argumentMap;
	for(auto argType: argumentTypes) {
		argumentMap[argType.arg] = argType;
	}
	for(int pass = 0; pass < N_ARGUMENT_PASSES; pass ++) {
		char**cargv = argv;
		int argc = _argc;
		for(int l=0; l<argc; l++, cargv++){
			if (!(*cargv)) continue;
			std::string cArg = std::string(*cargv);
			auto cArgTypeF = argumentMap.find(cArg);
			// Ensure argument type is familiar
			if (cArgTypeF == argumentMap.end()) {
				cout << m_error << "Invalid command-line argument \"" << cArg << "\". Aborting.\n";
				return false;
			}
			// Skip if for different pass
			auto cArgType = cArgTypeF->second;
			// Ensure we have enough arguments
			if(l >= argc-cArgType.nParam){
				argSyntaxError();
				return false;
			}
			// Accumulate arguments to list
			char**ccargv = cargv+1;
			std::vector<std::string> params;
			for(int y=0; y<cArgType.nParam; y++, ccargv++){
				params.push_back(std::string(*ccargv));
			}
			if (cArgType.pass == pass) {
				// Call argument processor
				if(!cArgType.apply(params, cfg, ml, features, trainseq, valseq))
					return false;
				// Remove argument (null-arguments are skipped above)
				ccargv = cargv;
				for(int y=0; y<=cArgType.nParam; y++, ccargv++)
					*ccargv = 0;
			}
			cargv += cArgType.nParam, argc -= cArgType.nParam;
		}
		// Add standard sequence classes if none have been explicitly registered
		if(pass == 0 && !nSequenceClasses()){
			if(!registerSeqClass(1,(char*)"Positive",true)){
				return false;
			}
			if(!registerSeqClass(-1,(char*)"Negative",false)){
				return false;
			}
		}
	}
	return true;
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
		case cSEQLO:cls=SEQLO::create(motifs,features);break;
		case cSEQDummy:cls=SEQDummy::create(motifs,features);break;
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
		if(cfg->classifier == cSEQSVM || cfg->classifier == cSEQLO || cfg->classifier == cSEQDummy)
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
			if(cfg->outSCVal.length() > 0)if(!saveVPairTable(cfg->outSCVal, vp.ptr, nvp))return false;
		}
		if(cfg->inFASTA.length() > 0 && cfg->outWig.length() > 0){
			cmdSection("FASTA scoring");
			cls.ptr->applyFASTA(cfg->inFASTA, cfg->outWig);
		}
		if(cfg->inFASTA.length() > 0 && cfg->outCoreSequence.length() > 0){
			cmdSection("FASTA scoring");
			cls.ptr->predictCoreSequence(cfg->inFASTA,cfg->outCoreSequence);
		}
		if(cfg->CAnalysisExportPath.length() > 0){
			cmdSection("Classifier analysis export");
			cls.ptr->exportAnalysisData(cfg->CAnalysisExportPath);
		}
	}
	cout << sepline;
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
	cout << sepline << " \033[1;34mMOCCA\033[0m\n Copyright, Bjørn Bredesen, 2013-2019\n bjorn@bjornbredesen.no\n" << sepline;
	
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
		print_licenses();
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

