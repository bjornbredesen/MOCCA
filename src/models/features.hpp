////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Feature set

/*
featureMotif_All
	Used to indicate that a motif index corresponds to all motifs.
*/
#define featureMotif_All -1

#define N_FEATURES 17

enum featureType{
	featureType_Invalid,
	featureType_nOcc,
	featureType_nPair,
	featureType_MDPA,
	featureType_MDP,
	featureType_MDM,
	featureType_MDDA,
	featureType_MDD,
	featureType_PEDI,
	featureType_nPairDH,
	featureType_nPairCosI,
	featureType_nOccPair,
	featureType_GC,
	featureType_nPair2D,
	featureType_nPairRgn,
};

/*
featureSetFeature
	Feature set features.
*/
typedef struct{
	featureType f;		// Type of feature
	int ia,ib,ic;		// Motif indexing parameters
	int icf;			// Cardinality factor index
	int cf;				// Cardinality factor
	double da,db;		// Floating point parameters
}featureSetFeature;

/*
featureSetInstFeature
	Instantiated features (including calibration and feature value)
*/
typedef struct{
	featureSetFeature*fsf;	// Base featureSetFeature
	int ia,ib,ic;			// Instantiated motif indexing parameters
	int icf;				// Cardinality factor index
}featureSetInstFeature;

/*
featureSet
	Holds a feature set.
*/
class featureSet{
private:
	// Private constructor
	featureSet();
public:
	autofree<featureSetFeature> features;		// Features
	int nfeatures;								// Number of features
	autofree<featureSetInstFeature> ifeatures;	// Instantiated features
	int nifeatures;								// Number of instantiated features
	bool instantiated;							// True if instantiated
	~featureSet();
	/*
	create
		Call to construct
	*/
	static featureSet*create();
	/*
	addFeature
		Adds a feature to the feature set.
		ia, ib, ic, da, db, cf and icf define parameters for the features.
	*/
	featureSetFeature*addFeature(featureType _f,int ia,int ib,int ic,double da,double db,int cf,int icf);
	/*
	uninstantiate
		Uninstantiates features.
	*/
	void uninstantiate();
	/*
	addInstFeature
		Adds and returns an instantiated feature.
		Note that it takes motif specific indices, and otherwise uses
		the settings from the feature that is instantiated.
	*/
	featureSetInstFeature*addInstFeature(featureSetFeature*bf,int ia,int ib,int ic,int icf);
	/*
	instantiateFeature
		Instantiates a feature for the specified motifList.
	*/
	bool instantiateFeature(featureSetFeature*f,motifList*ml);
	/*
	instantiate
		Instantiates the feature set for the specified motifList.
	*/
	bool instantiate(motifList*ml);
	bool skipUnusedMotifs(motifList*ml);
	/*
	printInfo
		Prints out information
	*/
	void printInfo();
	/*
	printInfoI
		Prints out instantiated feature information
	*/
	void printInfoI(motifList*ml);
	/*
	printInstFeatureName
		Prints out the name of an instantiated feature
	*/
	void printInstFeatureName(featureSetInstFeature*f,motifList*ml);
	/*
	*/
	std::vector<std::string> getInstFeatureNames(motifList*ml);
};


////////////////////////////////////////////////////////////////////////////////////
// Feature window

/*
featureWindow
	Used to extract features from a sequence based on motifs.
	Constructed for a particular combination of features and motifs.
*/
class featureWindow{
private:
	config*cfg;
	featureSet*features;				// Features used by the window
	motifList*motifs;					// Motifs used by the window
	motifOccContainer*occContainer;		// Motif occurrence container for the window
	motifWindow*mwin;
	double*fvec;
	
	// Private constructor
	featureWindow(motifWindow*mw,featureSet*fs);
public:
	~featureWindow();
	/*
	create
		Call to construct
		Note that featureWindow does not assume ownership of the
		memory for motifList and featureSet, so freeing should
		be done externally.
	*/
	static featureWindow*create(motifWindow*mw,featureSet*fs);
	/*
	getMDP
		Gets Mean Distance Proximal.
	*/
	double getMDP(int a,int b);
	/*
	getMDPA
		Gets Mean Distance Proximal, Any class for the given motif type.
	*/
	double getMDPA(int a);
	/*
	getMDM
		Gets Mean Distance Mean.
	*/
	double getMDM(int a,int b);
	/*
	getMDDA
		Gets Mean Distance Distal, Any class for the given motif type.
	*/
	double getMDDA(int a);
	/*
	getMDD
		Gets Mean Distance Distal.
	*/
	double getMDD(int a,int b);
	/*
	getNFeatures
		Returns the number of instantiated features.
	*/
	int getNFeatures();
	/*
	extractFeatures
		Extracts features based on a parsed window.
	*/
	double*extractFeatures(char*wseq,long long wpos,int wlen,bool doRead);
};


