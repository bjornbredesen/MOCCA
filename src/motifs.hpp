////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Motif list

enum motifType{
	motifType_Invalid,
	motifType_IUPAC,
};

/*
motifListMotif
	Motifs in a motif list.
*/
typedef struct{
	char*name;		// Motif name.
	int len;		// Motif sequence length.
	int index;		// Index of the motif in the list (for quick indexing)
	motifType type;	// Type of the motif.
	void*data;		// Extended data for the motif. Type specific.
	bool skip;		// True of the motif may be skipped.
}motifListMotif;

typedef struct{
	char*seq;
	int nmis;
}IUPACMotif;

/*
motifList
	List of motifs.
*/
class motifList{
private:
	// Private constructor
	motifList();
public:
	~motifList();
	/*
	create
		Call to construct.
	*/
	static motifList*create();
	motifListMotif*motifs;		// Motifs
	int nmotifs;				// Number of motifs
	int maxLen;				// Length of the longest motif
	/*
	addMotif
		Adds a motif to the list.
		Note that this makes copies of the name and sequence.
		Thus, the copy passed to this function can not
		be assumed to be freed with an instance.
	*/
	motifListMotif*addMotif(char*name,motifType type,void*data,int len);
	/*
	addFromXML
		Loads and adds motifs from an XML-file.
	*/
	bool addMotifsFromXML(char*path);
	/*
	addIUPACMotif
		Adds an IUPAC motif
	*/
	motifListMotif*addIUPACMotif(char*name,char*seq,int nmis,bool allowDuplicates=false);
	/*
	addKMers
		Adds k-mers as motifs 
	*/
	bool addKMers(int k,bool allowDuplicates=false);
	/*
	addRandom
		Adds randomly generated motifs
	*/
	bool addRandom(int n,int len);
	/*
	printInfo
		Prints out information
	*/
	void printInfo();
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif occurrence container

/*
motifOcc
	Motif occurrence data structure.
	It is an indexed buffer entry, such that when a new should be created,
	a vacant node may be looked up in constant time from an already
	allocated buffer, or if it is already full, the buffer may be expanded.
	Upon freeing of an entry, it will be placed as a free node.
*/
typedef struct motifOcc{
	long long start;
	motifListMotif*mot;	// Pointer to the motif specification. Should never be 0 for active occurrences.
	//-----------------------------------------------------------------------
	int iANext,iAPrev;	// Index of next/previous active motif. -1 for none.
	int iANextT,iAPrevT;	// Index of next/previous active motif of same type. -1 for none.
	int iFNext;			// Index of next free motif. -1 for none.
	bool active;			// True if active, or false if free.
	bool strand;			// True if reverse complement.
	bool skip;			// True if this occurrence should be skipped.
	//-----------------------------------------------------------------------
	void*extra_buffer;		// Will be freed with the occurrence using free().
}motifOcc;

/*
motifOccContainer
	This class holds motif occurrences in a multiply indexed table, such that
	it may be expanded rarely, and memory may be reused by an
	extra list of free entries.
*/
class motifOccContainer{
private:
	motifOcc*occTbl;			// Table of occurrences.
	int tblSize;				// Size of the table.
	int iAFirst,iALast;			// Indices of first and last active table entries.
	int*iAFirstT,*iALastT;		// Indices of first and last active table entries of specific types.
	int iFree;				// Index of first free table entry.
	motifList*motifs;			// Motifs.
	// Private constructor
	motifOccContainer(motifList*ml);
	/*
	growTable
		Expands the table by the given number of entries.
	*/
	bool growTable(int ngrow);
public:
	motifList*getMotifList();
	// Destruction
	~motifOccContainer();
	/*
	create
		Call to construct
	*/
	static motifOccContainer*create(int basesize,motifList*ml);
	int nOcc;			// Number of occurrences.
	int*nOccT;			// Type-specific numbers of occurrences.
	/*
	createMotifOcc
		Creates a motif occurrence.
		It also expands the table if it has to.
	*/
	motifOcc*createMotifOcc(long long start,motifListMotif*m,bool strand);
	/*
	freeMotifOcc
		Frees the motif occurrence.
		It does not actually free memory, but rather makes the memory
		available in the free table entry list.
	*/
	bool freeMotifOcc(motifOcc*o);
	/*
	flush
		Inactivates all occurrences.
	*/
	void flush();
	// Navigation
	motifOcc*getFirst();
	motifOcc*getLast();
	motifOcc*getNext(motifOcc*o);
	motifOcc*getPrev(motifOcc*o);
	motifOcc*getFirst(int t);
	motifOcc*getLast(int t);
	motifOcc*getNextSame(motifOcc*o);
	motifOcc*getPrevSame(motifOcc*o);
};

////////////////////////////////////////////////////////////////////////////////////
// IUPAC table

void initIUPACTbl();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif window

class motifFSM;

/*
motifWindow
	Parses motif occurrences from sequence windows.
*/
class motifWindow{
private:
	motifFSM*mFSM;
	motifWindow(motifList*_motifs);
	// Naive parsing
	bool initialize();
	/*
	motifMatchIUPAC
		IUPAC motif matching for naive parsing
	*/
	bool motifMatchIUPAC(char*seq,int seqlen,motifListMotif*mot,bool com);
public:
	long long wPos;
	int wLen;
	motifOccContainer*occContainer;
	motifList*motifs;
	
	// Construction
	~motifWindow();
	static motifWindow*create(motifList*_motifs);
	// Processing
	bool flush();
	bool readWindow(char*wseq,long long wpos,int wlen);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif pairing

// Comparison
double getDistance(motifOcc*a,motifOcc*b);

bool overlapping(motifOcc*a,motifOcc*b);

bool isMotifPair(motifOcc*a,motifOcc*b,int cutMin,int cutMax);

