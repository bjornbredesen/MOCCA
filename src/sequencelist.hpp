////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2017
// E-mail: Bjorn.Bredesen@ii.uib.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Sequence classes

/*
seqClass
	Classification class specification.
*/
typedef struct{
	double cls;		// Class identifier value.
	char*name;		// Name of the class.
	bool flag;		// True if this class should be treated as true.
}seqClass;

seqClass*getSeqClassByName(char*name);

/*
printSeqClasses
	Outputs a list of classes.
*/
void printSeqClasses();

/*
registerSeqClass
	Registers a classification class.
	Note: Should only be called in the first parameter parsing pass.
		Otherwise, pointers may be in use, and memory reallocation
		may invalidate these.
*/
seqClass*registerSeqClass(double cls,char*name,bool flag);

size_t nSequenceClasses();

////////////////////////////////////////////////////////////////////////////////////
// Sequence list

enum e_trainMode{
	train_Invalid,
	train_Full,
	train_Windows,
};

e_trainMode getTrainModeByName(char*name);

/*
seqListSeq
	Sequences in a sequence list.
*/
typedef struct{
	char*name;				// Sequence name
	char*buf;				// Sequence buffer
	int bufs;				// Sequence buffer length
	e_trainMode trainMode;	// Training mode for the sequence
	seqClass*cls;		// Corresponding class specification
}seqListSeq;

/*
seqList
	A datastructure to hold a collection of sequences.
	May be used both for training sequences and potentially
	other validation sequences.
*/
class seqList{
private:
	// Private constructor
	seqList();
	/*
	addSeq
		Call to add a sequence.
		Note that the structure assumes ownership of
		the buffer if successful, so as to not make
		unnecessary copies.
		Therefore, free the buffer manually on failure,
		and otherwise don't free it outside.
	*/
	seqListSeq*addSeq(char*name,char*buf,int bufs,seqClass*cls,e_trainMode tm);
public:
	int nseq;			// Number of sequences.
	seqListSeq*seq;		// Sequences.
	int nNeg,nPos;		// Numbers of negative and positive sequences.
	bool ownSeq;		// True if the sequence list owns the sequence memory.
	~seqList();
	/*
	create
		Call to construct.
	*/
	static seqList*create();
	/*
	create
		Call to construct with a preallocated, empty sequence table.
	*/
	static seqList*create(int n);
	/*
	loadFastaBatch
		Loads sequences from a fasta batch with the classification flag 'cls'.
	*/
	bool loadFastaBatch(char*path,seqClass*cls,e_trainMode tm);
	/*
	printInfo
		Prints out information
	*/
	void printInfo(char*title);
};


