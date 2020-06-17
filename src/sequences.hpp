////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Basic sequence streaming

/*
seqStream
	Base-class for sequence streams.
*/
class seqStream{
public:
	/*
	read
		Call to read 'len' bytes to a buffer.
	*/
	virtual int read(int len,char*dest)=0;
	/*
	setpos
		Sets the reading position.
		Does not need to be valid for all kinds of sequence streams.
	*/
	virtual bool setpos(long pos)=0;
	/*
	buffer
		Buffers the whole sequence stream to a buffer (which is
		allocated during buffering).
		Returns the number of bytes buffered.
		
		Note that this method does not need to be valid.
		For some streams, there is no implicit end, and for
		such streams it should never be called.
	*/
	long buffer(char*&dest);
	virtual ~seqStream(){}
};

/*
seqStreamBuffer
	A sequence stream that uses an external buffer.
	The object does not assume ownership of the
	buffer, so allocated memory should be
	freed externally.
*/
class seqStreamBuffer:public seqStream{
private:
	char*buf;		// Buffer pointer
	long bufsize;		// Buffer size
	long cursor;		// Current reading position
	// Private constructor
	seqStreamBuffer();
public:
	/*
	load
		Call to construct.
	*/
	static seqStreamBuffer*create();
	/*
	arm
		Arms the buffer stream with a buffer
	*/
	bool arm(char*b,long bs);
	int read(int len,char*dest);
	bool setpos(long pos);
};

/*
seqStreamRandomIid
	A sequence stream class for i.i.d. randomly generated sequences.
	
	This sequence stream has no end, and buffer()
	should never be called.
*/
class seqStreamRandomIid:public seqStream{
private:
	double rA, rT, rG;
	int nA, nT, nG, nC, nU;
public:
	seqStreamRandomIid();
	bool train(seqStream*input);
	int read(int len,char*dest);
	bool setpos(long pos);
};

typedef struct{
	unsigned int nA, nT, nG, nC, nU, total;
	double rA, rT, rG;
}MCProbability;

/*
seqStreamRandomMarkov
	A sequence stream class for randomly generated sequences by a Markov chain.
	
	This sequence stream has no end, and buffer()
	should never be called.
*/
class seqStreamRandomMC:public seqStream{
private:
	int order;
	int pseudo;
	bool addRC;
	int nprobs;
	int unsigned genstate;
	autofree<MCProbability> probs;
	int nspectrum;
	autofree<int> spectrum;
public:
	seqStreamRandomMC(int _order, int _pseudo = 1, bool _addRC = true);
	bool train(seqStream*input);
	int read(int len,char*dest);
	bool setpos(long pos);
};

////////////////////////////////////////////////////////////////////////////////////
// Sequence file streaming

/*
seqStreamFasta
	Sequence stream class for fasta files.
*/
class seqStreamFasta:public seqStream{
private:
	FILE*f;	// File stream handle.
	// Private constructor.
	seqStreamFasta(FILE*_f);
public:
	/*
	load
		Call to construct.
	*/
	static seqStreamFasta*load(char*path);
	virtual ~seqStreamFasta();
	bool setpos(long pos);
	int read(int len,char*dest);
	long long getLength();
};

/*
seqStreamFastaBatchBlock
	Enables reading a single fasta batch sequence block
	with a sequence stream.
*/
class seqStreamFastaBatchBlock:public seqStream{
private:
	FILE*f;			// File handle to stream from
	char*name;		// Sequence name
	// Private constructor.
	seqStreamFastaBatchBlock(FILE*_f,char*n);
public:
	/*
	load
		Call to construct.
		Assumed to only be called by seqStreamFastaBatch.
	*/
	static seqStreamFastaBatchBlock*load(FILE*_f);
	virtual ~seqStreamFastaBatchBlock();
	int read(int len,char*dest);
	// Reading position setting not implemented. 
	bool setpos(long pos);
	/*
	getName
		Returns the name of the batch sequence
		(NOT a copy, so don't free it).
	*/
	char*getName();
	long long getLength();
};

/*
seqStreamFastaBatch
	Class for reading a fasta sequence batch.
*/
class seqStreamFastaBatch{
private:
	FILE*f;		// File handle to stream from
	// Private constructor.
	seqStreamFastaBatch();
public:
	/*
	load
		Call to construct.
	*/
	static seqStreamFastaBatch*load(char*path);
	virtual ~seqStreamFastaBatch();
	/*
	getBlock
		Call to get the next fasta batch sequence as a stream.
		Returns 0 when there is no next block.
	*/
	seqStreamFastaBatchBlock*getBlock();
};

////////////////////////////////////////////////////////////////////////////////////
// Window streaming

/*
seqStreamWindow
	A class to simplify processing a sequence stream in windows.
*/
class seqStreamWindow{
private:
	char*buf[2];
	bool bufi;
	int winsize,winstep,winkeep;
	seqStream*sstr;
	bool end;
	int winvalid;
	// Private constructor
	seqStreamWindow(seqStream*ss,int wsize,int wstep);
public:
	/*
	create
		Call to construct
	*/
	static seqStreamWindow*create(seqStream*ss,int wsize,int wstep);
	virtual ~seqStreamWindow();
	/*
	get
		Call to get the next window.
		'dest' will be updated with a pointer to the window buffer.
		Returns the size of the window (which may be less than the specified
		size on the end of the sequence stream).
	*/
	int get(char*&dest);
};

/*
createBufferWindow
	Creates a sequence stream and window for a buffer.
	Returns false on failure.
*/
bool createBufferWindow(char*buf,int bufs,seqStreamBuffer*&ssb,seqStreamWindow*&ssw,int wsize,int wstep);

