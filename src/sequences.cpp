////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "aux.hpp"
#include "sequences.hpp"

////////////////////////////////////////////////////////////////////////////////////
// seqStream

long seqStream::buffer(char*&dest){
	long ind=0;
	int al=128;
	dest=(char*)malloc(sizeof(char)*al);
	if(!dest){
		outOfMemory();
		return 0;
	}
	int nr;
	for(;;){
		if(ind>al-128){
			al<<=1;
			dest=(char*)realloc(dest,sizeof(char)*al);
			if(!dest){
				outOfMemory();
				return 0;
			}
		}
		nr=read(128,&dest[ind]);
		ind+=nr;
		if(!nr){
			return ind;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamBuffer

seqStreamBuffer::seqStreamBuffer(){
	cursor=0;
	buf=0;
	bufsize=0;
}

seqStreamBuffer*seqStreamBuffer::create(){
	seqStreamBuffer*r=new seqStreamBuffer();
	if(!r){
		outOfMemory();
		return 0;
	}
	return r;
}

bool seqStreamBuffer::arm(char*b,long bs){
	cursor=0;
	buf=0;
	bufsize=0;
	if(!b||!bs){
		cmdError("Invalid buffer.");
		return false;
	}
	buf=b;
	bufsize=bs;
	return true;
}

int seqStreamBuffer::read(int len,char*dest){
	if(!buf||!bufsize)return 0;
	len=min(len,int(bufsize-cursor));
	if(len<=0)return 0;
	memcpy(dest,&buf[cursor],sizeof(char)*len);
	cursor+=len;
	return len;
}

bool seqStreamBuffer::setpos(long pos){
	if(pos<0||pos>=bufsize)return false;
	cursor=pos;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamRandom

seqStreamRandom::seqStreamRandom(double da,double dc,double dg){ // dt = 1-da-dc-dg.
	rA = da;
	rC1 = da+dc;
	rC2 = rC1+dg;
}

bool seqStreamRandom::train(seqStream*input){
	cmdTask task((char*)"Training i.i.d. background model");
	#define FGCSBUFSIZE 256
	autofree<char> buf((char*)malloc(FGCSBUFSIZE));
	// Extract occurrence frequencies per nucleotide
	int nNT[5]={0,0,0,0,0};
	long nti=0;
	for(long wind=0;;wind++){
		int nread=input->read(FGCSBUFSIZE,buf.ptr);
		if(!nread)break;
		char*c=buf.ptr;
		for(int y=0;y<nread;y++,c++){
			if(*c=='A')nNT[0]++;
			else if(*c=='C')nNT[1]++;
			else if(*c=='G')nNT[2]++;
			else if(*c=='T')nNT[3]++;
			else nNT[4]++;
		}
		nti+=nread;
		if(!(wind%100))task.setLongT(nti,(char*)"nt");
	}
	// Calculate weights
	int bptotal = nNT[0]+nNT[1]+nNT[2]+nNT[3];
	int bpinvalid = nNT[4];
	rA = double(nNT[0])/double(bptotal);
	rC1 = double(nNT[1])/double(bptotal);
	rC2 = double(nNT[2])/double(bptotal);
	return true;
}

int seqStreamRandom::read(int len,char*dest){
	double frv;
	char*d=dest;
	for(int x=0;x<len;x++,d++){
		frv=double(rand())/double(RAND_MAX);
		if(frv<rA)(*d)='A';
		else if(frv<rC1)(*d)='C';
		else if(frv<rC2)(*d)='G';
		else (*d)='T';
	}
	return len;
}

bool seqStreamRandom::setpos(long pos){
	return false;
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamFasta

seqStreamFasta::seqStreamFasta(FILE*_f){
	f=_f;
}

seqStreamFasta*seqStreamFasta::load(char*path){
	FILE*_f=fopen(path,"rb");
	if(!_f){
		cmdError("Could not open fasta file.");
		cout << "    (\"" << path << "\")\n\n";
		return 0;
	}
	seqStreamFasta*r=new seqStreamFasta(_f);
	if(!r){
		outOfMemory();
		fclose(_f);
		return 0;
	}
	return r;
}

seqStreamFasta::~seqStreamFasta(){
	if(f)fclose(f);
}

bool seqStreamFasta::setpos(long pos){
	return false;
}

int seqStreamFasta::read(int len,char*dest){
	char*b=dest,c;
	int r=0;
	for(;r<len;){
		c=(char)fgetc(f);
		if(c==EOF){
			*b=0;
			return r;
		}else if(c=='>'){
			for(;;){
				c=(char)fgetc(f);
				if(c==EOF){
					*b=0;
					return r;
				}else if(c=='\n'){
					break;
				}
			}
		}else if(c!='\n'){
			*(b++)=c;
			r++;
		}
	}
	return r;
}

long long seqStreamFasta::getLength(){
	long long bpt=0;
	#define FGCSBUFSIZE 256
	char buf[FGCSBUFSIZE];
	cmdTask taskp((char*)"Processed");
	taskp.setLongT(0,(char*)"nt");
	cmdTask::refresh();
	for(long wind=0;;wind++){
		int nread=read(FGCSBUFSIZE,buf);
		if(!nread)break;
		bpt+=(long long)nread;
		if(!(wind%100))taskp.setLongT((long int)bpt,(char*)"nt");
	}
	return bpt;
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamFastaBatchBlock

seqStreamFastaBatchBlock::seqStreamFastaBatchBlock(FILE*_f,char*n){
	name=n;
	f=_f;
}

seqStreamFastaBatchBlock*seqStreamFastaBatchBlock::load(FILE*_f){
	if(!_f)return 0;
	char c;
	for(;;){
		c=(char)fgetc(_f);
		if(c=='>'||feof(_f)){
			break;
		}
	}
	if(feof(_f)){
		return 0;
	}
	// Read name
	int bufS=16;
	char*buf=(char*)malloc(sizeof(char)*bufS);
	if(!buf){
		outOfMemory();
		return 0;
	}
	char*b=buf;
	for(;;){
		if(b+1>=buf+bufS){
			int ibuf=int(b-buf);
			bufS<<=1;
			buf=(char*)realloc(buf,sizeof(char)*bufS);
			if(!buf){
				outOfMemory();
				return 0;
			}
			b=&buf[ibuf];
		}
		c=(char)fgetc(_f);
		if(c=='\n'||c==EOF){
			*b=0;
			break;
		}
		*(b++)=c;
	}
	seqStreamFastaBatchBlock*r=new seqStreamFastaBatchBlock(_f,buf);
	if(!r){
		outOfMemory();
		free(buf);
		return 0;
	}
	return r;
}

seqStreamFastaBatchBlock::~seqStreamFastaBatchBlock(){
	if(name)free(name);
}

int seqStreamFastaBatchBlock::read(int len,char*dest){
	char*b=dest,c;
	int r=0;
	for(;r<len;){
		c=(char)fgetc(f);
		if(c=='>'||c==EOF){
			*b=0;
			if(c=='>'){
				fseek(f,-1,SEEK_CUR);
			}
			break;
		}else if(c!='\n'){
			*(b++)=c;
			r++;
		}
	}
	return r;
}

bool seqStreamFastaBatchBlock::setpos(long pos){
	return false;
}

char*seqStreamFastaBatchBlock::getName(){
	return name;
}

long long seqStreamFastaBatchBlock::getLength(){
	long long bpt=0;
	#define FGCSBUFSIZE 256
	char buf[FGCSBUFSIZE];

	cmdTask taskp((char*)"Processed");
	taskp.setLongT(0,(char*)"nt");
	cmdTask::refresh();

	for(long wind=0;;wind++){
		int nread=read(FGCSBUFSIZE,buf);
		if(!nread)break;
		bpt+=(long long)nread;
		if(!(wind%100))taskp.setLongT((long int)bpt,(char*)"nt");
	}
	return bpt;
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamFastaBatch

seqStreamFastaBatch::seqStreamFastaBatch(){
	f=0;
}

seqStreamFastaBatch*seqStreamFastaBatch::load(char*path){
	seqStreamFastaBatch*r=new seqStreamFastaBatch();
	if(!r){
		outOfMemory();
		return 0;
	}
	r->f=fopen(path,"rb");
	if(!r->f){
		cout << m_error << "Could not open file \"" << path << "\" for reading.\n";
		delete r;
		return 0;
	}
	return r;
}

seqStreamFastaBatch::~seqStreamFastaBatch(){
	if(f)fclose(f);
}

seqStreamFastaBatchBlock*seqStreamFastaBatch::getBlock(){
	if(!f||feof(f)){
		return 0;
	}
	seqStreamFastaBatchBlock*fbb=seqStreamFastaBatchBlock::load(f);
	if(!fbb)return 0;
	if(feof(f)){
		delete fbb;
		return 0;
	}
	return fbb;
}

////////////////////////////////////////////////////////////////////////////////////
// seqStreamWindow

seqStreamWindow::seqStreamWindow(seqStream*ss,int wsize,int wstep){
	sstr=ss;
	winsize=wsize;
	winstep=wstep;
	winkeep=max(wsize-wstep,0);
	buf[0]=buf[1]=0;
	bufi=false;
	end=false;
	winvalid=0;
}

seqStreamWindow*seqStreamWindow::create(seqStream*ss,int wsize,int wstep){
	if(!ss||wstep>=wsize||wsize<=0)return 0;
	seqStreamWindow*r=new seqStreamWindow(ss,wsize,wstep);
	if(!r){
		outOfMemory();
		return 0;
	}
	
	r->buf[0]=(char*)malloc(sizeof(char)*wsize);
	if(!r->buf[0]){
		outOfMemory();
		delete r;
		return 0;
	}
	r->buf[1]=(char*)malloc(sizeof(char)*wsize);
	if(!r->buf[1]){
		outOfMemory();
		delete r;
		return 0;
	}
	
	char*pb=r->buf[1];
	r->winvalid=ss->read(r->winkeep,&pb[wstep]);
	return r;
}

seqStreamWindow::~seqStreamWindow(){
	if(buf[0])free(buf[0]);
	if(buf[1])free(buf[1]);
}

int seqStreamWindow::get(char*&dest){
	if(end)return 0;
	// If it read less than what it should in the beginning, just return it and mark this point as the end.
	if(winvalid<winkeep){
		char*b=buf[1];
		dest=&b[winstep];
		end=true;
		return winvalid;
	}
	// Otherwise, try to read more.
	char*b=buf[bufi?1:0];
	char*pb=buf[bufi?0:1];
	dest=b;
	bufi=!bufi;
	memcpy(b,&pb[winstep],sizeof(char)*winkeep);
	int nr=sstr->read(winstep,&b[winkeep]);
	// If it read nothing, this is not a new window, so return nothing.
	if(!nr){
		return 0;
	}
	// If it read less than a full step, this is the last window.
	if(nr<winstep){
		end=true;
	}
	return winkeep+nr;
}

/*
createBufferWindow
	Creates a sequence stream and window for a buffer.
	Returns false on failure.
*/
bool createBufferWindow(char*buf,int bufs,seqStreamBuffer*&ssb,seqStreamWindow*&ssw,int wsize,int wstep){
	ssb=seqStreamBuffer::create();
	if(!ssb){
		return false;
	}
	if(!ssb->arm(buf,bufs)){
		delete ssb;
		return false;
	}
	ssw=seqStreamWindow::create(ssb,wsize,wstep);
	if(!ssw){
		delete ssb;
		return false;
	}
	return true;
}

