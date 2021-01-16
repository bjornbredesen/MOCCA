////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "vaux.hpp"
#include "sequences.hpp"
#include "sequencelist.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Sequence classes

std::vector<seqClass> seqClasses;

size_t nSequenceClasses(){
	return seqClasses.size();
}

seqClass*getSeqClassByName(std::string name){
	double cls;
	if(name[0]=='N'&&!name[1])cls=0.0;
	else if(name[0]=='+'&&!name[1])cls=1.0;
	else if(name[0]=='-'&&!name[1])cls=-1.0;
	else if(name[0]=='+'&&name[1]=='+'&&!name[2])cls=2.0;
	else if(name[0]=='-'&&name[1]=='-'&&!name[2])cls=-2.0;
	else cls=strtod(name.c_str(), 0);
	for(auto& c: seqClasses)
		if(c.cls == cls)
			return &c;
	cmdError("Class not found.");
	cout << t_indent << "Requested class: " << name << "\n";
	return 0;
}

seqClass*getSeqClassByValue(double cls){
	for(auto& c: seqClasses)
		if(c.cls == cls)
			return &c;
	cmdError("Class not found.");
	cout << t_indent << "Requested class: " << cls << "\n";
	return 0;
}

void printSeqClasses(){
	cmdSection("Classes");
	for(auto& c: seqClasses)
		cout << t_indent << c.name << t_indent << "Value: " << c.cls << t_indent << "Binary flag: " << (c.flag?"+":"-") << "\n";
}

seqClass*registerSeqClass(double cls,std::string name,bool flag){
	for(auto& c: seqClasses){
		if(c.cls==cls){
			cmdError("Class with the same identifier specified twice.");
			return 0;
		}else if(!c.name.compare(name)){
			cmdError("Class with the same name specified twice.");
			return 0;
		}
	}
	seqClass c = seqClass {
		cls,
		name,
		flag
	};
	seqClasses.push_back(c);
	return &seqClasses.back();
}

////////////////////////////////////////////////////////////////////////////////////
// Sequence list

e_trainMode getTrainModeByName(char*name){
	if(!strcmp(name,"win")){
		return train_Windows;
	}else if(!strcmp(name,"full")){
		return train_Full;
	}else{
		cmdError("Invalid training mode.");
		return train_Invalid;
	}
}

seqList::seqList(){
	nseq=0;
	seq=0;
	nPos=nNeg=0;
	ownSeq=true;
}

seqListSeq*seqList::addSeq(char*name,char*buf,int bufs,seqClass*cls,e_trainMode tm){
	if(!buf||!bufs)return 0;
	nseq++;
	seq=(seqListSeq*)realloc(seq,sizeof(seqListSeq)*nseq);
	if(!seq){
		outOfMemory();
		return 0;
	}
	seqListSeq*r=&seq[nseq-1];
	r->name=name;
	r->buf=buf;
	r->bufs=bufs;
	r->cls=cls;
	r->trainMode=tm;
	if(cls->flag){
		nPos++;
	}else{
		nNeg++;
	}
	return r;
}

seqList::~seqList(){
	if(seq){
		if(ownSeq){
			seqListSeq*s=seq;
			for(int l=0;l<nseq;l++,s++){
				if(s->name)free(s->name);
				if(s->buf)free(s->buf);
			}
		}
		free(seq);
	}
}

seqList*seqList::create(){
	seqList*r=new seqList();
	if(!r){
		outOfMemory();
		return 0;
	}
	return r;
}

seqList*seqList::create(int n){
	seqList*r=new seqList();
	if(!r){
		outOfMemory();
		return 0;
	}
	r->seq=(seqListSeq*)realloc(r->seq,sizeof(seqListSeq)*n);
	if(!r->seq){
		outOfMemory();
		delete r;
		return 0;
	}
	memset(r->seq,0,sizeof(seqListSeq)*n);
	return r;
}

bool seqList::loadFastaBatch(char*path,seqClass*cls,e_trainMode tm){
	if(!cls||tm==train_Invalid)return false;
	seqStreamFastaBatch*ssfb=seqStreamFastaBatch::load(path);
	if(!ssfb){
		return false;
	}
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb->getBlock());){
		char*sname=ssfbblk->getName();
		if(!sname)sname=(char*)"Unnamed";
		char*buf=0;
		int bufs=(int)ssfbblk->buffer(buf);
		if(!addSeq(cloneString(sname),buf,bufs,cls,tm)){
			free(buf);
			delete ssfb;
			return false;
		}
	}
	delete ssfb;
	return true;
}

bool seqList::addRandomIid(char*tpath,int nadd,int len,seqClass*cls,e_trainMode tm){
	if(!cls||tm==train_Invalid)return false;
	if(len<=0){
		cmdError("Invalid random sequence length requested.");
		return false;
	}else if(nadd<=0){
		cmdError("Invalid number of random sequences requested.");
		return false;
	}
	// Train generative model
	autodelete<seqStreamRandomIid> rss(new seqStreamRandomIid());
	if(!rss.ptr){
		outOfMemory();
		return false;
	}
	autodelete<seqStreamFastaBatch> ssfb(seqStreamFastaBatch::load(tpath));
	if(!ssfb.ptr){
		return false;
	}
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb.ptr->getBlock());){
		rss.ptr->train(ssfbblk);
	}
	// Generate and add random sequences
	for(int l=0;l<nadd;l++){
		autofree<char> buf((char*)malloc(len));
		if(!buf.ptr){
			outOfMemory();
			return false;
		}
		rss.ptr->read(len,buf.ptr);
		ostringstream oss;
		int seed=rand();
		srand(seed);
		oss << "Random (i.i.d.) " << nseq << "(seed=" << seed << ")";
		if(!addSeq(cloneString((char*)oss.str().c_str()),buf.disown(),len,cls,tm)){
			return false;
		}
	}
	return true;
}

bool seqList::addRandomMC(char*tpath,int nadd,int len,seqClass*cls,e_trainMode tm,int order){
	if(!cls||tm==train_Invalid)return false;
	if(len<=0){
		cmdError("Invalid random sequence length requested.");
		return false;
	}else if(nadd<=0){
		cmdError("Invalid number of random sequences requested.");
		return false;
	}
	// Train generative model
	autodelete<seqStreamRandomMC> rss(new seqStreamRandomMC(order));
	if(!rss.ptr){
		outOfMemory();
		return false;
	}
	autodelete<seqStreamFastaBatch> ssfb(seqStreamFastaBatch::load(tpath));
	if(!ssfb.ptr){
		return false;
	}
	for(seqStreamFastaBatchBlock*ssfbblk;(ssfbblk=ssfb.ptr->getBlock());){
		rss.ptr->train(ssfbblk);
	}
	// Generate and add random sequences
	for(int l=0;l<nadd;l++){
		autofree<char> buf((char*)malloc(len));
		if(!buf.ptr){
			outOfMemory();
			return false;
		}
		rss.ptr->read(len,buf.ptr);
		ostringstream oss;
		int seed=rand();
		srand(seed);
		oss << "Random (MC order " << order << ") " << nseq << "(seed=" << seed << ")";
		if(!addSeq(cloneString((char*)oss.str().c_str()),buf.disown(),len,cls,tm)){
			return false;
		}
	}
	return true;
}

bool seqList::addClone(seqListSeq*sls){
	autofree<char> buf((char*)malloc(sls->bufs));
	if(!buf.ptr){
		outOfMemory();
		return false;
	}
	memcpy(buf.ptr, sls->buf, sls->bufs);
	if(!addSeq(cloneString((char*)sls->name),buf.disown(),sls->bufs,sls->cls,sls->trainMode)){
		return false;
	}
	return true;
}

void seqList::printInfo(char*title){
	cmdSection(title);
	seqListSeq*s=seq;
	for(int l=0;l<nseq;l++,s++){
		cout << t_indent << "\"" << s->name << "\"" << t_indent << "Class: " << s->cls->name << t_indent << "Length: " << s->bufs << " bp" << t_indent << "Training mode: " << (s->trainMode==train_Windows?"Windows":"Full, normalized") << "\n";
	}
}

