////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "aux.hpp"
#include "sequences.hpp"
#include "sequencelist.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Sequence classes

std::vector<seqClass> seqClasses;

size_t nSequenceClasses(){
	return seqClasses.size();
}

seqClass*getSeqClassByName(char*name){
	double cls;
	if(name[0]=='N'&&!name[1])cls=0.0;
	else if(name[0]=='+'&&!name[1])cls=1.0;
	else if(name[0]=='-'&&!name[1])cls=-1.0;
	else if(name[0]=='+'&&name[1]=='+'&&!name[2])cls=2.0;
	else if(name[0]=='-'&&name[1]=='-'&&!name[2])cls=-2.0;
	else cls=strtod(name,0);
	for(auto& c: seqClasses)
		if(c.cls == cls)
			return &c;
	cmdError("Class not found.");
	cout << t_indent << "Requested class: " << name << "\n";
	return 0;
}

void printSeqClasses(){
	cmdSection("Classes");
	for(auto& c: seqClasses)
		cout << t_indent << c.name << t_indent << "Value: " << c.cls << t_indent << "Binary flag: " << (c.flag?"+":"-") << "\n";
}

seqClass*registerSeqClass(double cls,char*name,bool flag){
	for(auto& c: seqClasses){
		if(c.cls==cls){
			cmdError("Class with the same identifier specified twice.");
			return 0;
		}else if(!strcmp(c.name,name)){
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

void seqList::printInfo(char*title){
	cmdSection(title);
	seqListSeq*s=seq;
	for(int l=0;l<nseq;l++,s++){
		cout << t_indent << "\"" << s->name << "\"" << t_indent << "Class: " << s->cls->name << t_indent << "Length: " << s->bufs << " bp" << t_indent << "Training mode: " << (s->trainMode==train_Windows?"Windows":"Full, normalized") << "\n";
	}
}

