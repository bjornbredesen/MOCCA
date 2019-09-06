////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "./lib/rapidxml-1.13/rapidxml.hpp"
using namespace rapidxml;
#include "config.hpp"
#include "aux.hpp"
#include "motifs.hpp"

////////////////////////////////////////////////////////////////////////////////////
// Motif list

motifList::motifList(){
	motifs=0;
	nmotifs=0;
	maxLen=0;
}

motifList::~motifList(){
	if(motifs){
		motifListMotif*m=motifs;
		for(int l=0;l<nmotifs;l++,m++){
			if(m->name)free(m->name);
			if(m->data){
				switch(m->type){
					case motifType_IUPAC:{
						IUPACMotif*d=(IUPACMotif*)m->data;
						if(d->seq)free(d->seq);
						break;}
					default:
						cmdWarning("Invalid motif.");
				}
				free(m->data);
			}
		}
		free(motifs);
	}
}

motifList*motifList::create(){
	motifList*r=new motifList();
	if(!r){
		outOfMemory();
		return 0;
	}
	return r;
}

motifListMotif*motifList::addMotif(char*name,motifType type,void*data,int len){
	char*cname=cloneString(name);
	if(!cname)return 0;
	nmotifs++;
	motifs=(motifListMotif*)realloc(motifs,sizeof(motifListMotif)*nmotifs);
	if(!motifs){
		outOfMemory();
		free(cname);
		return 0;
	}
	motifListMotif*m=&motifs[nmotifs-1];
	memset(m,0,sizeof(motifListMotif));
	m->name=cname;
	m->index=nmotifs-1;
	m->len=len;
	m->type=type;
	m->data=data;
	maxLen=max(len,maxLen);
	return m;
}

bool motifList::addMotifsFromXML(char*path){
	FILE*f=fopen(path,"rb");
	if(!f){
		cout << m_error << "Could not open file \"" << path << "\".\n";
		return false;
	}
	fseek(f,0,SEEK_END);
	long fsize=ftell(f);
	fseek(f,0,SEEK_SET);
	autofree<char>_xmltxt;
	if(!_xmltxt.resize(fsize+1)){
		fclose(f);
		return false;
	}
	char*xmltxt=_xmltxt.ptr;
	if(!fread(xmltxt,1,fsize,f)){
		fclose(f);
		cmdError("File reading failed.");
		return false;
	}
	fclose(f);
	xmltxt[fsize]=0;
	xml_document<>*doc=new xml_document<>();
	doc->parse<0>(xmltxt);
	// Read motif definitions
	xml_node<>*node=doc->first_node("motifs");
	if(node){
		xml_node<>*m=node->first_node();
		while(m){
			xml_attribute<>*name=m->first_attribute("name");
			xml_attribute<>*seq=m->first_attribute("seq");
			xml_attribute<>*nmis=m->first_attribute("misses");
			if(name&&seq&&nmis){
				if(!addIUPACMotif(name->value(),seq->value(),(int)strtol(nmis->value(),0,10))){
					delete doc;
					return false;
				}
			}
			m=m->next_sibling();
		}
	}
	delete doc;
	return true;
}

motifListMotif*motifList::addIUPACMotif(char*name,char*seq,int nmis,bool allowDuplicates){
	if(!name||!seq){
		cmdError("addIUPACMotif: Null-pointer arguments.");
		return 0;
	}
	if(nmis<0){
		cmdError("addIUPACMotif: Negative number of allowed mismatches.");
		return 0;
	}
	// Skip motifs already in list (based on sequence)
	if(!allowDuplicates)for(int l=0;l<nmotifs;l++){
		if(motifs[l].type!=motifType_IUPAC)continue;
		IUPACMotif*mot=(IUPACMotif*)motifs[l].data;
		int sl=(int)strlen(seq);
		if(sl!=(int)strlen(mot->seq))continue;
		if(!strcmp(mot->seq,seq))return &motifs[l];
		bool eql=true;
		for(int x=0,y=sl-1;x<sl;x++,y--){
			char c=mot->seq[y];
			if(c=='T')c='A';
			else if(c=='A')c='T';
			else if(c=='G')c='C';
			else if(c=='C')c='G';
			// IUPAC extension
			else if(c=='R')c='Y';
			else if(c=='Y')c='R';
			else if(c=='K')c='M';
			else if(c=='M')c='K';
			else if(c=='B')c='V';
			else if(c=='V')c='B';
			else if(c=='D')c='H';
			else if(c=='H')c='D';
			//else if(c=='S')c='S';
			//else if(c=='W')c='W';
			//else if(c=='N')c='N';
			if(seq[x]!=c){eql=false;break;}
		}
		if(eql)return &motifs[l];
	}
	IUPACMotif*d=(IUPACMotif*)malloc(sizeof(IUPACMotif));
	if(!d){
		outOfMemory();
		return 0;
	}
	memset(d,0,sizeof(IUPACMotif));
	d->seq=cloneString(seq);
	if(!d->seq){
		outOfMemory();
		free(d);
		return 0;
	}
	d->nmis=nmis;
	int mlen=(int)strlen(seq);
	motifListMotif*r=addMotif(name,motifType_IUPAC,d,mlen);
	return r;
}

bool motifList::addKMers(int k,bool allowDuplicates){
	if(k<0){
		cmdError("Invalid k-mer arguments.");
		return false;
	}
	char nt[4]={'A','C','G','T'};
	int nkm=(int)pow(4,(double)k);
	autofree<char> kmer;
	kmer.resize(k+1);
	if(!kmer.ptr){
		outOfMemory();
		return false;
	}
	kmer[k]=0;
	for(int i=0;i<nkm;i++){
		for(int l=0;l<k;l++){
			kmer[l]=nt[(i>>(l*2))&3];
		}
		if(!addIUPACMotif(kmer.ptr,kmer.ptr,0,allowDuplicates)){
			return false;
		}
	}
	return true;
}

bool motifList::addRandom(int n,int len){
	#define nRNDMotifNT (15+12)
	if(n<=0||len<=0){
		cmdError("Invalid random motif arguments.");
		return false;
	}
	char nt[nRNDMotifNT]={'A','C','G','T',  'A','C','G','T',  'A','C','G','T',  'A','C','G','T',
		'R','Y','S','W','K','M','B','D','H','V','N'};
	autofree<char> kmer;
	kmer.resize(len+1);
	if(!kmer.ptr){
		outOfMemory();
		return false;
	}
	kmer[len]=0;
	int onmotifs = nmotifs;
	for(;;){
		for(int l=0;l<len;l++)
			kmer[l]=nt[rand()%nRNDMotifNT];
		if(!addIUPACMotif(kmer.ptr,kmer.ptr,0,false)){
			return false;
		}
		if(nmotifs>=onmotifs+n)
			break;
	}
	return true;
}

void motifList::printInfo(){
	motifListMotif*m=motifs;
	cmdSection("Motifs");
	for(int l=0;l<nmotifs;l++,m++){
		cout << t_indent << m->name;
		switch(m->type){
			case motifType_IUPAC:
				cout << " (IUPAC)\n";
				if(!m->data){
					cout << t_indent << t_indent << "Corrupted\n";
				}else{
					IUPACMotif*d=(IUPACMotif*)m->data;
					cout << t_indent << t_indent << "Sequence: \"" << d->seq << "\"\n";
					cout << t_indent << t_indent << "# mismatches: " << d->nmis << "\n";
				}
				break;
			default:
				cout << " (Invalid)\n";
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif occurrence container

motifOccContainer::motifOccContainer(motifList*ml){
	iAFirst=iALast=-1;
	tblSize=0;
	occTbl=0;
	iAFirstT=0;
	iALastT=0;
	
	motifs=ml;
	nOcc=0;
	nOccT=0;
}

bool motifOccContainer::growTable(int ngrow){
	if(!ngrow)return true;
	int ptblSize=tblSize;
	tblSize+=ngrow;
	occTbl=(motifOcc*)realloc(occTbl,sizeof(motifOcc)*tblSize);
	if(!occTbl){
		outOfMemory();
		return false;
	}
	motifOcc*o=&occTbl[ptblSize];
	iFree=ptblSize;
	for(int l=ptblSize;l<tblSize;l++,o++){
		o->iFNext=(l==tblSize-1)?-1:l+1;
		o->active=false;
	}
	return true;
}

motifList*motifOccContainer::getMotifList(){
	return motifs;
}

motifOccContainer::~motifOccContainer(){
	motifOcc*o=getFirst();
	while(o){
		if(o->extra_buffer){
			free(o->extra_buffer);
			o->extra_buffer=0;
		}
		o=getNext(o);
	}
	if(nOccT)free(nOccT);
	if(iAFirstT)free(iAFirstT);
	if(iALastT)free(iALastT);
	if(occTbl)free(occTbl);
}

motifOccContainer*motifOccContainer::create(int basesize,motifList*ml){
	if(!ml||!ml->nmotifs)return 0;
	motifOccContainer*r=new motifOccContainer(ml);
	if(!r){
		outOfMemory();
		return 0;
	}
	if(!r->growTable(basesize)){
		delete r;
		return 0;
	}
	// Allocate type-specific occurrence count table,
	r->nOccT=(int*)malloc(sizeof(int)*ml->nmotifs);
	if(!r->nOccT){
		outOfMemory();
		delete r;
		return 0;
	}
	memset(r->nOccT,0,sizeof(int)*ml->nmotifs);
	// Allocate type-specific occurrence index tables.
	r->iAFirstT=(int*)malloc(sizeof(int)*ml->nmotifs);
	if(!r->iAFirstT){
		outOfMemory();
		delete r;
		return 0;
	}
	r->iALastT=(int*)malloc(sizeof(int)*ml->nmotifs);
	if(!r->iALastT){
		outOfMemory();
		delete r;
		return 0;
	}
	// Initialize types as empty.
	for(int l=0;l<ml->nmotifs;l++){
		r->iAFirstT[l]=-1;
		r->iALastT[l]=-1;
	}
	return r;
}

motifOcc*motifOccContainer::createMotifOcc(long long start,motifListMotif*m,bool strand){
	if(!m)return 0;
	if(iFree==-1){
		if(!growTable(tblSize)){
			return 0;
		}
	}
	motifOcc*ret=&occTbl[iFree];
	if(ret->active){
		cmdError("Tried to create occurrence using active occurrence.");
		return 0;
	}
	ret->iANext=iAFirst;
	if(iAFirst!=-1)occTbl[iAFirst].iAPrev=iFree;
	iAFirst=iFree;
	if(iALast==-1)iALast=iFree;
	//
	ret->iANextT=-1;
	ret->iAPrevT=-1;
	int tind=m->index;
	if(iAFirstT){
		ret->iANextT=iAFirstT[tind];
		if(iAFirstT[tind]!=-1)occTbl[iAFirstT[tind]].iAPrevT=iFree;
		iAFirstT[tind]=iFree;
	}
	if(iALastT&&iALastT[tind]==-1)iALastT[tind]=iFree;
	nOccT[tind]++;
	nOcc++;
	//
	iFree=ret->iFNext;
	ret->active=true;
	ret->iAPrev=-1;
	ret->mot=m;
	ret->start=start;
	ret->strand=strand;
	ret->skip=false;
	ret->extra_buffer=0;
	return ret;
}

bool motifOccContainer::freeMotifOcc(motifOcc*o){
	if(!o->active){
		cmdError("Tried to free motif occurence twice.");
		return false;
	}
	if(o->iANext!=-1)occTbl[o->iANext].iAPrev=o->iAPrev;
	if(o->iAPrev!=-1)occTbl[o->iAPrev].iANext=o->iANext;
	int ind=int(o-occTbl);
	if(ind==iAFirst)iAFirst=o->iANext;
	if(ind==iALast)iALast=o->iAPrev;
	//
	int tind=o->mot->index;
	if(o->iANextT!=-1)occTbl[o->iANextT].iAPrevT=o->iAPrevT;
	if(o->iAPrevT!=-1)occTbl[o->iAPrevT].iANextT=o->iANextT;
	if(iAFirstT&&ind==iAFirstT[tind])iAFirstT[tind]=o->iANextT;
	if(iALastT&&ind==iALastT[tind])iALastT[tind]=o->iAPrevT;
	//
	o->iFNext=iFree;
	iFree=ind;
	o->active=false;
	nOcc--;
	nOccT[o->mot->index]--;
	if(o->iFNext!=-1&&occTbl[o->iFNext].active)cmdWarning("Next free motif occurrence pointer to an active motif occurrence.");
	if(o->extra_buffer){
		free(o->extra_buffer);
		o->extra_buffer=0;
	}
	return true;
}

void motifOccContainer::flush(){
	iAFirst=iALast=-1;
	for(int l=0;l<motifs->nmotifs;l++){
		iAFirstT[l]=-1;
		iALastT[l]=-1;
		nOccT[l]=0;
	}
	nOcc=0;
	motifOcc*o=occTbl;
	for(int l=0;l<tblSize-1;l++,o++){
		o->iFNext=l+1;
		o->active=false;
	}
	occTbl[tblSize-1].iFNext=-1;
	occTbl[tblSize-1].active=false;
	iFree=0;
}

motifOcc*motifOccContainer::getFirst(){
	return iAFirst==-1?0:&occTbl[iAFirst];
}

motifOcc*motifOccContainer::getLast(){
	return iALast==-1?0:&occTbl[iALast];
}

motifOcc*motifOccContainer::getNext(motifOcc*o){
	return o->iANext==-1?0:&occTbl[o->iANext];
}

motifOcc*motifOccContainer::getPrev(motifOcc*o){
	return o->iAPrev==-1?0:&occTbl[o->iAPrev];
}

motifOcc*motifOccContainer::getFirst(int t){
	return (t>=0&&t<motifs->nmotifs&&iAFirstT&&iAFirstT[t]!=-1)?&occTbl[iAFirstT[t]]:0;
}

motifOcc*motifOccContainer::getLast(int t){
	return (t>=0&&t<motifs->nmotifs&&iALastT&&iALastT[t]!=-1)?&occTbl[iALastT[t]]:0;
}

motifOcc*motifOccContainer::getNextSame(motifOcc*o){
	return o->iANextT==-1?0:&occTbl[o->iANextT];
}

motifOcc*motifOccContainer::getPrevSame(motifOcc*o){
	return o->iAPrevT==-1?0:&occTbl[o->iAPrevT];
}

////////////////////////////////////////////////////////////////////////////////////
// IUPAC table

#define iupacTblW ('Z'-'A'+1)
char iupacTbl[iupacTblW][iupacTblW];
#define iupac(X,Y) iupacTbl[X-'A'][Y-'A']
void initIUPACTbl(){
	memset(iupacTbl,0,sizeof(char)*iupacTblW*iupacTblW);
	iupac('A','A')=iupac('C','C')=iupac('G','G')=iupac('T','T')=1;
	iupac('N','A')=iupac('N','C')=iupac('N','G')=iupac('N','T')=1;
	iupac('R','A')=iupac('R','G')=1;
	iupac('Y','C')=iupac('Y','T')=1;
	iupac('S','G')=iupac('S','C')=1;
	iupac('W','A')=iupac('W','T')=1;
	iupac('K','G')=iupac('K','T')=1;
	iupac('M','A')=iupac('M','C')=1;
	iupac('B','C')=iupac('B','G')=iupac('B','T')=1;
	iupac('D','A')=iupac('D','G')=iupac('D','T')=1;
	iupac('H','A')=iupac('H','C')=iupac('H','T')=1;
	iupac('V','A')=iupac('V','C')=iupac('V','G')=1;
}

////////////////////////////////////////////////////////////////////////////////////
// Finite state machine

char FSMindNT[4]={'A','C','G','T'};
char FSMindNTC[4]={'T','G','C','A'};

/*
motifFSMNodeMotif
	Motif node motif.
	Used during construction to keep track of parsing states.
	Used after construction for registering motif occurrences.
*/
class motifFSMNodeMotif{
public:
	motifListMotif*mot;		// Motif.
	IUPACMotif*iumot;
	int nmis;				// Mismatches so far.
	int readpos;				// Reading position.
	bool com;				// Strand.
	motifFSMNodeMotif*next;	// Next node motif in the list.
	motifFSMNodeMotif(motifListMotif*m,int nm,int rp,bool c,motifFSMNodeMotif*_next){
		mot=m;
		iumot=(IUPACMotif*)mot->data;
		nmis=nm;
		readpos=rp;
		com=c;
		next=_next;
	}
	/*
	match
		Gives match state when appending a nucleotide.
		True if still a match, or otherwise false.
	*/
	bool match(int NTind,int&nnm,int&nrp){
		nrp=readpos+1;
		if(nrp>=mot->len)return false; // End of motif, so no longer a match.
		nnm=nmis;
		char m=iumot->seq[com?(mot->len-1-nrp):nrp];
		char s=com?FSMindNTC[NTind]:FSMindNT[NTind];
		if(!iupac(m,s)){
			nnm++;
		}
		return (nnm<=iumot->nmis); // Match depends on whether it has exceeded the allowed number of mismatches.
	}
};

/*
motifFSMNode
	Finite state machine node.
*/
class motifFSMNode{
public:
	motifFSMNodeMotif*mot;
	motifFSMNode*next[4];	// Next nodes, indexed for nucleotides.
	// Base constructor
	motifFSMNode(){
		memset(next,0,sizeof(motifFSMNode*)*4);
		mot=0;
	}
	// Constructor for node extension
	motifFSMNode(motifFSMNode*parent,motifList*motifs,int NTind){
		memset(next,0,sizeof(motifFSMNode*)*4);
		mot=0;
		// Add extended motifs where possible
		motifFSMNodeMotif*pmot=parent->mot;
		int nrp,nnm;
		while(pmot){
			if(pmot->match(NTind,nnm,nrp)){
				mot=new motifFSMNodeMotif(pmot->mot,nnm,nrp,pmot->com,mot);
			}
			pmot=pmot->next;
		}
		// Add any new motifs
		bool match;
		motifListMotif*bmot=motifs->motifs;
		for(int l=0;l<motifs->nmotifs;l++,bmot++){
			if(bmot->skip||bmot->type!=motifType_IUPAC||!bmot->data)continue;
			IUPACMotif*iumot=(IUPACMotif*)bmot->data;
			match=iupac(iumot->seq[0],FSMindNT[NTind]);
			if(iumot->nmis||match){
				mot=new motifFSMNodeMotif(bmot,match?0:1,0,false,mot);
			}
			match=iupac(iumot->seq[bmot->len-1],FSMindNTC[NTind]);
			if(iumot->nmis||match){
				mot=new motifFSMNodeMotif(bmot,match?0:1,0,true,mot);
			}
		}
	}
	~motifFSMNode(){
		motifFSMNodeMotif*m=mot,*pm;
		while(m){
			pm=m;
			m=m->next;
			delete pm;
		}
	}
	/*
	eql
		Returns true if this node is equivalent to a given node.
	*/
	bool eql(motifFSMNode*n){
		motifFSMNodeMotif*m1=mot,*m2;
		bool found;
		int nm1=0,nm2=0;
		// For every motif
		while(m1){
			found=false;
			m2=n->mot;
			// Search for equivalent
			while(m2){
				// I.e., same motif, strand, reading position and number of mismatches so far.
				if(m2->mot==m1->mot&&m2->com==m1->com&&m2->readpos==m1->readpos&&m2->nmis==m1->nmis){
					found=true;
					break;
				}
				m2=m2->next;
			}
			// If not found, then it is a mismatch.
			if(!found)return false;
			m1=m1->next;
			nm1++;
		}
		// Also make sure there are the same number of motifs (otherwise the other one might have some unmatched motifs although all in this one are just as long)
		m2=n->mot;
		while(m2){
			m2=m2->next;
			nm2++;
		}
		if(nm1!=nm2)return false;
		return true;
	}
	/*
	finalize
		Finishes construction by flushing away ineffective motifs.
		The final list of node motifs will correspond to registering occurrences.
	*/
	void finalize(){
		motifFSMNodeMotif*m=mot,*pm;
		mot=0;
		while(m){
			if(m->readpos==m->mot->len-1){
				mot=new motifFSMNodeMotif(m->mot,m->nmis,m->readpos,m->com,mot);
			}
			pm=m;
			m=m->next;
			delete pm;
		}
	}
	/*
	addMotifs
		Registers any motif occurrences for this node.
	*/
	inline void addMotifs(motifOccContainer*oc,long long pos){
		motifFSMNodeMotif*m=mot;
		while(m){
			if(m->readpos==m->mot->len-1){
				oc->createMotifOcc(pos-m->mot->len+1,m->mot,m->com);
			}
			m=m->next;
		}
	}
};

/*
motifFSMQueueEdge
	Edge for finite state graph extension queue.
*/
class motifFSMQueueEdge{
public:
	motifFSMQueueEdge*next;		// Pointer to next edge in queue.
	motifFSMNode*base;			// Base node.
	int nextind;					// Nucleotide index.
	motifFSMQueueEdge(motifFSMNode*b,int ni,motifFSMQueueEdge*n){
		base=b;
		nextind=ni;
		next=n;
	}
};

/*
motifFSM
	Class for motif finite state machine.
*/
class motifFSM{
private:
	motifFSMNode**nodes;			// Nodes.
	motifFSMNode*rootNode;			// Root node.
	motifFSMNode*state;				// Current state node.
	int nnodes;						// Number of nodes.
	motifFSMQueueEdge*edgeQueue;	// Edge queue.
	// Private constructor
	motifFSM(){
		state=0;
		nnodes=0;
		nodes=0;
		edgeQueue=0;
	}
	/*
	insertNode
		Tries to insert a node into the graph.
		The node given should always be a new one.
		If a corresponding node already is in the graph, it deletes the
		new one and returns the old one.
		Thus, the return is the new node or a corresponding one,
		ensuring uniqueness of nodes.
	*/
	motifFSMNode*insertNode(motifFSMNode*n){
		if(!n){
			outOfMemory();
			// A new node is always sent, so if it is a null-pointer, it is out of memory.
			return 0;
		}
		// If a corresponding node already exists, delete the new one and return the old one.
		for(int l=0;l<nnodes;l++){
			if(nodes[l]->eql(n)){
				delete n;
				return nodes[l];
			}
		}
		// Otherwise, insert into graph.
		nnodes++;
		nodes=(motifFSMNode**)realloc(nodes,sizeof(motifFSMNode*)*nnodes);
		nodes[nnodes-1]=n;
		for(int l=0;l<4;l++){
			edgeQueue=new motifFSMQueueEdge(n,l,edgeQueue);
		}
		return n;
	}
public:
	/*
	construct
		Call to construct.
	*/
	static motifFSM*construct(motifList*motifs){
		motifFSM*r=new motifFSM();
		if(!r){
			outOfMemory();
			return 0;
		}
		cmdTask task((char*)"Constructing motif Finite-State Machine");
		cmdTask::refresh();
		r->rootNode=r->insertNode(new motifFSMNode());
		if(!r->rootNode){
			delete r;
			return 0;
		}
		// Process edges in queue
		motifFSMQueueEdge*e;
		while(r->edgeQueue){
			e=r->edgeQueue;
			r->edgeQueue=r->edgeQueue->next;
			// Try to insert extended node (if not original it will be substituted by an old one via return).
			motifFSMNode*n=r->insertNode(new motifFSMNode(e->base,motifs,e->nextind));
			if(!n){
				delete r;
				return 0;
			}
			// Set base node's corresponding next-pointer
			e->base->next[e->nextind]=n;
			delete e;
		}
		for(int l=0;l<r->nnodes;l++){
			if(r->nodes[l])
				r->nodes[l]->finalize();
		}
		r->flush();
		cmdTask::wipe();
		cmdTaskComplete("Constructing motif Finite-State Machine");
		cout << t_indent << t_indent << "Nodes: " << r->nnodes << "\n";
		int nmotifsused=0;
		for(int l=0;l<motifs->nmotifs;l++)
			if(!motifs->motifs[l].skip&&motifs->motifs[l].type==motifType_IUPAC)nmotifsused++;
		cout << t_indent << t_indent << "Motifs: " << nmotifsused << "\n";
		return r;
	}
	~motifFSM(){
		if(nodes){
			for(int l=0;l<nnodes;l++){
				if(nodes[l])
					delete nodes[l];
			}
			free(nodes);
		}
		motifFSMQueueEdge*e=edgeQueue,*pe;
		while(e){
			pe=e;
			e=e->next;
			delete pe;
		}
	}
	/*
	flush
		Resets parsing to base state.
	*/
	inline void flush(){
		state=rootNode;
	}
	/*
	feed
		Feeds a nucleotide to the finite state machine.
	*/
	void feed(char nt,motifOccContainer*oc,long long pos){
		if(!state){
			cout << "Error: Illegal state.\n";
			return;
		}
		int nti;
		switch(nt){
			case 'A':nti=0;break;
			case 'C':nti=1;break;
			case 'G':nti=2;break;
			case 'T':nti=3;break;
			case 'N':flush();return;
			default:cout << "Warning: Unrecognized character, '" << nt << "', fed to finite state machine.\n";flush();return;
		}
		state=state->next[nti];
		if(!state){
			cout << "Error: State transition resulted in dead end. Resetting.\n";flush();return;
		}
		state->addMotifs(oc,pos);
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif window

motifWindow::motifWindow(motifList*_motifs){
	occContainer=0;
	mFSM=0;
	motifs=_motifs;
}

bool motifWindow::initialize(){
	occContainer=motifOccContainer::create(128,motifs);
	if(!occContainer){
		return false;
	}
	int nIUPAC=0;
	motifListMotif*m=motifs->motifs;
	for(int l=0;l<motifs->nmotifs;l++,m++){
		switch(m->type){
			case motifType_IUPAC:
				nIUPAC++;
				break;
			default:
				cmdError("Invalid motif");
				return false;
		}
	}
	if(nIUPAC){
		if(getConfiguration()->useFSM){
			mFSM=motifFSM::construct(motifs);
			if(!mFSM){
				return 0;
			}
		}
	}
	return true;
}

bool motifWindow::motifMatchIUPAC(char*seq,int seqlen,motifListMotif*mot,bool com){
	IUPACMotif*iumot=(IUPACMotif*)mot->data;
	char*mots=iumot->seq;
	char s,m;
	if(com){
		mots+=mot->len-1;
	}
	int nmatch=0;
	int nmis=iumot->nmis;
	for(int l=0;l<min(mot->len,seqlen);l++,seq++){
		s=*seq;
		m=*mots;
		if(com)switch(s){
			case 'A':s='T';break;
			case 'T':s='A';break;
			case 'C':s='G';break;
			case 'G':s='C';break;
		}
		if(iupac(m,s)){
			nmatch++;
		}else{
			nmis--;
			if(nmis<0)return false;
		}
		if(com){
			mots--;
		}else{
			mots++;
		}
	}
	return (nmatch+iumot->nmis>=mot->len);
}

motifWindow::~motifWindow(){
	if(occContainer)delete occContainer;
}

motifWindow*motifWindow::create(motifList*_motifs){
	if(!_motifs){
		cmdError("Motif list null-pointer.");
	}
	if(!_motifs->nmotifs){
		cmdError("Empty motif list.");
		return 0;
	}
	motifWindow*r=new motifWindow(_motifs);
	if(!r){
		outOfMemory();
		return 0;
	}
	if(!r->initialize()){
		delete r;
		return 0;
	}
	return r;
}

bool motifWindow::flush(){
	if(!occContainer){
		cmdError("Occurrence container not created.");
		return false;
	}
	occContainer->flush();
	if(mFSM)mFSM->flush();
	wPos=0;
	wLen=0;
	return true;
}

bool motifWindow::readWindow(char*wseq,long long wpos,int wlen){
	if(!occContainer){
		cmdError("Occurrence container not created.");
		return false;
	}
	bool wskip=false;
	if(wLen){
		if(wpos<wPos||wpos>=wPos+wLen){
			flush();
		}else{
			wskip=true;
			motifOcc*o=occContainer->getFirst(),*co;
			while(o){
				co=o;
				o=occContainer->getNext(o);
				if(co->start<wpos){
					occContainer->freeMotifOcc(co);
				}
			}
		}
	}
	int wstartbase=int(wPos+(long long)wLen-wpos); // Start at the end of the previous window, localized to the new one.
	wPos=wpos;
	wLen=wlen;
	if(mFSM){
		// Parse IUPAC with FSM
		int wstart=0;
		if(wskip){
			wstart=wstartbase;
			if(wstart<0)wstart=0;
		}
		// Run through the sequence.
		char*b=wseq+wstart;
		for(unsigned int rcoord=wstart;int(rcoord)<int(wlen);b++,rcoord++){
			mFSM->feed(*b,occContainer,wpos+(long long)rcoord);
		}
	}else{
		// Parse IUPAC with naive parsing
		// Run through the motifs.
		motifListMotif*m=motifs->motifs;
		for(int mtype=0;mtype<motifs->nmotifs;mtype++,m++){
			if(m->skip||m->type!=motifType_IUPAC||!m->data)continue;
			// Skip any sequence portion which should already have been processed.
			int wstart=0;
			if(wskip){
				wstart=wstartbase-m->len+1;
				if(wstart<0)wstart=0;
			}
			// Run through the sequence.
			char*b=wseq+wstart;
			// The sequence should be scanned until (and including): rcoord=wlen-m->slen,
			// since at that point it will read the final, full motif occurrence for this window.
			// I.e., it will parse up until (and including): c=wlen-1
			for(unsigned int rcoord=wstart;int(rcoord)<=int(wlen)-int(m->len);b++,rcoord++){
				// Check if there is a match, either forward, or in reverse complementary.
				//int match=0;
				if(motifMatchIUPAC(b,wlen-rcoord,m,false)){
					motifOcc*o=occContainer->createMotifOcc(wpos+rcoord,m,false);
					if(!o){
						return false;
					}
				}
				if(motifMatchIUPAC(b,wlen-rcoord,m,true)){
					motifOcc*o=occContainer->createMotifOcc(wpos+rcoord,m,true);
					if(!o){
						return false;
					}
				}
			}
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Motif pairing

double getDistance(motifOcc*a,motifOcc*b){
	if(!a||!b){ cmdWarning("Null-pointer."); return 0; }
	switch(getConfiguration()->distanceMode){
		case dmBetween:{
			long long o1a = a->start;
			long long o1b = o1a+a->mot->len;
			long long o2a = b->start;
			long long o2b = o2a+b->mot->len;
			if(o2a-o1b<o1a-o2b)
				return double(o1a-o2b);
			else
				return double(o2a-o1b);
			}
		case dmCenters:{
			double d_gamma=(double(b->start)+double(b->mot->len)/2.0)
								-(double(a->start)+double(a->mot->len)/2.0);
			if(d_gamma<0)d_gamma=-d_gamma;
			return d_gamma;
			}
		default:return 0;
	}
}

bool overlapping(motifOcc*a,motifOcc*b){
	if(!a||!b){ cmdWarning("Null-pointer."); return false; }
	long long o1a = a->start;
	long long o1b = o1a+a->mot->len;
	long long o2a = b->start;
	long long o2b = o2a+b->mot->len;
	if(o1b < o2a || o2b < o1a)return false;
	else return true;
}

bool isMotifPair(motifOcc*a,motifOcc*b,int cutMin,int cutMax){
	if(!a||!b){ cmdWarning("Null-pointer."); return false; }
	if(!getConfiguration()->motifPairsCanOverlap)if(overlapping(a,b))return false;
	double dist=getDistance(a,b);
	return dist>=cutMin&&dist<=cutMax;
}

