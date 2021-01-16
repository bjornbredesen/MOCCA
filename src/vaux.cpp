////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "vaux.hpp"

/*
cloneString
	Returns a clone of a 0-terminated string, or 0 on failure.
*/
char*cloneString(char*sstr){
	if(!sstr){cmdError("Can't clone null-string.");return 0;}
	size_t sl=strlen(sstr)+1;
	char*r=(char*)malloc(sl);
	if(!r){
		outOfMemory();
		return 0;
	}
	memcpy(r,sstr,sl);
	return r;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Task

// Linked list of currently active tasks and subtasks
cmdTask*cmdTasks=0,*cmdTaskLast=0;

// Task cursor animation
char*cmdTaskC[]={
	(char*)"-",
	(char*)"\\",
	(char*)"|",
	(char*)"/",
};
int cmdTaskI=0;

bool cmdTaskSilent=false;

cmdTask::cmdTask(char*n){
	cmdTaskI=0;
	if(n)strcpy(name,n);
	else name[0]=0;
	text[0]=0;
	Next=0;
	if(cmdTaskLast){
		cmdTaskLast->Next=this;
		Prev=cmdTaskLast;
		cmdTaskLast=this;
	}else{
		cmdTaskLast=cmdTasks=this;
		Prev=0;
	}
	refresh();
}

cmdTask::~cmdTask(){
	if(Prev){
		Prev->Next=0;
		cmdTaskLast=Prev;
	}else{
		cmdTaskLast=cmdTasks=0;
	}
	if(cmdTaskSilent)return;
	wipe();
	refresh();
}

/*
refresh
	Call to redraw.
*/
void cmdTask::refresh(){
	if(!cmdTasks)return;
	if(cmdTaskSilent)return;
	cmdTaskI++;
	if((cmdTaskI>>5)>=4)cmdTaskI=0;
	// Write the status.
	cout << "\r \033[0;34m" << cmdTaskC[(cmdTaskI>>5)&3] << "\033[0m \033[1;24m[\033[0m ";
	cmdTask*t=cmdTasks;
	while(t){
		// Print name and status
		if(t->name[0]){
			cout << "\033[1;34m" << t->name << "\033[0m";
			if(t->text[0])cout << ": " << t->text;
		}else if(t->text[0]){
			cout << t->text;
		}
		t=t->Next;
		if(t)cout << " | ";
	}
	cout << " \033[1;24m]\033[0m" << flush;
}

/*
wipe
	Call to wipe out the progress indication.
*/
void cmdTask::wipe(){
	if(cmdTaskSilent)return;
	cout << "\r\033[K" << flush;
}

// Progress indication methods.
void cmdTask::setText(char*t){
	if(cmdTaskSilent)return;
	strcpy(text,t);
	refresh();
}

void cmdTask::setPercent(double p){
	if(cmdTaskSilent)return;
	sprintf(text,"%lf %%",p);
	refresh();
}

void cmdTask::setInt(int v){
	if(cmdTaskSilent)return;
	sprintf(text,"%d",v);
	refresh();
}

void cmdTask::setIntFraction(int a,int b){
	if(cmdTaskSilent)return;
	sprintf(text,"%d / %d",a,b);
	refresh();
}

void cmdTask::setLong(long l){
	if(cmdTaskSilent)return;
	sprintf(text,"%ld",l);
	refresh();
}

void cmdTask::setLongT(long l,char*t){
	if(cmdTaskSilent)return;
	sprintf(text,"%ld %s",l,t);
	refresh();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Timing

timer::timer(char*n){
	gettimeofday(&timerStart,0);
	name=n?n:(char*)"Unnamed timer";
}

timer::~timer(){
	timeval end;
	gettimeofday(&end,0);
	int sec = int(end.tv_sec-timerStart.tv_sec);
	int min=int(sec)/60;
	int hr=int(min)/60;
	sec-=min*60;
	min-=hr*60;
	printf(" > \033[1;24mTimer\033[0m(%s): %d:%02d:%02d\n", name, hr, min, sec);
}

