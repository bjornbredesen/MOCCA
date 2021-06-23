////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#include "common.hpp"
#include "vaux.hpp"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

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

////////////////////////////////////////////////////////////////////////////////////
// Console output

bool cmdColorsEnabled = true;

int winColorMap[] = { 8, 12, 10, 14, 9, 5, 11, 15 };

void cmdColor(std::string text, int col) {
	if (cmdColorsEnabled) {
		#ifdef WINDOWS
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
		GetConsoleScreenBufferInfo(hConsole, &csbi);
		WORD oldAttr = csbi.wAttributes;
		SetConsoleTextAttribute(hConsole, winColorMap[col]);
		cout << text;
		SetConsoleTextAttribute(hConsole, oldAttr);
		#elif __EMSCRIPTEN__
		cout << "%<span class=\"tc" << (col&7) << "\"%>" << text << "%</span%>";
		#else
		cout << "\033[0;24;" << int(col+30) << "m" << text << "\033[0m";
		#endif
	} else
		cout << text;
}

void cmdBold(std::string text) {
	if (cmdColorsEnabled) {
		#ifdef WINDOWS
		cout << text;
		#elif __EMSCRIPTEN__
		cout << "%<b%>" << text << "%</b%>";
		#else
		cout << "\033[1;24m" << text << "\033[0m";
		#endif
	} else
		cout << text;
}

void cmdBoldUnderline(std::string text) {
	if (cmdColorsEnabled) {
		#ifdef WINDOWS
		cout << text;
		#elif __EMSCRIPTEN__
		cout << "%<b%>%<u%>" << text << "%</u%>%</b%>";
		#else
		cout << "\033[1;4m" << text << "\033[0m";
		#endif
	} else
		cout << text;
}

void cmdBoldColor(std::string text, int col) {
	if (cmdColorsEnabled) {
		#ifdef WINDOWS
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
		GetConsoleScreenBufferInfo(hConsole, &csbi);
		WORD oldAttr = csbi.wAttributes;
		SetConsoleTextAttribute(hConsole, winColorMap[col]);
		cout << text;
		SetConsoleTextAttribute(hConsole, oldAttr);
		#elif __EMSCRIPTEN__
		cout << "%<span class=\"tc" << (col&7) << "\"%>%<b%>" << text << "%</b%>%</span%>";
		#else
		cout << "\033[1;24;" << int(col+30) << "m" << text << "\033[0m";
		#endif
	} else
		cout << text;
}

void cmdBoldLine(std::string text) {
	if (cmdColorsEnabled) {
		#ifdef WINDOWS
		cout << text;
		#else
		cmdBold(text + std::string("\n"));
		#endif
	} else
		cout << text;
}

void cmdSetColorsEnabled(bool state) {
	cmdColorsEnabled = state;
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

#ifdef __EMSCRIPTEN__
#define EMSCRIPTEN_NCNT 100
int emscripten_sleep_cnt = 0;
#endif

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
	cout << "\r ";
	cmdColor(std::string(cmdTaskC[(cmdTaskI>>5)&3]), cmdColBlue);
	cmdBold(" [ ");
	cmdTask*t=cmdTasks;
	while(t){
		// Print name and status
		if(t->name[0]){
			cmdColor(std::string(t->name), cmdColBlue);
			if(t->text[0])cout << ": " << t->text;
		}else if(t->text[0]){
			cout << t->text;
		}
		t=t->Next;
		if(t)cout << " | ";
	}
	cmdBold(" ]");
	cout << flush;
	#ifdef __EMSCRIPTEN__
	cout << "\n";
	emscripten_sleep_cnt++;
	if (emscripten_sleep_cnt == EMSCRIPTEN_NCNT) {
		emscripten_sleep(0);
		emscripten_sleep_cnt = 0;
	}
	#endif
}

/*
wipe
	Call to wipe out the progress indication.
*/
void cmdTask::wipe(){
	if(cmdTaskSilent)return;
	#ifdef __EMSCRIPTEN__
	cout << "\r";
	emscripten_sleep_cnt++;
	if (emscripten_sleep_cnt == EMSCRIPTEN_NCNT) {
		emscripten_sleep(0);
		emscripten_sleep_cnt = 0;
	}
	#else
	cout << "\r\033[K" << flush;
	#endif
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
	timerStart = time(0);
	name=n?n:(char*)"Unnamed timer";
}

timer::~timer(){
	int sec = time(0) - timerStart;
	int min=int(sec)/60;
	int hr=int(min)/60;
	sec-=min*60;
	min-=hr*60;
	cout << " > ";
	cmdBold("Timer");
	printf("(%s): %d:%02d:%02d\n", name, hr, min, sec);
}

