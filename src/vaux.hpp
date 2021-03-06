////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Console output

//#define t_indent " \033[35m-\033[0m "
#define t_indent "   "

//----------------------------------------------

void cmdSetColorsEnabled(bool state);

#define cmdColGray 0
#define cmdColRed 1
#define cmdColGreen 2
#define cmdColYellow 3
#define cmdColBlue 4
#define cmdColMag 5
#define cmdColCyan 6
#define cmdColWhite 7

void cmdColor(std::string text, int col);
void cmdBold(std::string text);
void cmdBoldUnderline(std::string text);
void cmdBoldColor(std::string text, int col);
void cmdBoldLine(std::string X);

//----------------------------------------------

#define cmdSepline() cmdColor("-----------------------------------------------------------------\n", cmdColMag);

#define cmdLogo() { cmdBoldColor(" MOCCA", cmdColBlue); cout << "\n Copyright, Bjørn Bredesen, 2013-2021\n bjorn@bjornbredesen.no\n"; }

#define cmdSection(txt) cmdBoldUnderline(txt + std::string("\n"))

#define cmdError(errt) { cout << "\r "; cmdBold("! "); cmdBoldColor("Error", cmdColRed); cout << ": " << errt << "\n"; }

#define cmdWarning(warnt) { cout << "\r "; cmdBold("! "); cmdBoldColor("Warning", cmdColGreen); cout << ": " << warnt << "\n"; }

#define cmdTaskComplete(compt) { cmdTask::wipe(); cmdBoldColor(" * ", cmdColBlue); cmdBold(compt); cout << t_indent; cmdBoldColor("Complete", cmdColGreen); }

#define outOfMemory() cmdError("Out of memory.")

#define argSyntaxError() cmdError("Invalid argument syntax.")

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Div

template<typename T> class autofree{
private:
	autofree(const autofree<T>&other){}
public:
	T*ptr;
	inline autofree(){
		ptr=0;
	}
	inline autofree(T*p){
		ptr=p;
	}
	inline autofree(size_t s){
		ptr=(T*)malloc(sizeof(T)*s);
	}
	inline ~autofree(){
		if(ptr)free(ptr);
	}
	inline T*resize(size_t s){
		ptr=(T*)realloc(ptr,sizeof(T)*s);
		if(!ptr){outOfMemory();return 0;}
		return ptr;
	}
	inline void fill(size_t s,int v){
		if(ptr)memset(ptr,v,sizeof(T)*s);
	}
	inline T*disown(){
		T*r=ptr;
		ptr=0;
		return r;
	}
	inline autofree<T>& operator =(autofree<T>&o){
		if(ptr)free(ptr);
		ptr=o.disown();
		return *this;
	}
	inline T& operator [](int i){
		return ptr[i];
	}
	inline operator T* (){
		return ptr;
	};
};

template<typename T> class autodelete{
private:
	autodelete(const autodelete<T>&other){}
public:
	T*ptr;
	inline autodelete(){
		ptr=0;
	}
	inline autodelete(T*p){
		ptr=p;
	}
	inline ~autodelete(){
		if(ptr)delete ptr;
	}
	inline T*disown(){
		T*r=ptr;
		ptr=0;
		return r;
	}
	inline autodelete<T>& operator =(autodelete<T>&o){
		if(ptr)delete ptr;
		ptr=o.disown();
		return *this;
	}
	inline operator T* (){
		return ptr;
	};
};

template<typename T> class deletevector{
private:
	deletevector(const deletevector<T>&other){}
public:
	std::vector<T*>v;
	inline deletevector(){}
	inline ~deletevector(){
		for(T*e:v)if(e)delete e;
	}
	inline void push_back(T*e){v.push_back(e);};
	inline T* operator [](int i){
		return v[i];
	}
};

char*cloneString(char*sstr);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Task

class cmdTask{
private:
	char name[512];
	char text[512];
	cmdTask*Next,*Prev;
public:
	cmdTask(char*n);
	~cmdTask();
	/*
	refresh
		Call to redraw.
	*/
	static void refresh();
	/*
	wipe
		Call to wipe out the progress indication.
	*/
	static void wipe();
	// Progress indication methods.
	void setText(char*t);
	void setPercent(double p);
	void setInt(int v);
	void setIntFraction(int a,int b);
	void setLong(long l);
	void setLongT(long l,char*t);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Timing

class timer{
private:
	char*name;
	long long timerStart;
public:
	timer(char*n);
	~timer();
};

