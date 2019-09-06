////////////////////////////////////////////////////////////////////////////////////
// MOCCA
// Copyright, Bjørn Bredesen, 2019
// E-mail: bjorn@bjornbredesen.no
////////////////////////////////////////////////////////////////////////////////////
// General

#pragma once

////////////////////////////////////////////////////////////////////////////////////
// Console output

// Some general output styling
#define sepline "\033[35m-----------------------------------------------------------------\033[0m\n"
#define t_success "\033[1;32mSuccess\033[0m"
#define m_error "\r \033[1;24m!\033[0m \033[1;31mError\033[0m: "
#define m_warning "\r \033[1;24m!\033[0m \033[1;32mWarning\033[0m: "

#define t_indent " \033[35m-\033[0m "
#define t_star " \033[1;34m*\033[0m "
#define t_bold(X) " \033[1;24m" << X << "\033[0m "


#define cmdSection(txt) cout << "\033[1;4m" << txt << "\033[0m:\n"

#define cmdError(errt) cout << m_error << errt << "\n"

#define cmdWarning(warnt) cout << m_warning << warnt << "\n"

#define cmdTaskComplete(compt) cmdTask::wipe();cout << t_star << "\033[1;24m" << compt << "\033[0m" << t_indent << "\033[1;32mComplete\033[0m\n"

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
	timeval timerStart;
public:
	timer(char*n);
	~timer();
};

