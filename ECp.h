#pragma once
#include "ZZn.h"
#include "ZZ.h"

class CECp
{
public:
	CZZn CoefA,CoefB;
	CECp(void);
public:
	CECp(char * A, char* B);
	~CECp(void);
	bool operator==(CECp& b);
};



class CECPoint
{
	CECp* curve;
public:
	bool infinity;
	CZZn x,y;
	bool isValid();
	CECPoint(void);
	CECPoint(CECp* c);
	CECPoint(CECp* c, char* coord_x, char* coord_y);
	CECPoint(CECp* c, CZZn coord_x, CZZn coord_y);
	CECPoint(const CECPoint& p);
	CECPoint operator+(const CECPoint& p);
	void operator+=(CECPoint& p);
	CECPoint operator*(CZZ factor);
	CECPoint operator*(CZZn factor);
	bool operator==(const CECPoint& p);
	void xToStr(char* Buf) { x.toString(Buf);}
	void yToStr(char* Buf) { y.toString(Buf);}
	bool isInfinity(void){return infinity;}
	~CECPoint(void);
private:
	CZZn Three;
	CZZn Two;
	CZZn tmp_grad, tmp_rx, tmp_ry;
};

#define NUM_MULTIPOINTS 1
class CECMultiPoints
{
	CECPoint points[NUM_MULTIPOINTS];
public:
	CECMultiPoints(CECp* curve);
	void set(const CECPoint& pt, int index);
	CECPoint& get(int index);
	void MultipleAdd(CECPoint* pts[NUM_MULTIPOINTS]);
private:
	CECp* curve;
	CZZn Three;
	CZZn Two;
	CZZn rx[NUM_MULTIPOINTS], ry[NUM_MULTIPOINTS], grad[NUM_MULTIPOINTS];
};
