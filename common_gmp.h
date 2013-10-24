#ifndef COMMON_GMP
#define COMMON_GMP
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>

#include "ZZ.h"
#include "ECp.h"
#include "md5.h"

#define MAXSIZE 262144
#define DP_BITS 20
#define DP_MASK ((2<<DP_BITS)-1)
#ifdef WIN32
#include "windows.h"
#else
#include <sys/time.h>
long GetTickCount();
#endif

using namespace std;

void ucheck_gmp(char *infile);
void ugenerate_gmp(char *infile);
void check_gmp(char *infile);
void generate_gmp(char *infile);

struct Cell{
	CECPoint p;
	CZZ gi,ki;
};

struct InputCell{
	CZZ x,g,k;
	int y;
};

int cellcompare(InputCell* elem1,InputCell* elem2 );
void innum(CZZ& value ,FILE * f);
void otnum(CZZ& value ,FILE * f);
void MD5_String(char * src, char * dst);
#endif
