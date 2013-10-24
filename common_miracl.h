#define _CRT_SECURE_NO_DEPRECATE
#define _WIN32_WINNT 0x0500
#include "windows.h"
#include "winuser.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#define MAXSIZE 99999
#include "ecn.h"    // Elliptic Curve Class
using namespace std;



struct Cell{
	epoint* p;
	big gi,ki;
};

struct InputCell{
	big x,g,k;
	int y;
} *PInputCell;

int cellcompare(  const void * elem1, const void * elem2 ){
	int r = compare(((InputCell*)elem1)->x,((InputCell*)elem2)->x);
	if (r==0) {
		if (( ((InputCell*)elem1)->y) >(((InputCell*)elem2)->y) ) r = 1;
		else if (( ((InputCell*)elem1)->y) < (((InputCell*)elem2)->y) ) r = -1;
	}
	return r;
}
