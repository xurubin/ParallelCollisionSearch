#define USE_GMP

#ifndef USE_GMP
#include "common_miracl.h"
#include "unknownorder_miracl.h"
#include "knownorder_miracl.h"
#else
#include "common_gmp.h"
#define generate(x) generate_gmp(x)
#define check(x) check_gmp(x)
#define ugenerate(x) ugenerate_gmp(x)
#define ucheck(x) ucheck_gmp(x)
#endif
#ifdef WIN32
#include "windows.h"
#endif

int main(int argc,char **argv)
{
    argv++; argc--;
	if (argc!=2) {
		printf("Wrong number of arguments\r\n");
		printf("-g(G) *.ecs or -c *.dat\r\n");
		printf("For unknown orders, use -gu(GU) *.ecs or -cu *.dat\r\n");
		return 1;
	}
	if (strcmp(argv[0],"-g")==0) generate(argv[1]);
	else if (strcmp(argv[0],"-G")==0) {
#ifdef WIN32
		ShowWindow(FindWindowA("ConsoleWindowClass",0),SW_HIDE);
#endif
		generate(argv[1]);
	}
	else if (strcmp(argv[0],"-c")==0) check(argv[1]);

	else if (strcmp(argv[0],"-gu")==0) ugenerate(argv[1]);
	else if (strcmp(argv[0],"-GU")==0) {
#ifdef WIN32
		ShowWindow(FindWindowA("ConsoleWindowClass",0),SW_HIDE);
#endif
		ugenerate(argv[1]);
	}
	else if (strcmp(argv[0],"-cu")==0) ucheck(argv[1]);

}
