#include "common_gmp.h"

int cellcompare(InputCell* elem1,InputCell* elem2 ){
	int r = elem1->x.compare(elem2->x);
	if (r==0) {
		if (( (elem1)->y) >((elem2)->y) ) r = 1;
		else if (( (elem1)->y) < ((elem2)->y) ) r = -1;
	}
	return r;
}

void innum(CZZ& value ,FILE * f)
{
	char buf[1024];
	if (NULL != fgets(buf, sizeof(buf), f))
	{
		while(buf[strlen(buf)-1] == '\r' ||buf[strlen(buf)-1] == '\n')
			buf[strlen(buf)-1] = '\0';
		value.SetFromStr(buf);
	}
	else 
		value.SetFromInt(0);
}
void otnum(CZZ& value ,FILE * f)
{
	char buf[1024];
	value.toString(buf);
	fprintf(f, "%s\n", buf);
	fflush(f);
}

void MD5_String(char * src, char * dst)
{
	char mapping[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	md5_state_t state;
	char digest[128];
	char c;
	md5_init(&state);
	md5_append(&state, (const md5_byte_t*)src, strlen(src));
	md5_finish(&state, (md5_byte_t*)digest);
	for(int i=15; i>=0; i--)
	{
		c = digest[i];
		digest[2*i+0] = mapping[(c>>4) & 0xF];
		digest[2*i+1] = mapping[(c>>0) & 0xF];
	}
	digest[32] = '\0';
	strcpy(dst, digest);
}
#ifndef WIN32
long GetTickCount()
{
	timeval ts;
	gettimeofday(&ts,0);
	return (long)(ts.tv_sec * 1000 + (ts.tv_usec / 1000));

}
#endif
