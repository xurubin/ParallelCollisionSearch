#pragma once
#include "gmp.h"
#include <stdlib.h>

class CZZ
{
public:
	static void RandInit(char* seed);

	CZZ(void);

	CZZ(int i);
	CZZ(char* d);
	CZZ(const CZZ& d);
	CZZ(mpz_t& d);

	CZZ operator+(const CZZ& p)
	{ mpz_add(t, value, p.value); return CZZ(t); }
	void operator+=(const CZZ& p)
	{ mpz_add(value, value, p.value); }
	CZZ operator-(const CZZ& p)
	{ mpz_sub(t, value, p.value); return CZZ(t); }
	void operator-=(const CZZ& p)
	{ mpz_sub(value, value, p.value); }
	CZZ operator*(const CZZ& p)
	{ mpz_mul(t, value, p.value); return CZZ(t); }
	CZZ operator/(const CZZ& p)
	{ mpz_fdiv_q(t, value, p.value);	return CZZ(t); }
	CZZ operator%(const CZZ& p)	
	{ mpz_fdiv_r(t, value, p.value); 	return CZZ(t); }
	bool operator>(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) > 0;}
	bool operator<(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) < 0;}
	bool operator>=(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) >= 0;}
	bool operator<=(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) <= 0;}
	bool operator==(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) == 0;}
	bool operator!=(const CZZ& p)	
	{ return mpz_cmp(this->value, p.value) != 0;}
	CZZ& operator=(const CZZ& p)
	{	mpz_set(value, const_cast<CZZ&>(p).value); return *this; }

	int getBit(int pos);
	int getBitCount();
	bool isZero();
	int Last32Bits()
	{return value->_mp_d[0];}
	~CZZ(void);

	int toString(char * buf);
	char * toString(void);

	void SetFromStr(char * str);
	void SetFromInt(int v);

	int compare(CZZ& p);

	void Rand(CZZ& upperbound);

public:
	mpz_t value, t;
private:
	static gmp_randstate_t zz_rand_state;


};
