#pragma once
#include "gmp.h"

class CZZn
{
public:
	static void setModulus(char* modulus);
	static void RandInit(char* strSeed);
	static void MultipleInvert(CZZn* nums[], int count);
	static mpz_t T0, T1;

	CZZn(void);
	CZZn(unsigned int data);
	CZZn(char* data);
	CZZn(const CZZn& data);

	~CZZn(void);

	void clone(const CZZn& d)
	{ mpz_set(value,d.value); }
	
	CZZn& operator+=(const CZZn& d)
	{ mpz_add(value, value, d.value); 
	  if (mpz_cmp(value,Modulus)>=0) 
		  mpz_sub(value,value,Modulus);
	  return *this;
	}

	CZZn& operator-=(const CZZn& d)
	{ if (mpz_cmp(value,d.value)<0) 
		mpz_add(value,value,Modulus);
	  mpz_sub(value, value, d.value); 
	  return *this;
	}

	CZZn& operator*=(const CZZn& d);

	CZZn& operator/=(const CZZn& d);
	
	CZZn operator+(const CZZn& d)
	{ return CZZn(*this) += d;} 
	CZZn operator-(const CZZn& d)
	{ return CZZn(*this) -= d;} 
	
	CZZn operator*(const CZZn& d)
	{ return CZZn(*this) *= d;} 
	CZZn operator/(const CZZn& d)
	{ return CZZn(*this) /= d;} 
	
	CZZn& operator=(const CZZn& d)
	{ mpz_set(value, const_cast<CZZn&>(d).value); return *this;	}
	bool operator==(const CZZn& d)
	{  return (mpz_cmp(d.value, value) == 0);}

	int getBit(int pos) {
		return mpz_tstbit(value, pos);
	}
	int getBitCount() {
		return (int)mpz_sizeinbase(value, 2);
	}
	bool isZero()
	{ return (mpz_sgn(value) == 0); }
	void negate()
	{ mpz_sub(value, Modulus, value); }
	void clear()
	{mpz_set_ui(value, 0);	}
	/// Misc..///////////
	int toString(char * buf);

	void Rand();

	int Last32Bits();

private:
	static mpz_t Modulus, M_1, RSquared;
	static int R_Bits;
	static size_t Modulus_BitLen;
	static bool ModulusInited;
	static gmp_randstate_t rand_state;

	CZZn(mpz_t data);

	mpz_t value;
	mpz_t temp;

	void EnterMontgomery(mpz_t& dest, mpz_t& src);
	void MontgomeryRedc(mpz_t& data);

};
