#include "ZZ.h"
#include <string.h>
gmp_randstate_t CZZ::zz_rand_state;


void CZZ::RandInit( char* strSeed )
{
	gmp_randinit_mt(zz_rand_state);
	mpz_t seed;
	mpz_init_set_str(seed, strSeed, 10);
	gmp_randseed(zz_rand_state, seed);
	mpz_clear(seed);
}

CZZ::CZZ( void )
{
	mpz_init(value); mpz_init(t);
}

CZZ::CZZ( int i )
{
	mpz_init(t);mpz_init_set_ui(value, i);
}

CZZ::CZZ( char* d )
{
	mpz_init(t);mpz_init_set_str(value, d, 16);
}

CZZ::CZZ(const CZZ& d )
{
	mpz_init(t);
	mpz_init_set(value ,const_cast<CZZ&>(d).value);
}

CZZ::CZZ( mpz_t& d )
{
	mpz_init(t);mpz_init_set(value ,d);
}
void CZZ::Rand( CZZ& upperbound )
{
	mpz_urandomm(value, zz_rand_state, upperbound.value);
}

int CZZ::compare(CZZ& p )
{
	return mpz_cmp(value, p.value);
}

void CZZ::SetFromInt( int v )
{
	mpz_set_ui(value, v);
}

void CZZ::SetFromStr( char * str )
{
	if (0 != mpz_set_str(value, str, 16))
		exit(999);
}

char * CZZ::toString( void )
{
	int strlen = ((int)mpz_sizeinbase(value, 2) + 7)/ 8 * 2 + 1;
	char * str = (char *)malloc(strlen);
	mpz_set(t, value);
	mpz_get_str(str, -16, t);
	return str;
}

int CZZ::toString( char * buf )
{
	mpz_set(t, value);
	mpz_get_str(buf, -16, t);
	return strlen(buf);
}

bool CZZ::isZero()
{
	return (mpz_sgn(value) ==0);
}

int CZZ::getBitCount()
{
	return (int)mpz_sizeinbase(value, 2);
}

CZZ::~CZZ( void )
{
	mpz_clear(value); 
	mpz_clear(t);
}

int CZZ::getBit( int pos )
{
	return mpz_tstbit(value, pos);
}