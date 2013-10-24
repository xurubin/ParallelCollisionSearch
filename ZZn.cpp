
#include "ZZn.h"
#include <string.h>
//#define Montgomery_Redc

/////////////////Global functions //////////////
void CZZn::RandInit(char* strSeed)
{
	gmp_randinit_mt(rand_state);
	mpz_t seed;
	mpz_init_set_str(seed, strSeed, 10);
	gmp_randseed(rand_state, seed);
	mpz_clear(seed);
}

void CZZn::setModulus(char* modulus)
{
	if (!ModulusInited) 
	{
		mpz_init(Modulus);
		mpz_init(M_1);
		mpz_init(RSquared);
		ModulusInited = true;
	}
	mpz_set_str(Modulus,modulus,16);
	Modulus_BitLen = mpz_sizeinbase(Modulus,2);

	R_Bits = ( ((Modulus_BitLen - 1)>>5) + 1) << 5; //Nearest DWORDs
	mpz_t R; mpz_init_set_ui(R,1);
	mpz_mul_2exp(R, R, R_Bits);
	mpz_mul(RSquared,R,R); mpz_mod(RSquared,RSquared,Modulus);
	mpz_set(M_1, Modulus);
	mpz_invert(M_1, M_1, R);
	mpz_sub(M_1, R, M_1);

	/*
	mpz_t t;
	char buf[200];
	mpz_init_set(t, M_1);
	mpz_get_str(buf, -16, t);
	*/
}
mpz_t CZZn::M_1;
mpz_t CZZn::Modulus;
mpz_t CZZn::RSquared;
int CZZn::R_Bits = 0;
bool CZZn::ModulusInited = false;
size_t CZZn::Modulus_BitLen = 0;
gmp_randstate_t CZZn::rand_state;
//////////////////Constructors/Destructor ////////
CZZn::CZZn(void)
{
	mpz_init(value);
	mpz_init2(temp, Modulus_BitLen*2);
}
CZZn::CZZn(unsigned int data)
{
	mpz_init_set_ui(value,data);
#ifdef Montgomery_Redc
	EnterMontgomery(value,value);
#endif
	mpz_init2(temp, Modulus_BitLen*2);
}
CZZn::CZZn(char* data)
{
	mpz_init_set_str(value,data,16);
#ifdef Montgomery_Redc
	EnterMontgomery(value,value);
#endif
	mpz_init2(temp, Modulus_BitLen*2);
}
CZZn::CZZn(const CZZn& data)
{
	mpz_init_set(value,const_cast<CZZn&>(data).value);
	mpz_init2(temp, Modulus_BitLen*2);
}
CZZn::CZZn(mpz_t data)
{
	mpz_init_set(value,data);
	mpz_init2(temp, Modulus_BitLen*2);
}
CZZn::~CZZn(void)
{
	mpz_clear(temp);
	mpz_clear(value);
}

/////////////Arithmetics Operations/////

void CZZn::EnterMontgomery(mpz_t& dest, mpz_t& src)
{
	mpz_mul_2exp(dest,src, R_Bits);
	mpz_mod(dest, dest, Modulus);
}

void CZZn::MontgomeryRedc(mpz_t& data)
{
	mpz_mul(temp, data, M_1);
	mpz_tdiv_r_2exp( temp, temp, R_Bits);
	mpz_addmul(data, temp, Modulus);
	mpz_tdiv_q_2exp(data, data, R_Bits);
	if (mpz_cmpabs(data, Modulus)>=0)
		mpz_sub(data, data, Modulus);
}

CZZn& CZZn::operator*=(const CZZn& d)
{ 
#ifdef Montgomery_Redc
	mpz_mul(value, value, d.value);
	MontgomeryRedc(value);
#else
	mpz_mul(temp, value, d.value); 
	mpz_mod(value,temp,Modulus); 
#endif
	return *this;
}

CZZn& CZZn::operator/=(const CZZn& d2)
{ 
#ifdef Montgomery_Redc
	mpz_set(temp, d.value);
	MontgomeryRedc(temp);
	mpz_invert(temp, temp, Modulus);
	mpz_mul(value, temp, value);
	MontgomeryRedc(value);
#else
	CZZn& d = const_cast<CZZn&>(d2);
	mpz_invert(temp, d.value, Modulus);
	mpz_mul(temp, temp, value);
	mpz_mod(value,temp,Modulus); 
#endif
	return *this;
}

/////////////Misc.../////////////////
int CZZn::toString(char * buf)
{
	mpz_t t;
	mpz_init_set(t, value);
#ifdef Montgomery_Redc
	MontgomeryRedc(t);	
#endif	
	mpz_get_str(buf, -16, t);
	mpz_clear(t);
	return strlen(buf);
}

int CZZn::Last32Bits()
{
	return value->_mp_d[0];
}

void CZZn::Rand()
{
	mpz_urandomm(value, rand_state, Modulus);
}

mpz_t CZZn::T0;
mpz_t CZZn::T1;
void CZZn::MultipleInvert(CZZn* nums[], int count )
{
	if (count == 1)
	{
		mpz_invert(nums[0]->value,nums[0]->value,Modulus);
		return;
	}
	else if( count == 2)
	{
		mpz_t& tmp0 = nums[0]->temp;
		mpz_t& tmp1 = nums[1]->temp;

		mpz_mul(tmp0, nums[0]->value, nums[1]->value);
		mpz_mod(tmp1, tmp0, Modulus);

		mpz_invert(tmp0, tmp1, Modulus);

		mpz_mul(tmp1, tmp0, nums[0]->value);
		mpz_mul(tmp0, tmp0, nums[1]->value);
		mpz_mod(nums[0]->value, tmp0, Modulus);
		mpz_mod(nums[1]->value, tmp1, Modulus);
		return;
	}

	if (T0 == 0)
	{
		mpz_init(T0);
		mpz_init(T1);
	}

	mpz_t& tmp0 = nums[0]->temp;
	for(int i=1;i<count;i++)
	{
		mpz_mul(tmp0, nums[i]->value, (i == 1) ? nums[i-1]->value : nums[i-1]->temp);
		mpz_mod(nums[i]->temp, tmp0, Modulus);
	}
	mpz_invert(tmp0, nums[count-1]->temp, Modulus);
	mpz_t& tmp1 = nums[count-1]->temp;
	mpz_t& tmp2 = nums[count-2]->temp;
	for(int i=count-1;i>=1;i--)
	{
		mpz_mul(tmp1, tmp0, (i == 1) ? nums[i-1]->value : nums[i-1]->temp);
		mpz_mul(tmp2, tmp0, nums[i]->value); //tmp2 only available from this line onwards
		mpz_mod(tmp0, tmp2, Modulus);
		mpz_mod(nums[i]->value, tmp1, Modulus);
	}
	mpz_set(nums[0]->value, tmp0);

/*
//Recover two inverses, 
//before the invokation tmp##x should have (ab)-1
//after this nums[x] = a^-1 and nums[y] = b^-1
#define RECOVER_TWO(x,y,t) 		mpz_mul(tmp##y, tmp##x, nums[x]->value); \
							mpz_mul(tmp##t, tmp##x, nums[y]->value); \
							mpz_mod(nums[x]->value, tmp##t, Modulus); \
							mpz_mod(nums[y]->value, tmp##y, Modulus);

	int i;
	if (T0 == 0)
	{
		mpz_init(T0);
		mpz_init(T1);
	}
	mpz_t& tmp0 = nums[0]->temp;
	mpz_t& tmp1 = nums[1]->temp;

	mpz_set(tmp1, nums[0]->value);
	for(i=1;i<count;i++)
	{
		mpz_mul(tmp0, tmp1, nums[i]->value);
		mpz_mod(tmp1, tmp0, Modulus);
	}
	mpz_invert(tmp0, tmp1, Modulus);
	//Now tmp0 = (num[0]*num[1]*..num[count-1])^(-1)
	if (count == 2)
	{
		RECOVER_TWO(0, 1, 0);
	}
	else if (count == 3)
	{
		mpz_t& tmp2 = nums[2]->temp;
		mpz_mul(tmp2, tmp0, nums[0]->value);
		mpz_mod(tmp1, tmp2, Modulus);

		mpz_mul(nums[0]->value, nums[2]->value, tmp0);
		mpz_mod(tmp0, nums[0]->value, Modulus);

		mpz_mul(tmp2, nums[1]->value, tmp0);
		mpz_mod(nums[0]->value, tmp2, Modulus);

		RECOVER_TWO(1, 2, 0);
	}
	else if (count == 4)
	{
		mpz_t& tmp2 = nums[2]->temp;
		mpz_t& tmp3 = nums[3]->temp;

		mpz_mul(tmp3, tmp0, nums[0]->value);
		mpz_mod(tmp2, tmp3, Modulus);
		mpz_mul(tmp3, tmp2, nums[1]->value);
		mpz_mod(tmp2, tmp3, Modulus);
		
		mpz_mul(tmp3, tmp0, nums[2]->value);
		mpz_mod(tmp0, tmp3, Modulus);
		mpz_mul(tmp3, tmp0, nums[3]->value);
		mpz_mod(tmp0, tmp3, Modulus);

		RECOVER_TWO(0, 1, 3);
		RECOVER_TWO(2, 3, 1);
	}
	else if (count == 6)
	{
		mpz_t& tmp2 = nums[2]->temp;
		mpz_t& tmp3 = nums[3]->temp;
		mpz_t& tmp4 = nums[4]->temp;
		mpz_t& tmp5 = nums[5]->temp;
		
		//nums[5]  = ef
		mpz_mul(tmp4, nums[4]->value, nums[5]->value);
		mpz_mod(tmp5, tmp4, Modulus);
		//nums[3]  = cd
		mpz_mul(tmp4, nums[2]->value, nums[3]->value);
		mpz_mod(tmp3, tmp4, Modulus);
		//nums[1]  = ab
		mpz_mul(tmp4, nums[0]->value, nums[1]->value);
		mpz_mod(tmp1, tmp4, Modulus);
		
		//nums[4] = ab*cd*tmp0 = ef ^ -1
		//nums[2] = ab*ef*tmp0 = cd ^ -1
		mpz_mul(tmp2, tmp0, tmp1);
		mpz_mod(T0, tmp2, Modulus); //T0 = tmp0*ab
		mpz_mul(T1, T0, tmp3);
		mpz_mod(tmp4, T1, Modulus); //tmp4 = tmp0*ab*cd
		mpz_mul(T1, T0, tmp5);
		mpz_mod(tmp2, T1, Modulus); //tmp2 = tmp0*ab*ef

		//nums[0] = cd*ef*tmp0 = ab ^ -1
		mpz_mul(tmp1, tmp0, tmp3);
		mpz_mod(tmp0, tmp1, Modulus);
		mpz_mul(tmp1, tmp0, tmp5);
		mpz_mod(tmp0, tmp1, Modulus);

		RECOVER_TWO(0, 1, 3);
		RECOVER_TWO(2, 3, 1);
		RECOVER_TWO(4, 5, 3);
	}
*/
}

