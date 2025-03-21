#ifndef _XIS
#define _XIS

/*
  Generating schemes for different types of random variables
    -> +/-1 random variables
    -> k-valued random variables
  Fast range-summable random variables
  For details see the papers:
	1) "Fast Range-Summable Random Variables for Efficient Aggregate Estimation" by F. Rusu and A. Dobra
	2) "Pseudo-Random Number Generation for Sketch-Based Estimations" by F. Rusu and A. Dobra
*/

#include "gen_scheme.h"
#include "range_sum.h"
#include "RM7_range_sum.h"

using namespace std;



class Xi
{
  public:
    virtual double element(unsigned int j) = 0;
    virtual double interval_sum(unsigned int alpha, unsigned int beta) = 0;

    virtual ~Xi();
};


/*
+/-1 random variables that are 3-wise independent
*/

class Xi_BCH3 : public Xi
{
  protected:
    unsigned int seeds[2];

  public:
    Xi_BCH3(unsigned int I1, unsigned int I2);
    virtual ~Xi_BCH3();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



class Xi_EH3 : public Xi
{
  protected:
    unsigned int seeds[2];

  public:
    Xi_EH3(unsigned int I1, unsigned int I2);
    virtual ~Xi_EH3();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



class Xi_CW2 : public Xi
{
  protected:
    unsigned int seeds[2];

  public:
    Xi_CW2(unsigned int I1, unsigned int I2);
    virtual ~Xi_CW2();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



/*
B-valued random variables that are 2-wise independent
*/

class Xi_CW2B : public Xi
{
  protected:
    unsigned int seeds[2];
    unsigned int buckets_no;

  public:
    Xi_CW2B(unsigned int I1, unsigned int I2, unsigned int B);
    virtual ~Xi_CW2B();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



/*
+/-1 random variables that are 4-wise independent
*/

class Xi_CW4 : public Xi
{
  protected:
    unsigned int seeds[4];

  public:
    Xi_CW4(unsigned int I1, unsigned int I2);
    virtual ~Xi_CW4();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



/*
B-valued random variables that are 4-wise independent
*/

class Xi_CW4B : public Xi
{
  protected:
    unsigned int seeds[4];
    unsigned int buckets_no;

  public:
    Xi_CW4B(unsigned int I1, unsigned int I2, unsigned int B);
    virtual ~Xi_CW4B();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



/*
+/-1 random variables that are 5-wise independent
*/

class Xi_BCH5 : public Xi
{
  protected:
    unsigned int seeds[3];

  public:
    Xi_BCH5(unsigned int I1, unsigned int I2);
    virtual ~Xi_BCH5();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



class Xi_RM7 : public Xi
{
  protected:
    unsigned int seeds[33];

  public:
    Xi_RM7(unsigned int I1, unsigned int I2);
    virtual ~Xi_RM7();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};





/*
Random variables for dyadic mapping as an alternative to fast range-summation
*/

class Xi_Dyadic_Map_EH3 : public Xi
{
  protected:
    unsigned int seeds[2];
    unsigned int dom_size;

  public:
    Xi_Dyadic_Map_EH3(unsigned int Dom_size, unsigned int I1, unsigned int I2);
    virtual ~Xi_Dyadic_Map_EH3();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};



class Xi_Dyadic_Map_BCH5 : public Xi
{
  protected:
    unsigned int seeds[3];
    unsigned int dom_size;

  public:
    Xi_Dyadic_Map_BCH5(unsigned int Dom_size, unsigned int I1, unsigned int I2);
    virtual ~Xi_Dyadic_Map_BCH5();

    virtual double element(unsigned int j);
    virtual double interval_sum(unsigned int alpha, unsigned int beta);
};


Xi::~Xi()
{
}




Xi_BCH3::Xi_BCH3(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_BCH3::~Xi_BCH3()
{
}


double Xi_BCH3::element(unsigned int j)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  double res = BCH3(i0, i1, j);
  return res;
}


double Xi_BCH3::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  double res = BCH3_Range(alpha, beta, i1, i0);
  return res;
}




Xi_EH3::Xi_EH3(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_EH3::~Xi_EH3()
{
}


double Xi_EH3::element(unsigned int j)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  double res = EH3(i0, i1, j);
  return res;
}


double Xi_EH3::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  double res = EH3_Range(alpha, beta, i1, i0);
  return res;
}




Xi_CW2::Xi_CW2(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_CW2::~Xi_CW2()
{
}


double Xi_CW2::element(unsigned int j)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];

  double res = CW2(a, b, j);
  return res;
}


double Xi_CW2::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];

  double res = 0;

  for (unsigned int k = alpha; k <= beta; k++)
    res += CW2(a, b, k);

  return res;
}





Xi_CW2B::Xi_CW2B(unsigned int I1, unsigned int I2, unsigned int B)
{
  buckets_no = B;

  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_CW2B::~Xi_CW2B()
{
}


double Xi_CW2B::element(unsigned int j)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];

  double res = CW2B(a, b, j, buckets_no);
  return res;
}


double Xi_CW2B::interval_sum(unsigned int alpha, unsigned int beta)
{
  return -1;
}







Xi_BCH5::Xi_BCH5(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[2] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_BCH5::~Xi_BCH5()
{
}


double Xi_BCH5::element(unsigned int j)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];
  unsigned int i2 = seeds[2];

  double res = BCH5(i0, i1, i2, j);
  return res;
}


double Xi_BCH5::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];
  unsigned int i2 = seeds[2];

  double res = 0;

  for (unsigned int k = alpha; k <= beta; k++)
    res += BCH5(i0, i1, i2, k);

  return res;
}




Xi_RM7::Xi_RM7(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  for (int i = 0; i < 33; i++)
  {
    seeds[i] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

    I1 = 36969*(I1 & 0177777) + (I1>>16);
    I2 = 18000*(I2 & 0177777) + (I2>>16);
  }
}


Xi_RM7::~Xi_RM7()
{
}


double Xi_RM7::element(unsigned int j)
{
  unsigned int i0 = seeds[32];
  unsigned int I[32];
  for (int k = 0; k < 32; k++)
    I[k] = seeds[k];

  double res = RM7(i0, I, j);
  return res;
}


double Xi_RM7::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[32];
  unsigned int I[32];
  for (int k = 0; k < 32; k++)
    I[k] = seeds[k];

  double res = Interval(alpha, beta, I, i0);
  return res;
}




Xi_CW4::Xi_CW4(unsigned int I1, unsigned int I2)
{
  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[2] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[3] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_CW4::~Xi_CW4()
{
}


double Xi_CW4::element(unsigned int j)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];
  unsigned long c = seeds[2];
  unsigned long d = seeds[3];

  double res = CW4(a, b, c, d, j);
  return res;
}


double Xi_CW4::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];
  unsigned long c = seeds[2];
  unsigned long d = seeds[3];

  double res = 0;

  for (unsigned int k = alpha; k <= beta; k++)
    res += CW4(a, b, c, d, k);

  return res;
}




Xi_CW4B::Xi_CW4B(unsigned int I1, unsigned int I2, unsigned int B)
{
  buckets_no = B;

  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[2] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[3] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_CW4B::~Xi_CW4B()
{
}


double Xi_CW4B::element(unsigned int j)
{
  unsigned long a = seeds[0];
  unsigned long b = seeds[1];
  unsigned long c = seeds[2];
  unsigned long d = seeds[3];

  double res = CW4B(a, b, c, d, j, buckets_no);
  return res;
}


double Xi_CW4B::interval_sum(unsigned int alpha, unsigned int beta)
{
  return -1;
}




Xi_Dyadic_Map_EH3::Xi_Dyadic_Map_EH3(unsigned int Dom_size, unsigned int I1, unsigned int I2)
{
  dom_size = Dom_size;

  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_Dyadic_Map_EH3::~Xi_Dyadic_Map_EH3()
{
}


double Xi_Dyadic_Map_EH3::element(unsigned int j)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  unsigned int front_mask = 0;
  unsigned int map;

  double res = 0;

  for (int k = 0; k < dom_size; k++)
  {
    map = front_mask ^ (j >> k);
    res += EH3(i0, i1, map);

    if ((j >> k) & 1U == 1U)
      j ^= (1 << k);

    front_mask = (front_mask >> 1) ^ (1U << dom_size);
  }

  res += EH3(i0, i1, front_mask);
  return res;
}


double Xi_Dyadic_Map_EH3::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];

  double res = Dyadic_Range_EH3(alpha, beta, i0, i1, dom_size);
  return res;
}




Xi_Dyadic_Map_BCH5::Xi_Dyadic_Map_BCH5(unsigned int Dom_size, unsigned int I1, unsigned int I2)
{
  dom_size = Dom_size;

  const unsigned int k_mask = 0xffffffff;

  seeds[0] = ((I1 << 16)^(I2 & 0177777)) & 1UL;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[1] = ((I1 << 16)^(I2 & 0177777)) & k_mask;

  I1 = 36969*(I1 & 0177777) + (I1>>16);
  I2 = 18000*(I2 & 0177777) + (I2>>16);

  seeds[2] = ((I1 << 16)^(I2 & 0177777)) & k_mask;
}


Xi_Dyadic_Map_BCH5::~Xi_Dyadic_Map_BCH5()
{
}


double Xi_Dyadic_Map_BCH5::element(unsigned int j)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];
  unsigned int i2 = seeds[2];

  unsigned int front_mask = 0;
  unsigned int map;

  double res = 0;
  
  for (int k = 0; k < dom_size; k++)
  {
    map = front_mask ^ (j >> k);
    res += BCH5(i0, i1, i2, map);

    if ((j >> k) & 1U == 1U)
      j ^= (1 << k);

    front_mask = (front_mask >> 1) ^ (1U << dom_size);
  }

  res += BCH5(i0, i1, i2, front_mask);
  return res;
}


double Xi_Dyadic_Map_BCH5::interval_sum(unsigned int alpha, unsigned int beta)
{
  unsigned int i0 = seeds[0];
  unsigned int i1 = seeds[1];
  unsigned int i2 = seeds[2];

  double res = Dyadic_Range_BCH5(alpha, beta, i0, i1, i2, dom_size);
  return res;
}




#endif
