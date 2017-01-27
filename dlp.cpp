//the naive algrithm for discrete logrithm

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

void generator (mpz_t g, mpz_t n);

int main()
{
    mpz_t n,g,t,l,gl,glmodn,tmodn;
    mpz_init(n);
    mpz_init(g);
    mpz_init(t);
    mpz_init(gl);
    mpz_init(glmodn);
    mpz_init(tmodn);
    mpz_init(l);
    
    mpz_set_str(n,"59",10);
    generator (g,n);
    mpz_set_ui(gl,1);
    //mpz_set_str(t,"19",10);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    mpz_urandomm (t,state,n); if(mpz_cmp_ui(t,0)==0){mpz_add_ui(t,t,1);}//randomly choose t from 1 to n-1
    
    mpz_mod(tmodn,t,n);
    while(mpz_cmp(l,n)<0)
    {
        if(mpz_cmp(glmodn,tmodn)==0)
        {
            gmp_printf("n = %Zd, g = %Zd, t = %Zd, l = %Zd \n",n,g,t,l);
            mpz_clear(n);
            mpz_clear(g);
            mpz_clear(t);
            mpz_clear(gl);
            mpz_clear(glmodn);
            mpz_clear(tmodn);
            mpz_clear(l);
            return 0;
        }
        mpz_add_ui(l,l,1);
        mpz_mul(gl,gl,g);
        mpz_mod(glmodn,gl,n);
    }
}


void generator (mpz_t g, mpz_t n)
{
    mpz_set_ui(g,2);
    bool isgenerator =  true;
    mpz_t pm,n1const,n1,n1i,i,sqrtN1;
    mpz_init(pm);
    mpz_init(n1const); //n1 = n-1
    mpz_init(n1);
    mpz_init(n1i); //n1i = (n-1)/i
    mpz_init(i);
    mpz_init(sqrtN1);
    mpz_sub_ui(n1,n,1);
    mpz_sub_ui(n1const,n,1);
    mpz_sqrt(sqrtN1,n1);
    while( mpz_cmp(g,n)<0)
    {
        mpz_set_ui (i ,2);
        mpz_set(n1,n1const);
        while(mpz_cmp_ui(n1,1)>0&&isgenerator)
        {
            int count=0;
            while(mpz_divisible_p (n1,i)!=0&&isgenerator)
            {
                count++;
                if(count==1)
                {
                    mpz_div(n1i,n1const,i);
                    mpz_powm(pm,g,n1i,n);
                    if (mpz_cmp_ui(pm,1)==0){isgenerator = false;}
                }
                mpz_div(n1,n1,i);
            }
            mpz_nextprime(i,i);
        }
        if(isgenerator)
        {
            mpz_clear(pm);
            mpz_clear(n1);
            mpz_clear(n1const);
            mpz_clear(n1i);
            mpz_clear(i);
            mpz_clear(sqrtN1);
            return;
        }
        mpz_nextprime(g,g);
        isgenerator = true;
    }
    
}
