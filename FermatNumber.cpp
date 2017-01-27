//Using Pollard Rho method to find a nontrivial factor of a Fermat number n


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

void F(mpz_t n,mpz_t x,long int k);
int main()
{
    int counter = 0;
    long int k = 6;
    long int k2 = pow(2,k);
    mpz_t n,s,g,V,U,UsubV;
    mpz_init(n); //input an composite number n
    mpz_init(s);
    mpz_init(g);
    mpz_init(V);
    mpz_init(U);
    mpz_init(UsubV);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    //Fermat number n
    mpz_ui_pow_ui(n,2,k2);
    mpz_add_ui(n,n,1);
    if(mpz_divisible_ui_p (n,2)!=0){gmp_printf("The non-trivial factor of %Zd is 2 \n",n);}
    do{
        mpz_urandomm (s,state,n);
        mpz_set(U,s);
        mpz_set(V,s);
        do{
            counter++;
            F(n,U,k);
            F(n,V,k);
            F(n,V,k);
            mpz_sub(UsubV, U,V);
            mpz_gcd (g, UsubV, n);
        } while(mpz_cmp_ui(g,1)==0);
    }while(mpz_cmp(g,n)==0);
    
    gmp_printf("The non-trivial factor of %Zd is %Zd \n",n,g);
    cout<<"counter "<<counter<<endl;
    mpz_clear(n);
    mpz_clear(s);
    mpz_clear(g);
    mpz_clear(V);
    mpz_clear(U);
    mpz_clear(UsubV);
    return 0;
    
}

void F(mpz_t n,mpz_t x,long int k)
{
    k=k+2;
    k=pow(2,k);
    mpz_pow_ui(x,x,k);
    mpz_mod(x,x,n);
}