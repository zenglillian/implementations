
//  pollard.cpp
//pollard rho method for discrete logarithms


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

void modinv(mpz_t ans,mpz_t a, mpz_t b);
void generator (mpz_t g, mpz_t n);
void fxab(mpz_t a,mpz_t b,mpz_t x, mpz_t p,mpz_t t,mpz_t g);

int main()
{
    mpz_t p,g,t,l,ai,bi,xi,aj,bj,xj,p1,ajsubai,bisubbj,ainv,d,p1d,tmodp,glmodp,m;
    mpz_init(p);
    mpz_init(g);
    mpz_init(t);
    mpz_init(l);
    mpz_init(ai);
    mpz_init(bi);
    mpz_init(xi);
    mpz_init(aj);
    mpz_init(bj);
    mpz_init(xj);
    mpz_init(p1);
    mpz_init(ajsubai);
    mpz_init(bisubbj);
    mpz_init(ainv);
    mpz_init(d);
    mpz_init(p1d);
    mpz_init(tmodp);
    mpz_init(glmodp);
    mpz_init(m);
    mpz_set_str(p,"23",10);
    generator(g,p);
    mpz_set_str(t,"11",10);
    gmp_printf("p = %Zd, g = %Zd, t = %Zd\n",p,g,t);
    mpz_sub_ui(p1,p,1);
    mpz_set_ui(ai,0);
    mpz_set_ui(bi,0);
    mpz_set_ui(xi,1);
    mpz_set_ui(aj,0);
    mpz_set_ui(bj,0);
    mpz_set_ui(xj,1);
    while(mpz_cmp_ui(d,1)<0|| mpz_cmp(xj,xi)!=0){
        fxab(ai,bi,xi,p,t,g);
        fxab(aj,bj,xj,p,t,g);
        fxab(aj,bj,xj,p,t,g);
        mpz_sub(ajsubai,aj,ai);
        mpz_sub(bisubbj,bi,bj);
        mpz_gcd(d,ajsubai,p1);
        gmp_printf("ai=%Zd, bi=%Zd, xi=%Zd, aj=%Zd, bj=%Zd, xj=%Zd, ajsubai=%Zd, d=%Zd \n",ai,bi,xi,aj,bj,xj,ajsubai,d);
        
    }
        if(mpz_cmp_ui(d,1)==0)
        {
            modinv(ainv,ajsubai,p1);
            mpz_mul(l,ainv,bisubbj);
            mpz_mod(l,l,p1);
            gmp_printf("ainv=%Zd \n",ainv);
            gmp_printf("p = %Zd, g = %Zd, t = %Zd, l = %Zd \n",p,g,t,l);
            mpz_clear(p);
            mpz_clear(g);
            mpz_clear(t);
            mpz_clear(l);
            mpz_clear(ai);
            mpz_clear(bi);
            mpz_clear(xi);
            mpz_clear(aj);
            mpz_clear(bj);
            mpz_clear(xj);
            mpz_clear(p1);
            mpz_clear(ajsubai);
            mpz_clear(bisubbj);
            mpz_clear(ainv);
            mpz_clear(d);
            mpz_clear(p1d);
            mpz_clear(tmodp);
            mpz_clear(glmodp);
            mpz_clear(m);
            return 0;
        }
  	    mpz_div(p1d,p1,d);
            modinv(ainv,ajsubai,p1d);
            mpz_mul(l,ainv,bisubbj);
            mpz_mod(l,l,p1d);
   	 mpz_mod(tmodp,t,p);
	gmp_printf("ainv=%Zd, l0=%Zd, p1d=%Zd, tmodp=%Zd \n",ainv,l,p1d,tmodp);
        cout<<"-----------------------"<<endl;
    while(mpz_cmp(m,d)<0 && mpz_cmp(l,p)<0) //since we will always find a answer, do we have to check this?
    {
        mpz_add(l,l,p1d);
        mpz_powm(glmodp,g,l,p);
            gmp_printf("m=%Zd, l=%Zd, glmodp=%Zd \n",m,l,glmodp);
        if(mpz_cmp(glmodp,tmodp)==0)
        {
cout<<"--------------------------------------"<<endl;
            gmp_printf("ainv=%Zd \n",ainv);
            gmp_printf("p = %Zd, g = %Zd, t = %Zd, l = %Zd \n",p,g,t,l);
            mpz_clear(p);
            mpz_clear(g);
            mpz_clear(t);
            mpz_clear(l);
            mpz_clear(ai);
            mpz_clear(bi);
            mpz_clear(xi);
            mpz_clear(aj);
            mpz_clear(bj);
            mpz_clear(xj);
            mpz_clear(p1);
            mpz_clear(ajsubai);
            mpz_clear(bisubbj);
            mpz_clear(ainv);
            mpz_clear(d);
            mpz_clear(p1d);
            mpz_clear(tmodp);
            mpz_clear(glmodp);
                        mpz_clear(m);
            return 0;
        }
mpz_add_ui(m,m,1);
    }

    
}
void modinv(mpz_t ans,mpz_t a, mpz_t b)
{
    mpz_t  a1, a3, b1, b3, t1, t3, q, iter;
    /* Step X1. Initialise */
    mpz_init_set_ui(a1,1);
    mpz_init_set(a3,a);
    mpz_init(b1);
    mpz_init_set(b3,b);
    mpz_init(t1);
    mpz_init(t3);
    mpz_init(q);
    /* Remember odd/even iterations */
    mpz_init_set_ui(iter,1);
    /* Step X2. Loop while b3 != 0 */
    while (mpz_cmp_ui(b3,0) != 0)
    {
        /* Step X3. Divide and "Subtract" */
        mpz_div(q,a3,b3);
        mpz_mod(t3,a3,b3);
        mpz_mul(t1,q,b1);
        mpz_add(t1,t1,a1);//t1 = a1 + q * b1;
        /* Swap */
        mpz_set(a1,b1);
        mpz_set(b1,t1);
        mpz_set(a3,b3);
        mpz_set(b3,t3);
        mpz_neg(iter,iter);
    }
    /* Make sure a3 = gcd(a,b) == 1 */
    if (mpz_cmp_ui(a3,1) != 0)
    {mpz_set_ui(ans,0);}  /* Error: No inverse exists */
    else if (mpz_cmp_ui(iter,0)<0) {mpz_sub(ans,b,a1);}
    else{ mpz_set(ans,a1);}
    mpz_clear(a1);
    mpz_clear(a3);
    mpz_clear(b1);
    mpz_clear(b3);
    mpz_clear(t1);
    mpz_clear(t3);
    mpz_clear(q);
    mpz_clear(iter);
    return;
}

void generator(mpz_t g, mpz_t n)
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
        mpz_add_ui(g,g,1);
        isgenerator = true;
    }
}

void fxab(mpz_t a,mpz_t b,mpz_t x, mpz_t p,mpz_t t,mpz_t g)
{
    mpz_t p1,p2,p3;
    mpz_init(p1);
    mpz_init(p2);
    mpz_init(p3);
    mpz_sub_ui(p1,p,1);
    mpz_div_ui(p2,p,3);
    mpz_mul_ui(p3,p2,2);
    if(mpz_cmp(x,p2)<0)
    {
        mpz_add_ui(a,a,1);
        mpz_mod(a,a,p1);
        mpz_mul(x,x,t);
        mpz_mod(x,x,p);
    }
    else if(mpz_cmp(x,p3)<0)
    {
        mpz_mul_ui(a,a,2);
        mpz_mod(a,a,p1);
        mpz_mul_ui(b,b,2);
        mpz_mod(b,b,p1);
        mpz_mul(x,x,x);
        mpz_mod(x,x,p);
    }
    else
    {
        mpz_add_ui(b,b,1);
        mpz_mod(b,b,p1);
        mpz_mul(x,x,g);
        mpz_mod(x,x,p);
    }
}
