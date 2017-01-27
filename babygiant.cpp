
//baby steps giant steps algorithm

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;


class hash1
{
private:
    size_t tableSize;
    struct item{
        mpz_t gimodp,i;
        item* next;
    };
    item* HashTable;
    
public:
    hash1(mpz_t size);
    size_t getindex(mpz_t key);
    void AddItem(mpz_t gimodp,mpz_t i);
    bool search(mpz_t i,mpz_t key);
    void clear();
};


void modinv(mpz_t ans,mpz_t a, mpz_t b);
void generator (mpz_t g, mpz_t p);

int main()
{
    mpz_t p,b,b1,h,g,t,i,j,gi,thj,gimodp,thjmodp,l;
    mpz_init(p);
    mpz_init(b);
    mpz_init(h);
    mpz_init(g);
    mpz_init(t);
    mpz_init(i);
    mpz_init(j);
    mpz_init(gi);
    mpz_init(thj);
    mpz_init(gimodp);
    mpz_init(thjmodp);
    mpz_init(l);
    mpz_init(b1); //b-1
    
    mpz_set_str(p,"59",10);
    generator (g,p);
    mpz_set_str(t,"52",10);
    //set limits
    mpz_sqrt(b,p);if(mpz_perfect_square_p (p)==0){mpz_add_ui(b,b,1);}
    modinv(h,g,p);
    mpz_powm(h,h,b,p);
    //construct list
    hash1 Hashy(b);
    mpz_set_ui(gi,1);
    mpz_mod(gimodp,gi,p);
    Hashy.AddItem(gimodp,i);
    mpz_sub_ui(b1,b,1);
    while(mpz_cmp(i,b1)<0)
    {
        mpz_add_ui(i,i,1);
        mpz_mul(gi,gi,g);
        mpz_mod(gimodp,gi,p);
        Hashy.AddItem(gimodp,i);
    }
    //sort and find intersection
    mpz_set(thj,t);
    mpz_mod(thjmodp,thj,p);
    
    while(mpz_cmp(j,b)<0)
    {
        if( Hashy.search(i,thjmodp))
        {
            
            mpz_mul(l,j,b);
            mpz_add(l,l,i);
            gmp_printf("p: %Zd, g: %Zd, t: %Zd, l: %Zd \n",p,g,t,l);
            mpz_clear(p);
            mpz_clear(b);
            mpz_clear(h);
            mpz_clear(g);
            mpz_clear(t);
            mpz_clear(i);
            mpz_clear(j);
            mpz_clear(gi);
            mpz_clear(thj);
            mpz_clear(gimodp);
            mpz_clear(thjmodp);
            mpz_clear(l);
            mpz_clear(b1);
            Hashy.clear();
            return 0;
        }
        mpz_add_ui(j,j,1);
        mpz_mul(thj,thj,h);
        mpz_mod(thjmodp,thj,p);
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
        mpz_add_ui(g,g,1);
        isgenerator = true;
    }
}


hash1::hash1(mpz_t size)
{ //initialize the array
    size_t k;
    tableSize = mpz_get_ui(size);
    HashTable = new item[tableSize];
    for( k = 0;k<tableSize;k++)
    {
        mpz_init_set_ui(HashTable[k].gimodp,-1);
        mpz_init_set_ui(HashTable[k].i,-1);
        HashTable[k].next = NULL;
    }
    
}

size_t hash1::getindex(mpz_t key)
{
    size_t index;
    mpz_t idx;
    mpz_init(idx);
    mpz_mod_ui(idx,key,tableSize);
    index = mpz_get_ui (idx);
    return index;
}

void hash1::AddItem(mpz_t gimodp,mpz_t i)
{
    size_t index = getindex(gimodp);
    if(mpz_cmp_ui(HashTable[index].gimodp,-1)==0)
    {
        mpz_init_set(HashTable[index].gimodp,gimodp);
        mpz_init_set(HashTable[index].i,i);
    }
    else
    {
        item* Ptr;
        Ptr= &HashTable[index];
        item* p;
        p = new item;
        mpz_init_set(p->gimodp,gimodp);
        mpz_init_set(p->i,i);
        p->next = NULL;
        while(Ptr->next != NULL)
        {
            Ptr = Ptr->next;
        }
        Ptr->next = p;
    }
    return;
}

bool hash1::search(mpz_t i, mpz_t key)
{
    size_t index = getindex(key);
    item* Ptr;
    Ptr= &HashTable[index];
    while(Ptr != NULL)
    {
        if(mpz_cmp(Ptr->gimodp,key)==0)
        {
            mpz_set(i,Ptr->i);
            return true;
        }
        Ptr = Ptr->next;
    }
    return false;
}

void hash1::clear()
{
    delete []HashTable;
}












