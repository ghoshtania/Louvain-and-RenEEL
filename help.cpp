//
//  help.cpp
//  
//
//  Created by Tania Ghosh on 7/26/22.
//

#include <climits>
#include <cstdlib>
#include <ctime>
#include <iostream>
using namespace std;
#include"body.h"
double ran2(long *idum);


int randint(int K,unsigned int *seed)
{
    return ((int)(rand_r(seed)/(RAND_MAX*1.0+1)*K));
    
}

double rand1(unsigned int* seed)
{   //printf(RAND_MAX);
    return (rand_r(seed)/(RAND_MAX*1.0+1));

}

// double absol(double x)
// {    if (x<0)  return (-x); return (x);    }

double comp(struct insert I)
{
    double ans=0;
    int i;
    for (i=0;i<I.com;i++)
    {

        ans+=2*I.wl[i][0]-I.dl[i]*1.0*I.dl[i]/(2*I.n);
    }
    return (ans);
}


void lcopy(int *p, int *q, int N)
{
    int i;
    for (i=0;i<N;i++)
        p[i]=q[i];
}

void erout(FILE *fo,char *message)
{
    fprintf(fo,message);
    exit(0);
}

