#ifndef _SIMPLEFUNC_H_
#define _SIMPLEFUNC_H_

#define random(x) (rand()%((x)+1)) // including x
#define randomf(a, b) (((double)rand()/RAND_MAX)*((b)-(a)) + (a))
#define equal(a,b) (fabs((a)-(b))<=FLT_EPSILON)

#endif