#include "f2c.h"
#include "config.h"

#define log10e 0.43429448190325182765

#ifdef KR_headers
double log();
double igraphd_lg10(x) doublereal *x;
#else
#undef abs
#include "math.h"
double igraphd_lg10(doublereal *x)
#endif
{
return( log10e * log(*x) );
}
