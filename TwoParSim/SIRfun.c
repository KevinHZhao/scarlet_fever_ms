#include <math.h>
#include <R.h>
/***Compile-time definitions***/
#define Time_step 0.0005 /*in units of years*/
#define TWO_PI 6.283185307179586476925286766559012
#define Real double
#define NDIM 2 /*dimension of the dynamical system*/

static double parms[4];
#define mu parms[0]
#define gam parms[1]
#define R0 parms[2]
#define a parms[3]
/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=4;
odeparms(&N, parms);
}
/***Seasonally forced transmission rate***/
Real Seasonal_beta(Real beta0, Real time)
{
    Real c2pt; /*cos(2*pi*t)*/
    c2pt = cos(TWO_PI * time);
    return (beta0 * (1 + a * c2pt));
}

void SIRmap(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
    Real beta0, nonlin_term;
    beta0 = R0 * (gam + mu); /*mean transmission rate*/
    nonlin_term = Seasonal_beta(beta0, *t) * y[0] * y[1];
    ydot[0] = mu - nonlin_term - mu * y[0];
    ydot[1] = nonlin_term - (mu + gam) * y[1];
}
