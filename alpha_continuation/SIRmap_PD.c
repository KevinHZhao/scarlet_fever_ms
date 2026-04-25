/***
Poincare map of the seasonally forced SIR model
to be called from XPPAUT for bifurcation analysis.
Under MacOSX, compile this function via:
gcc −dynamiclib −m32 −o SIRmap.so SIRmap.c
***/
#include <math.h>
/***Compile−time definitions***/
#define Time_step 0.0005 /*in units of years*/
#define TWO_PI 6.283185307179586476925286766559012
#define Real double
#define NDIM 2 /*dimension of the dynamical system*/
/********************************************/
/***FUNCTIONS CALLED BY THE MAIN ROUTINE***/
/********************************************/
/***Euler integrator***/
void Euler(Real *x, Real *dx, Real dt, int ndim)
{
    int i;
    for (i = 0; i < ndim; i++)
    {
        x[i] = x[i] + dx[i] * dt;
    }
}
/***Seasonally forced transmission rate***/
Real Seasonal_beta(Real beta0, Real alpha, Real p, Real time)
{
    // macpan forcing and sinusoidal phase
    Real s1, s2, s3, c1, c2, c3, phase;

    // //Segment 1:
    // s1=-0.0212876629313472;
    // s2=-0.00369216037279813;
    // s3=-0.01727736427829;
    // c1=-0.0663677942156747;
    // c2=0.0354446708440654;
    // c3=-0.0145366333626863;
    // phase=0.528321113831112;
    // //Segment 2:
    // s1=-0.0658573339519255;
    // s2=-0.0112968068359809;
    // s3=0.0110982969184059;
    // c1=-0.0572067656816494;
    // c2=0.00631968906397068;
    // c3=-0.0319569596194759;
    // phase=0.790596156049784;
    // //Segment 3:
    // s1=-0.0454449388736814;
    // s2=0.00943873651368273;
    // s3=0.0326132286775226;
    // c1=-0.0726726084855097;
    // c2=-0.0397831090190311;
    // c3=-0.0120308400694489;
    // phase=0.731643947017079;
    // //Segment 4:
    // s1=-0.012960027196855;
    // s2=0.0304611251271464;
    // s3=-0.00665194190407292;
    // c1=-0.0754377370333189;
    // c2=-0.0638552187924994;
    // c3=0.0250430519014666;
    // phase=0.66708410123755;
    // //Segment 5:
    // s1=-0.00904979587420617;
    // s2=0.0363912092538879;
    // s3=-0.0011500302364437;
    // c1=-0.0947121923563323;
    // c2=-0.0016817215015505;
    // c3=-0.00285495726850529;
    // phase=0.576941610514705;
    // //Segment 6:
    // s1=0.0479184905983673;
    // s2=0.0232695054116134;
    // s3=-0.0158373311943926;
    // c1=-0.0880768617378953;
    // c2=-0.00201508367910589;
    // c3=0.00177814749059823;
    // phase=0.867635921977533;
    //Segment 7:
    s1=0.0685145633625019;
    s2=0.0228342467365722;
    s3=-0.011860788348295;
    c1=-0.0785965823251489;
    c2=0.033112843091219;
    c3=-0.0121522708187906;
    phase=0.506174260007292;
    // //Segment 8:
    // s1=0.058653057872436;
    // s2=0.017215946135317;
    // s3=-0.00783117477331741;
    // c1=-0.0595610990529564;
    // c2=0.0239382970989211;
    // c3=-0.0108237137549129;
    // phase=0.499543148265808;
    // //Segment 9:
    // s1=0.0390955367604266;
    // s2=0.0198214228266663;
    // s3=-0.00199062629266614;
    // c1=-0.0151900612246536;
    // c2=0.00841574338720668;
    // c3=0.0127493097840919;
    // phase=0.295087491703159;
    // //Segment 10:
    // s1=0.0276009980814211;
    // s2=0.0193580086583179;
    // s3=-0.0128026751861481;
    // c1=-0.0218680164031158;
    // c2=0.00754602797337115;
    // c3=0.00504566934889123;
    // phase=0.809798122269995;
    // //Segment 11:
    // s1=-0.0126918014831283;
    // s2=0.0182747739616139;
    // s3=0.000286622628035738;
    // c1=-0.0243480879748356;
    // c2=0.0111615558666263;
    // c3=-0.0028041865339011;
    // phase=0.568416961798229;


    Real mcpn = s1*sin(TWO_PI*time)+s2*sin(2*TWO_PI*time)+s3*sin(3*TWO_PI*time)+c1*cos(TWO_PI*time)+c2*cos(2*TWO_PI*time)+c3*cos(3*TWO_PI*time);

    // Sinusoidal forcing
    Real c2pt; /*cos(2*pi*t)*/
    c2pt = cos(TWO_PI * (time-phase));

    return (beta0 * (1 + alpha * ((1-p)*10*mcpn + p*c2pt)));
}
/************************/
/***THE MAIN ROUTINE***/
/************************/
/***The function SIR map is what XPPAUT calls***/
void SIRmap(Real *in, Real *out, int nin, int nout, Real *var, Real *con)
/*
in=initial and parameter values we get
from the ode file (s,i,R0,alpha,gamma,mu)
out=what we are returning (sp,ip):
calculated values of S and I after one year
nin=dimension of in[]
nout=dimension of out[]
*/
{
    /*define starting values in log base 10*/
    Real s = in[0], i = in[1];
    Real x[NDIM], dx[NDIM]; /*for Euler integrator*/
    /*converting back to the original values,
not in log*/
    s = pow(10, s);
    i = pow(10, i);
    /*define parameter values*/
    Real R0 = in[2], alpha = in[3], gamma = in[4], mu = in[5], p = in[6];
    Real ds, di;
    Real beta0, nonlin_term;
    Real time; /*in units of years*/
    long istep, nsteps;
    /*number of steps in a year*/
    nsteps = (int)(1 / Time_step + 0.5);
    /*integrating for one year*/
    for (istep = 0; istep < nsteps; istep++)
    {
        time = (Real)(istep)*Time_step;
        /*compute the vector field*/
        beta0 = R0 * (gamma + mu); /*mean transmission rate*/
        nonlin_term = Seasonal_beta(beta0, alpha, p, time) * s * i;
        ds = mu - nonlin_term - mu * s;
        di = nonlin_term - (mu + gamma) * i;
        /*integrate using euler's method*/
        x[0] = s;
        x[1] = i;
        dx[0] = ds;
        dx[1] = di;
        Euler(x, dx, Time_step, NDIM);
        s = x[0];
        i = x[1];
    }
    s = log10(s);
    i = log10(i);
    out[0] = s; /*Output in log_10*/
    out[1] = i;
}
