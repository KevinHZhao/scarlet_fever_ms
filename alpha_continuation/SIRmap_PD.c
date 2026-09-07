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

  // // Segment 1:
  // s1=-0.0172596047072509;
  // s2=-0.0032750035825358;
  // s3=-0.0169163749764281;
  // c1=-0.0661042367339843;
  // c2=0.0363710448520568;
  // c3=-0.0167933204037997;
  // phase=0.525133470780461;
  // // Segment 2:
  // s1=-0.0702751736778501;
  // s2=-0.0121177716268683;
  // s3=0.0112020114679481;
  // c1=-0.0588818800195153;
  // c2=0.00684702224883206;
  // c3=-0.0323928115154836;
  // phase=0.790566891905994;
  // // Segment 3:
  // s1=-0.0450813495559104;
  // s2=0.0097149537350211;
  // s3=0.0334527705884202;
  // c1=-0.0740105493557181;
  // c2=-0.0417912858023739;
  // c3=-0.0136868820352668;
  // phase=0.733182753652705;
  // // Segment 4:
  // s1=-0.0124211931654728;
  // s2=0.0312170002868357;
  // s3=-0.00632325889777525;
  // c1=-0.0772566125388097;
  // c2=-0.0644040065731635;
  // c3=0.0248612206933521;
  // phase=0.666874670438105;
  // // Segment 5:
  // s1=-0.00962744998084198;
  // s2=0.0361226554429802;
  // s3=-0.00203455079807741;
  // c1=-0.0946602192024968;
  // c2=-0.00263802308542256;
  // c3=-0.00152774994463603;
  // phase=0.581288668752062;
  // // Segment 6:
  // s1=0.0480282963776851;
  // s2=0.0230387761735651;
  // s3=-0.01534218810509;
  // c1=-0.0863627253832895;
  // c2=-0.00255487419017575;
  // c3=0.000923555973662662;
  // phase=0.871312223108843;
  // Segment 7:
  s1=0.0697647253272409;
  s2=0.0226416278979439;
  s3=-0.0127066665901618;
  c1=-0.0784230260201418;
  c2=0.0334556348686397;
  c3=-0.0126583474838829;
  phase=0.506482686864681;
  // // Segment 8:
  // s1=0.072981794200076;
  // s2=0.0365097392605001;
  // s3=-0.0155751191451615;
  // c1=-0.0598445655529398;
  // c2=0.0194869730606118;
  // c3=-0.0111536072520092;
  // phase=0.210657597857167;
  // // Segment 9:
  // s1=0.0562395842901831;
  // s2=0.00325352360681964;
  // s3=-0.00186115006397254;
  // c1=-0.0656351897010105;
  // c2=0.0286350769973658;
  // c3=-0.0114728378711725;
  // phase=0.474502658857982;
  // // Segment 10: # all p1
  // s1=0.044877061610554;
  // s2=0.0185609989287175;
  // s3=-0.00487110286831302;
  // c1=-0.00873239374125142;
  // c2=0.00535212227596913;
  // c3=0.0118521680070813;
  // phase=0.265147493586845;
  // // Segment 11: all p1
  // s1=0.0333820121123095;
  // s2=0.0298712223272543;
  // s3=-0.0142904637785188;
  // c1=-0.0285340288230717;
  // c2=0.0141638755571451;
  // c3=0.00347773849219099;
  // phase=0.812618911466807;
  // // Segment 12:
  // s1=0.058848328947128;
  // s2=0.0579931636636079;
  // s3=-0.00650741790427339;
  // c1=-0.0174804304809751;
  // c2=0.0319387676817976;
  // c3=-0.0134064461450929;
  // phase=0.810430465356286;
  // // Segment 13:
  // s1=0.0310281412891829;
  // s2=0.024100353802362;
  // s3=-0.0409784618318223;
  // c1=-0.0457979351738254;
  // c2=0.0372572537667035;
  // c3=-0.0246827764449057;
  // phase=0.7558637402789;
  // // Segment 14: all p1
  // s1=0.016825987386511;
  // s2=0.036149796983199;
  // s3=-0.0243642899376438;
  // c1=-0.0235042230762647;
  // c2=0.0359387144748583;
  // c3=-0.0336943553354116;
  // phase=0.537464751264087;
  // // Segment 15: all p1
  // s1=0.0274731308365419;
  // s2=-0.0211948021889934;
  // s3=0.0246271055418431;
  // c1=-0.0268581041717614;
  // c2=0.0339615626354796;
  // c3=-0.0350802270116981;
  // phase=0.458292771646705;
  // // Segment 16: all p1
  // s1=0.0289497246211529;
  // s2=-0.0466254501989807;
  // s3=0.0462223348968031;
  // c1=-0.0181422173373687;
  // c2=0.0364219181619443;
  // c3=-0.016092807746204;
  // phase=0.427847004573772;
  // // Segment 17: all p1
  // s1=-0.0111255551060104;
  // s2=-0.0409529069051457;
  // s3=0.050991801925046;
  // c1=-0.00831759658514021;
  // c2=0.0248545311264066;
  // c3=-0.0206349371716427;
  // phase=0.620646372652932;



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
