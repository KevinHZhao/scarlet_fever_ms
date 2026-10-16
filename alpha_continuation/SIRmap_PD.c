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

  // Segment 1:
  // b0=3.39278853365536;
  // s1=-0.0655276082824137;
  // s2=0.00824144803439916;
  // s3=-0.018399456165716;
  // c1=-0.0606437664618818;
  // c2=-0.00991868273890906;
  // c3=-0.00322056385920645;
  // phase=0.609576972544188;
  // Segment 2:
  // b0=3.19664899494418;
  // s1=-0.0866906527989088;
  // s2=-0.00156862833537588;
  // s3=-0.0470677563090881;
  // c1=-0.0207787995885545;
  // c2=0.00903252019960112;
  // c3=-0.00626153349203608;
  // phase=0.592650411453936;
  // Segment 3:
  // b0=2.83844268693297;
  // s1=-0.0894017094607646;
  // s2=-0.0457661231138116;
  // s3=-0.0415316794380626;
  // c1=-0.00246710921226798;
  // c2=0.0163327623460009;
  // c3=-0.00795714732048603;
  // phase=0.121322810209603;
  // Segment 4:
  // b0=2.45227400753782;
  // s1=-0.07628424370166;
  // s2=-0.0569062667716558;
  // s3=-0.0119888032601803;
  // c1=-0.0038296172955699;
  // c2=0.044668811259658;
  // c3=-0.0251310527834427;
  // phase=0; // trying this to avoid crazy behaviour
  // Segment 5:
  // b0=2.27547442507189;
  // s1=-0.0917938998950198;
  // s2=-0.033522570948086;
  // s3=-0.000547733253459533;
  // c1=0.0195790736372547;
  // c2=0.0117719073257379;
  // c3=-0.00275474026153355;
  // phase=0.209320899361362;
  // Segment 6:
  // b0=2.04044007815409;
  // s1=-0.10787818701876;
  // s2=-0.0151421311211967;
  // s3=-0.0185170266360679;
  // c1=0.000519891647789314;
  // c2=0.0269727354208919;
  // c3=-0.0119281275676661;
  // phase=0; //trying this to avoid crazy behaviour
  // Segment 7:
  // b0=1.73573443494335;
  s1=-0.103582173677835;
  s2=-0.0309469743380368;
  s3=-0.017645829897115;
  c1=0.0127940377074972;
  c2=0.0288483637220368;
  c3=0.00460465475563706;
  phase=0.183772284157415;
  // Segment 8:
  // b0=1.41625798524637;
  // s1=-0.0970418864744058;
  // s2=-0.0169068509072327;
  // s3=-0.020335693964556;
  // c1=-0.0091917664220574;
  // c2=0.035944823441271;
  // c3=0.00346362353372541;
  // phase=0.59016043802734;
  // Segment 9:
  // b0=1.21954403753408;
  // s1=-0.0893325997905696;
  // s2=-0.0303197086463664;
  // s3=-0.0111572101620274;
  // c1=0.0142004521007489;
  // c2=0.0037010591646133;
  // c3=-0.00567931672665393;
  // phase=0.166167168724369;
  // Segment 10:
  // b0=1.07448149973967;
  // s1=-0.042680804409166;
  // s2=-3.82976067914358e-05;
  // s3=0.0107175892128128;
  // c1=-0.0174641031838146;
  // c2=0.0189395299017891;
  // c3=0.0100516052369638;
  // phase=0.675400535628241;
  // Segment 11:
  // b0=0.989538872628786;
  // s1=-0.0446766241083989;
  // s2=-0.00793686410310266;
  // s3=-0.00658240406045131;
  // c1=0.00508729340508013;
  // c2=0.0390109892879331;
  // c3=0.0116111355083008;
  // phase=0.584030626916746;
  // Segment 12:
  // b0=0.935028079172278;
  // s1=-0.0433614216168957;
  // s2=-0.0187542663534447;
  // s3=-0.0241175532057328;
  // c1=-0.0234431481622157;
  // c2=0.0481565639209478;
  // c3=-0.00689565013479413;
  // phase=0.537985778211923;
  // Segment 13:
  // b0=0.859293352075038;
  // s1=-0.0548943261460712;
  // s2=-0.0290493508060133;
  // s3=-0.0432148845071454;
  // c1=-0.00242434445468044;
  // c2=0.0234776323500825;
  // c3=0.0136077298347442;
  // phase=0.590902252146415;
  // Segment 14:
  // b0=0.797188237132598;
  // s1=-0.0570726550884129;
  // s2=-0.0391729574453816;
  // s3=-0.0494420549265408;
  // c1=0.0309477808675647;
  // c2=0.0335969610465321;
  // c3=0.0303667947840144;
  // phase=0.14627398508778;
  // Segment 15:
  // b0=0.767711875734533;
  // s1=-0.0342242164178104;
  // s2=-0.0307404273134962;
  // s3=-0.0466738869851473;
  // c1=0.0122540730267082;
  // c2=0.0448035588163104;
  // c3=0.0120254127679637;
  // phase=0.931005193206847;
  // Segment 16:
  // b0=0.752454375213868;
  // s1=-0.0248560157748592;
  // s2=-0.0361895938677943;
  // s3=-0.0424388655465249;
  // c1=0.0221399636151382;
  // c2=0.0451023885489337;
  // c3=0.0322384813565753;
  // phase=0.75063924459045;
  // Segment 17:
  // b0=0.747235396981138;
  // s1=-0.00753360678023388;
  // s2=-0.0410122716976199;
  // s3=-0.0407462563020793;
  // c1=0.019938275982936;
  // c2=0.0273392345897244;
  // c3=0.0252026144964448;
  // phase=0.311283678569734;


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
