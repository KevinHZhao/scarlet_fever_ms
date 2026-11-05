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
  //  b0=3.39287184232719;
  // s1=-0.0655194901406117;
  // s2=0.00823990602912866;
  // s3=-0.0183904888399301;
  // c1=-0.0606363544810858;
  // c2=-0.00991779849784657;
  // c3=-0.00321795004458909;
  // phase=0.60958722352732;
  // Segment 2:
  //   b0=3.19659775932005;
  // s1=-0.0866782879607374;
  // s2=-0.00157029015261022;
  // s3=-0.0470438366633819;
  // c1=-0.020776251987272;
  // c2=0.00902727486323518;
  // c3=-0.00625840788558551;
  // phase=0.592657420884138;
  // Segment 3:
  //   b0=2.85717798360749;
  // s1=-0.0885971606801825;
  // s2=-0.0411492984721109;
  // s3=-0.0451675674790941;
  // c1=-0.000901197621920214;
  // c2=0.013480405769896;
  // c3=-0.00668767691371222;
  // phase=0.118158768379424;
  // Segment 4:
  //   b0=2.4939022207431;
  // s1=-0.0787169714464135;
  // s2=-0.0585880898217885;
  // s3=-0.0135916251612717;
  // c1=-0.00475501508015189;
  // c2=0.0426436922231834;
  // c3=-0.0235735013469128;
  // phase=0.177906455506854;
  // Segment 5:
  //   b0=2.28383361033798;
  // s1=-0.0948784123073904;
  // s2=-0.0371229430526528;
  // s3=0.00406169196460477;
  // c1=0.0258221803482324;
  // c2=0.0112255859397309;
  // c3=-0.00387197124126681;
  // phase=0.218969267061612;
  // Segment 6:
  //   b0=2.10135446402971;
  // s1=-0.100572465900293;
  // s2=-0.017282749766981;
  // s3=-0.0173239990537396;
  // c1=0.00060493151162415;
  // c2=0.0230871148737975;
  // c3=-0.00834560064935748;
  // phase=0.210790984023636;
  // Segment 7:
  //   b0=1.73745015637891;
  // s1=-0.103567542959195;
  // s2=-0.030939458175431;
  // s3=-0.0176403221349948;
  // c1=0.012788021858463;
  // c2=0.0288400217593975;
  // c3=0.00460265799136907;
  // phase=0.183774771170995;
  // Segment 8:
  //   b0=1.41900959013377;
  // s1=-0.0970199669468875;
  // s2=-0.0169040062933452;
  // s3=-0.0203265560658011;
  // c1=-0.00918483746593198;
  // c2=0.0359361284382128;
  // c3=0.00346549767517426;
  // phase=0.590168551553583;
  // Segment 9:
  //   b0=1.22224001531413;
  // s1=-0.0893149413887708;
  // s2=-0.030310477925732;
  // s3=-0.0111517638987014;
  // c1=0.0141991197196617;
  // c2=0.003703580827222;
  // c3=-0.00567592591391861;
  // phase=0.166191107754234;
  // Segment 10:
  //   b0=1.0773169769846;
  // s1=-0.0426762121837302;
  // s2=-4.45548378600785e-05;
  // s3=0.0107146968009237;
  // c1=-0.0174594511243675;
  // c2=0.0189404453451323;
  // c3=0.010050750201946;
  // phase=0.675400328625866;
  // Segment 11:
  //   b0=0.992605639988173;
  // s1=-0.0445951321898915;
  // s2=-0.00789993840260852;
  // s3=-0.00663938597252768;
  // c1=0.00511683426751281;
  // c2=0.0389264724103934;
  // c3=0.0115580043916191;
  // phase=0.583989287934496;
  // Segment 12:
  //   b0=0.938205221629681;
  // s1=-0.0431993485536462;
  // s2=-0.0187039597021814;
  // s3=-0.0242023484520815;
  // c1=-0.0233897267335901;
  // c2=0.0480107841337464;
  // c3=-0.00694931243185736;
  // phase=0.538058782133647;
  // Segment 13:
  //   b0=0.862845433811439;
  // s1=-0.0548927577150481;
  // s2=-0.0290501339658798;
  // s3=-0.0432144952297974;
  // c1=-0.00242225353132945;
  // c2=0.023479438014991;
  // c3=0.0136125471558432;
  // phase=0.590906752036159;
  // Segment 14:
  //   b0=0.802685937309344;
  // s1=-0.0468121553030855;
  // s2=-0.0383880427487022;
  // s3=-0.0499381138419779;
  // c1=0.0237573162578808;
  // c2=0.029257574411763;
  // c3=0.0345596706265015;
  // phase=0.609253394869702;
  // Segment 15:
  //   b0=0.786151685669694;
  s1=-0.0627830366423996;
  s2=-0.0387211585999606;
  s3=-0.0460442248511711;
  c1=0.035629360324041;
  c2=0.0423549532614044;
  c3=0.0158225540971261;
  phase=0.151847648786739;
  // Segment 16:
  //   b0=0.773845064962346;
  // s1=-0.0284881047987553;
  // s2=-0.0280723424180925;
  // s3=-0.0489187143265773;
  // c1=0.00605674033632728;
  // c2=0.0448938747476485;
  // c3=0.0158094395087542;
  // phase=0.746003261624047;
  // Segment 17:
  //   b0=0.756941432579632;
  // s1=-0.0248509376135688;
  // s2=-0.0361838983343244;
  // s3=-0.0424311774442514;
  // c1=0.0221345675485541;
  // c2=0.0450958494440633;
  // c3=0.0322332297849662;
  // phase=0.750639663561295;
  // Segment 18:
  //   b0=0.751952724120715;
  // s1=-0.00753341180540277;
  // s2=-0.0410036224065759;
  // s3=-0.0407385306022414;
  // c1=0.019937262970521;
  // c2=0.0273353565112439;
  // c3=0.0252007309677575;
  // phase=0.311284519518978;


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
