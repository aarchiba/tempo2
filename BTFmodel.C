#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

/* Based on BTXmodel.C */

/* Should support:
 * PB
 * PBDOT
 * FBAn - n>0, Fourier cosine coefficient
 * FBBn - n>0, Fourier sine coefficient
 * BTFSPAN - the length of one period of the fundamental, in days
 *
 * Orbital period is:
 * 1/(PB + PBDOT*(T-T0))/86400 + sum_n (PBAn * cos(w*(T-T0))+PBBn * sin(w*(T-T0)))
 * where w=2*pi/BTFSPAN
 * FIXME: better to make orbital frequency vary sinusoidally
 * */

double BTFmodel_i(pulsar *psr,int p,int ipos,int param,int arr)
{
  double torb;
  double tt0;
  double orbits;
  double pb;     /* Orbital period (sec) */
  double pbdot;
  double xpbdot;
  double ecc;    /* Orbital eccentricity */
  double edot;
  double asini;
  double xdot;
  double omdot;
  double omega;
  double gamma;
  int    norbits;
  double phase;
  double ep,dep,bige,tt,som,com;
  double alpha,beta,sbe,cbe,q,r,s;

  int i;
  double wspan;

  /* FIXME: PB must be set or else tempo2 doesn't believe it's a binary */
  tt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_t0].val[0])*SECDAY;

  pb     = psr[p].param[param_pb].val[0] * SECDAY;
  pbdot  = psr[p].param[param_pbdot].val[0];
  edot   = 0.0;
  ecc    = psr[p].param[param_ecc].val[0] + edot*tt0;

  if (psr[p].param[param_a1].paramSet[0]) {
      asini = psr[p].param[param_a1].val[0];
  } else {
      asini = 0;
  }
  if (psr[p].param[param_a1dot].paramSet[0]) {
      asini += psr[p].param[param_a1dot].val[0]*tt0;
  }

  if (ecc < 0.0 || ecc > 1.0)
    {
      printf("BTFmodel: problem with eccentricity = %Lg\n",psr[p].param[param_ecc].val[0]);
      exit(1);
    }

  if (psr[p].param[param_omdot].paramSet[0] == 1) omdot = psr[p].param[param_omdot].val[0];
  else omdot  = 0.0;
  omega  = (psr[p].param[param_om].val[0] + omdot*tt0/(SECDAY * 365.25))/(180.0/M_PI);
  if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
  else gamma  = 0.0;

  torb = 0.0;
  /* FIXME: use longdoubles in these calculations? */
  orbits = tt0/pb - pbdot/(pb*pb)*tt0*tt0/2.;
  if (psr[p].param[param_btfspan].paramSet[0]) {
      wspan = 2*M_PI/(SECDAY*psr[p].param[param_btfspan].val[0]);
      for (i=1;i<=MAX_BTF_TERMS;i++) {
          double a, b;
          a = psr[p].param[param_fban].val[i-1]*cos(wspan*i*tt0)/(wspan*i);
          b = psr[p].param[param_fbbn].val[i-1]*sin(wspan*i*tt0)/(wspan*i);
          if(0 && psr[p].param[param_fban].paramSet[i-1])
              printf("FBA%d adjustment: %g\tFBB%d adjustment: %g\n",
                      i,a,i,b);
          orbits += a+b;
      }
  } else {
      printf("Error: BTFSPAN not set; model cannot work\n");
      exit(2);
  }
  norbits = (int)orbits;
  if (orbits < 0.0) norbits--;
  
  phase = 2.0*M_PI * (orbits-norbits);

  /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
  ep = phase + ecc*sin(phase)*(1.0+ecc*cos(phase));

  /* This line is wrong in the original tempo: should be inside the do loop */
  /*  denom = 1.0 - ecc*cos(ep);*/
  
  do {
    dep = (phase - (ep-ecc*sin(ep)))/(1.0 - ecc*cos(ep));
    ep += dep;
  } while (fabs(dep) > 1.0e-12);
  bige = ep;

  tt = 1.0-ecc*ecc;
  som = sin(omega);
  com = cos(omega);

  alpha = asini*som;
  beta = asini*com*sqrt(tt);
  sbe = sin(bige);
  cbe = cos(bige);
  q = alpha * (cbe-ecc) + (beta+gamma)*sbe;
  r = -alpha*sbe + beta*cbe;
  s = 1.0/(1.0-ecc*cbe);

  torb = -q+(2*M_PI/pb)*q*r*s + torb;

  //printf("tt0:\t%g\tphase:\t%d\t%g\ttorb:\t%g\n", tt0, norbits, (orbits-norbits), torb);
  /* torb is the time correction to move the pulses to the binary barycenter */
  if (param==-1) return torb;

  /* FIXME: make sure I know what this does */
  /* FIXME: I think it needs to be able to return the partial derivative
   * of torb with respect to each parameter
   * Unfortunately if these are wrong, convergence is slow and the
   * uncertainties are wrong, but failure is not otherwise obvious */
  if (param==param_pb)
    return -2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb) * SECDAY;  /* fctn(12+j) */
  else if (param==param_a1)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt));                /* fctn(9+j) */
  else if (param==param_ecc)
    return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt; /* fctn(10+j) */
  else if (param==param_om)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe);          /* fctn(13+j) */
  else if (param==param_t0)
    return -2.0*M_PI/pb*r*s*SECDAY;                           /* fctn(11+j) */
  else if (param==param_pbdot)
    return 0.5*(-2.0*M_PI*r*s/pb*SECDAY*tt0/(SECDAY*pb))*tt0; /* fctn(18+j) */
  else if (param==param_a1dot)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt))*tt0;            /* fctn(24+j) */
  else if (param==param_omdot)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe)*tt0;      /* fctn(14+j) */
  else if (param==param_edot)                            
    return (-(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt)*tt0; /* fctn(25+j) */
  else if (param==param_gamma) 
    return sbe;                                               /* fctn(15+j) */
  else if (param==param_fban) {
    /* The formula for pb_sec is -2.0*M_PI*r*s/pb_sec*tt0/pb_sec
     * dphase/dpb_sec is tt0/pb_sec
     * so this formula is -2.0*M_PI*r*s/pb_sec*dphase/dpb_sec
     * The formula for pbdot_sec is 0.5*(-2.0*M_PI*r*s/pb_sec*tt0/pb_sec)*tt0
     * The formula for t0_sec is -2.0*M_PI/pb_sec*r*s
     * If you have some parameter P that affects ony the (zero to two pi) 
     * orbital phase
     * then its partial should be -r*s*fb[0]*dphase/dP 
     * Or should it? should that really be fb[0] or the binary frequency
     * at the moment of the observation? I think the latter, but the rest
     * of the code ignores this, assuming that pb never changes enough to
     * matter.
     * Check the factors of 2 and M_PI though.
     * phase is 2*M_PI*sum(fb[i]*pow(tt0,i+1)/factorial(i+1))
     * dphase/dfb[i] = 2*M_PI*pow(tt0,i+1)/factorial(i+1)
     * */
      /* FIXME: how do I test this? */
      return -2.0*M_PI*r*s*cos(wspan*(arr+1)*tt0)/(wspan*(arr+1));  
  } else if (param==param_fbbn) {
      return -2.0*M_PI*r*s*sin(wspan*(arr+1)*tt0)/(wspan*(arr+1));  
  }
  return 0.0;
}
double BTFmodel(pulsar *psr,int p,int ipos,int param,int arr)
{
    double h, v, l, r, d;
    int report_value = 0;
    if (param==-1) {
        int i,j;
        if (report_value)
            printf("Value call\n");
        if (report_value)
            for (i=0;i<MAX_PARAMS;i++)
                for (j=0;j<psr[p].param[i].aSize;j++)
                    if (psr[p].param[i].paramSet[j]) 
                        printf("\t%s:\t%Lg\n", 
                                psr[p].param[i].label[j], 
                                psr[p].param[i].val[j]);


        v = BTFmodel_i(psr,p,ipos,param,arr);
        if (report_value)
            printf("returned %g\n", v);
        return v;
    }

    d = BTFmodel_i(psr,p,ipos,param,arr);
    printf("Deriv. of obs. %d w.r.t. %s:\t%g",
            ipos, psr[p].param[param].label[arr], d);
    v = psr[p].param[param].val[arr];
    h = psr[p].param[param].err[arr];
    if (h==0) {
        h=v*1e-8;
        if (v==0) 
            h=1e-16;
    }
    psr[p].param[param].val[arr] = v-h;
    l = BTFmodel_i(psr,p,ipos,-1,arr);
    psr[p].param[param].val[arr] = v+h;
    r = BTFmodel_i(psr,p,ipos,-1,arr);
    psr[p].param[param].val[arr] = v;

    printf("\tnumerical:\t%g\n", (r-l)/(2*h));
    return d;
}




void updateBTF(pulsar *psr,double val,double err,int pos,int arr)
{
  //printf("Updating parameter %s by %g (err %g)\n", psr->param[pos].label[arr], val, err);
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val;
      psr->param[param_pb].err[0]  = err;
    }
  else if (pos==param_pbdot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_fban || pos==param_fbbn)
    {
      /* I think it's correct to not update the PB/PBDOT */
      psr->param[pos].val[arr] += val;
      psr->param[pos].err[arr]  = err;
    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_om)
    {
      psr->param[pos].val[0] += val*180.0/M_PI;
      psr->param[pos].err[0]  = err*180.0/M_PI;
    }
  else if (pos==param_omdot)
    {
      psr->param[pos].val[0] += val*(SECDAY*365.25)*180.0/M_PI;
      psr->param[pos].err[0]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
  else if (pos==param_a1dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
}
