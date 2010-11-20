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

/* Based on BTmodel.C */
/* Should agree with bnrybtx.f, but with only one orbit (i.e. no planets) */

/* Should support:
 * FBn - derivatives of orbital period
 * A1
 * XDOT
 * XDOTn - for n>1; derivatives of orbital size (?)
 *
 * Should support on input:
 * PB - binary period
 * PBDOT - binary period derivative
 *
 * */

double BTXmodel(pulsar *psr,int p,int ipos,int param,int arr)
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
  double fac;
  double fb[MAX_BTX_DERIVS];
  double x[MAX_BTX_DERIVS+2];

  /* FIXME: PB must be set or else tempo2 doesn't believe it's a binary */
  tt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_t0].val[0])*SECDAY;

  pb     = psr[p].param[param_pb].val[0] * SECDAY;
  edot   = 0.0;
  ecc    = psr[p].param[param_ecc].val[0] + edot*tt0;

  /* Handle input parameters */
  /* load fb and x */
  printf("Setting fb to");
  for (i=0; i<MAX_BTX_DERIVS; i++) {
      if (psr[p].param[param_fbn].paramSet[i]) {
          fb[i] = psr[p].param[param_fbn].val[i];
          printf("\t%g",fb[i]);
      } else {
          fb[i] = 0;
          printf("\t--",fb[i]);
      }
  }
  printf("\n");
  /* convert PB and PBDOT to FB0 and FB1 if necessary */
  if (psr[p].param[param_pb].paramSet[0] && !psr[p].param[param_fbn].paramSet[0]) {
      fb[0] = 1./pb;
      printf("Setting fb[0] to %g from pb\n", fb[0]);
  }
  if (psr[p].param[param_pbdot].paramSet[0] && !psr[p].param[param_fbn].paramSet[1]) {
      fb[1] = -(psr[p].param[param_pbdot].val[0]*SECDAY)*fb[0]*fb[0];
      printf("Setting fb[1] to %g from pbdot\n", fb[1]);
  }

  if (psr[p].param[param_a1].paramSet[0]) {
      x[0] = psr[p].param[param_a1].val[0];
  } else {
      x[0] = 0;
  }
  if (psr[p].param[param_a1dot].paramSet[0]) {
      x[1] = psr[p].param[param_a1dot].val[0];
  } else {
      x[1] = 0;
  }
  for (i=0; i<MAX_BTX_DERIVS; i++) {
      if (psr[p].param[param_xdotn].paramSet[i]) {
          x[i+2] = psr[p].param[param_xdotn].val[i];
      } else {
          x[i+2] = 0;
      }
  }

  if (ecc < 0.0 || ecc > 1.0)
    {
      printf("BTXmodel: problem with eccentricity = %Lg\n",psr[p].param[param_ecc].val[0]);
      exit(1);
    }

  asini = x[0] + x[1]*tt0;
  fac = 1.0;
  for (i=0; i<MAX_BTX_DERIVS; i++) {
      fac /= i+2;
      asini += x[i+2]*fac*pow(tt0,i+2);
  }


  if (psr[p].param[param_omdot].paramSet[0] == 1) omdot = psr[p].param[param_omdot].val[0];
  else omdot  = 0.0;
  omega  = (psr[p].param[param_om].val[0] + omdot*tt0/(SECDAY * 365.25))/(180.0/M_PI);
  if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
  else gamma  = 0.0;

  torb = 0.0;
  /* FIXME: use longdoubles in these calculations? */
  orbits = 0;
  fac = 1;
  for (i=0;i<MAX_BTX_DERIVS;i++) {
      fac /= i+1;
      orbits += fac*fb[i]*pow(tt0,i+1);
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

  printf("tt0:\t%g\tphase:\t%d\t%g\ttorb:\t%g\n", tt0, norbits, (orbits-norbits), torb);
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
  else if (param==param_fbn) {
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
      fac = 1.0;
      for (i=0;i<arr;i++) fac/=i+2;
      return -2.0*M_PI*r*s*fb[0]*pow(tt0,arr+1)*fac;  
  } else if (param==param_xdotn) {
      /* FIXME: work out this partial based on a1 and a1dot */
    return 0;                                               
  }
  return 0.0;
}

void updateBTX(pulsar *psr,double val,double err,int pos,int arr)
{
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val;
      psr->param[param_pb].err[0]  = err;
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
  else if (pos==param_pbdot)
    {
        /* FIXME: use the fbn instead */
        /* FIXME: also add other new quantities */
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
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
