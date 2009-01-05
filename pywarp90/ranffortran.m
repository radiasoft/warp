/* Created by David P. Grote, January 1, 2004 */
/* $Id: ranffortran.m,v 1.3 2009/01/05 19:19:51 dave Exp $ */

/* This is needed since the ranf in ranlib is single precision */
extern double Ranf(void);
%'double '+fname(py_ifelse(f90 or f90f,1,'w','')+'ranf')+'(void)'
{
  return Ranf();
}

%py_ifelse(machine,'T3E','extern void RANSET(double *x);','extern void Seedranf(double *x);')
%'void '+fname('seedranf')+'(double* x)'
{
%py_ifelse(machine,'T3E','RANSET(x);','Seedranf(x);')
}

extern void Mixranf();
%'void '+fname('mixranf')+'(int* x)'
{
  double seed[2];
  Mixranf(x,&seed);
}

