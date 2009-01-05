/* Created by David P. Grote, January 1, 2004 */
/* $Id: ranffortran.m,v 1.2 2009/01/05 19:10:17 dave Exp $ */

/* This is needed since the ranf in ranlib is single precision */
extern double Ranf();
%'double '+fname(py_ifelse(f90 or f90f,1,'w','')+'ranf')+'()'
{
  return Ranf();
}

%py_ifelse(machine,'T3E','extern void RANSET();','extern void Seedranf();')
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

