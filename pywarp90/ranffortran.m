/* Created by David P. Grote, January 1, 2004 */
/* $Id: ranffortran.m,v 1.3 2009/01/05 19:19:51 dave Exp $ */

/* This is needed since the ranf in ranlib is single precision */
extern double Ranf(void);
%'double '+fname('wranf')+'(void)'
{
  return Ranf();
}

extern void Seedranf(double *x);
%'void '+fname('seedranf')+'(double* x)'
{
Seedranf(x);
}

/* Mixranf is never used
extern void Mixranf();
%'void '+fname('mixranf')+'(int* x)'
{
  double seed[2];
  Mixranf(x,&seed);
}
*/

