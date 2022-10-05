/*==============================================================================*
 * MATHS
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef MATHS_H
#define MATHS_H

#include <random>

double det4(double (*M)[4]);
double det3(double (*M)[3]);
void matVecMultiply(double *y, double**A, double *x, int nu);
void vecVecSubtract(double *, double *v1, double *v2, int nu);

const double PI = acos(-1.0);
double rnd();

#endif /* !MATHS_H */
