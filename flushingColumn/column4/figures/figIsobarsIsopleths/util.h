#ifndef __UTIL_H
#define __UTIL_H

#include <complex>
#include <vector>
#include <stdio.h>
#include <math.h>

#include "Function.hpp"


using namespace std;
using std::complex;
using std::vector;


#ifndef UTIL_PI
#define UTIL_PI        3.14159265358979323846264338328
#endif

inline double Abs(double a)  { return (((a) > 0) ? (a) : -(a)); }

inline double Abs(complex<double> a)  { return (sqrt(a.real()*a.real()+a.imag()*a.imag())); }

inline double Conjg(double x) { return x; }

inline complex<double> Conjg(complex<double> x) { complex<double> y(x.real(),-x.imag()); return y; }

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

inline double Max(double a, double b, double c)  { return max(max(a,b),c); }
inline double Max(double a, double b, double c, double d)  { return max(max(max(a,b),c),d); }
inline int    Max(int    a, int    b, int    c)  { return max(max(a,b),c); }
inline double Min(double a, double b, double c)  { return min(min(a,b),c); }
inline int    Min(int    a, int    b, int    c)  { return min(min(a,b),c); }

inline double Arg(complex<double> x) {
  double a = 0.;
  if (x.real() == 0.)
    a = 0.5*UTIL_PI;
  else
    a = atan(x.imag()/x.real());
  return a;
}

template<class T_>
inline void Scale(T_ a, const vector<T_> &x, vector<T_> &y)
{
  typename std::vector<T_>::const_iterator j = x.begin();
  for (typename std::vector<T_>::iterator i=y.begin(); i!=y.end();)
    *(i++) = a * *(j++);
}

template<class T_>
inline void Scale(T_ a, vector<T_> &x)
{
  for (typename std::vector<T_>::iterator i=x.begin(); i!=x.end();)
    *(i++) *= a;
}

template<class T_>
inline void Xpy(int n, T_ *x, T_ *y)
{
  for (int i=0; i<n; i++)
    y[i] += x[i];
}

template<class T_>
inline void Xpy(const vector<T_> &x, vector<T_> &y)
{
  typename std::vector<T_>::const_iterator j = x.begin();
  for (typename std::vector<T_>::iterator i=y.begin(); i!=y.end();)
    *(i++) += *(j++);
}

template<class T_>
inline void Axpy(int &n, T_ a, T_ *x, T_ *y)
{
  for (int i=0; i<n; i++)
    y[i] += a*x[i];
}

template<class T_>
inline void Axpy(T_ a, const vector<T_> &x, vector<T_> &y)
{
  typename std::vector<T_>::const_iterator j = x.begin();
  for (typename std::vector<T_>::iterator i=y.begin(); i!=y.end();)
    *(i++) += a * *(j++);
}

template<class T_>
inline T_ Dot(int n, T_ *x, T_ *y)
{
  T_ d = T_(0);
   
  for (int i=0; i<n; i++)
    d += x[i]*y[i];
  return d;
}


template<class T_>
inline T_ Dot(const vector<T_> &x, const vector<T_> &y)
{
  typedef typename std::vector<T_>::const_iterator iter;
  T_ d = T_(0);
  iter j = x.begin();
  for (iter i=y.begin(); i!=y.end();)
    d += *(j++) * *(i++);
  return d;
}


inline complex<double> CDot(int n, const complex<double> *x, const complex<double> *y)
{
  complex<double> d=0;
  for (int i=0; i<n; i++)
    d += x[i]*Conjg(y[i]);
  return d;
}


inline complex<double> CDot(const vector<complex<double> > &x, vector<complex<double> > &y)
{
  complex<double> d=0;
  for (unsigned int i=0; i<x.size(); i++)
    d += x[i]*Conjg(y[i]);
  return d;
}

inline double Nrm2(int n, double *x) { return sqrt(Dot(n,x,x)/n); }
 
inline double Nrm2(const vector<double> &x) { return sqrt(Dot(x,x)); }



template<class T_>
inline void Clear(vector<T_> &v)
{
  for (typename std::vector<T_>::iterator i=v.begin(); i!=v.end();)
    *(i++) = T_(0);
}


inline char itoc(int i)
{
  static char buf[10];
  sprintf(buf,"%d",i);
  return buf[0];
}

inline void closeToOne(valarray<double> &x)
{
  double sum_x = 0;
  sum_x=x.sum();
  x/=sum_x;
}

inline double interpLinOnEquispacedVector(const double &x0, const valarray<double> &x, const valarray<double> &y)
{
  double delta_x = x[1]-x[0];
  int i = (int)floor((x0-x[0])/delta_x);
  double y0=y[i]+(y[i+1]-y[i])/delta_x*(x0-x[i]);
  return y0;
}


#endif
