#ifndef __NUMERICALRECIPES_HPP
#define __NUMERICALRECIPES_HPP

#include <vector>
#include <stdio.h>
#include <math.h>

#include "Function.hpp"
#include "Matrix.hpp"

using namespace std;
using std::vector;

class Function;
 
class NumericalRecipes {

public:
  NumericalRecipes() {;}
  ~NumericalRecipes() {;}

  valarray<double>*fvec_p; 

  double (NumericalRecipes::*aFunc_)(valarray<double> &x, Function &vecfunc);

  //! Compute forward-difference approximation to Jacobian (from Numerical Recipes).
  void fdjac(valarray<double> &x, valarray<double> &fvec,
	     Matrix<double> &df, Function &vecfunc);

  //! Returns f=1/2(F*F) at point x (from Numerical Recipes, modified).
  double fmin(valarray<double> &x, Function &vecfunc);

  //! Given an n-dimensional point (xold), the value of the function (fold) and 
  //! gradient there (g), and a direction (p), finds a new point (x) along the 
  //! direction p from xold where the function func has decreased "sufficiently" 
  //! (from Numerical Recipes).      
  void lnsrch(valarray<double> &xold, const double fold,
	      valarray<double> &g, valarray<double> &p,
	      valarray<double> &x, double &f, const double stpmax,
	      bool &check, Function & vecfunc);

  //! Solve the set of n linear equations A*X = B (from Numerical Recipes).   
  void lubksb(Matrix<double> &a, valarray<int> &indx,
	      valarray<double> &b);

  //! Given a matrix, this routine replaces it by the LU decomposition of a rowwise 
  //! permutation of itself (from Numerical Recipes).    
  void ludcmp(Matrix<double> &a, valarray<int> &indx, double &d);

  //! Given an initial guess (x) for a root in n dimensions, find the root by a 
  //! globally convergent Newton's method.    
  void newt(valarray<double> &x, bool &check,
	    Function &vecfunc);
  //!function for the secant method (from Numerical Recipes(rtsec), modified).  
  int secantMethod( Function &func, double &x1,double &x2,
		    const double &precision,double &root); 

};

#endif
