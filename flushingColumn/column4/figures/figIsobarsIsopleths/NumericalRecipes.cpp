#include "NumericalRecipes.hpp"

void NumericalRecipes::fdjac(valarray<double> &x, valarray<double> &fvec,
			     Matrix<double> &df,
			     Function &vecfunc)
{
  const double EPS=1.0e-8; //approximate square root of the machine precision
  int i,j;
  double h,temp;
  int n=x.size(); 
  valarray<double> f(n);
  for (j=0;j<n;j++) 
    { 
      temp=x[j];
      h=EPS*fabs(temp);
      if (h == 0.0) h=EPS;
      x[j]=temp+h; //trick to reduce finite precision error
      h=x[j]-temp;
      vecfunc.function2(x,f);
      // 
      x[j]=temp;
      for (i=0;i<n;i++) 
	{
	  df[i][j]=(f[i]-fvec[i])/h; //forward difference formula
	}
    }
}

double NumericalRecipes::fmin(valarray<double> &x, Function &vecfunc)
{
  int i; 
  double sum;
  valarray<double> &fvec=*fvec_p;
  vecfunc.function2(x,fvec);
  int n=x.size();
  sum=0.0;
  for (i=0;i<n;i++) 
    {
      sum += (fvec[i]*fvec[i]);
    }
  return 0.5*sum;
} 


void NumericalRecipes::lnsrch(valarray<double> &xold, const double fold,
			      valarray<double> &g, valarray<double> &p,
			      valarray<double> &x, double &f, const double stpmax,
			      bool &check, Function & vecfunc)
{
  const double ALF=1.0e-6, TOLX = 1.0e-12;
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double fold2 = 0;
  double rhs1,rhs2,slope,sum,temp,test,tmplam;
  int n=xold.size();
  check=false; 
  sum=0.0; 
  for (i=0;i<n;i++) 
    {
      sum += p[i]*p[i]; 
    }
  sum=sqrt(sum);
  if (sum > stpmax)
    {
      for (i=0;i<n;i++) 
	{
	  p[i] *= stpmax/sum; //scale if attempted step is too big
	}
    }
  slope=0.0;
  for (i=0;i<n;i++) 
    {
      slope += g[i]*p[i];
    }
  if (slope >= 0.0)
    {
      cerr <<"Roundoff problem in lnsrch." << endl;
      exit(1);
    }
  //
  test=0.0; // for compute lambda_min(alamin)
  for (i=0;i<n;i++) 
    {
      temp=fabs(p[i])/max(fabs(xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;//minimum step lenght
  alam=1.0; //always try full Newton step first
  for (;;) 
    {
      for (i=0;i<n;i++) 
	{
	  x[i]=xold[i]+alam*p[i];
	}
      f=(this->*aFunc_)(x, vecfunc);
	
      if (alam < alamin) //convergence on Delta_x
	{
	  for (i=0;i<n;i++) 
	    {
	      x[i]=xold[i];
	    }
	  check=true;
	  return;
	} 
      else if (f <= fold+ALF*alam*slope) //sufficient function decrease
	{
	  return;
	}
      else //backtrack
	{
	  if (alam == 1.0) //first time
	    {
	      tmplam = -slope/(2.0*(f-fold-slope));
	    }
	  else //subsequent backtracks
	    {
	      rhs1=f-fold-alam*slope;
	      rhs2=f2-fold-alam2*slope;

	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/
		 (alam2*alam2))/(alam-alam2);
	      if (a == 0.0) 
		{
		  tmplam = -slope/(2.0*b);
		}
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc < 0.0)  
		    { 
		      tmplam=0.5*alam; 
		    } 
		  else if (b <= 0.0)  
		    { 
		      tmplam=(-b+sqrt(disc))/(3.0*a); 
		    } 
		  else  
		    { 
		      tmplam=-slope/(b+sqrt(disc)); 
		    } 
		}
	      if (tmplam>0.5*alam) 
		{
		  tmplam=0.5*alam; //lambda <= 0.5*lambda1 (alam)
		}
	    }
	}
      alam2=alam;
      f2 = f;
      alam=max(tmplam,0.1*alam); //lambda >= 0.1*lambda1
    }
}
  

void NumericalRecipes::lubksb(Matrix<double> &a, valarray<int> &indx,
			      valarray<double> &b)
{
  int i,ii=0,ip,j;
  double sum;
  int n=a.nRows();
  for (i=0;i<n;i++) 
    {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii != 0)
	{
	  for (j=ii-1;j<i;j++) 
	    {
	      sum -= a[i][j]*b[j];
	    }
	}
      else if (sum != 0.0)
	{
	  ii=i+1;
	}
      b[i]=sum;
    }
  for (i=n-1;i>=0;i--) 
    {
      sum=b[i];
      for (j=i+1;j<n;j++) 
	{
	  sum -= a[i][j]*b[j];
	}
      b[i]=sum/a[i][i];
    }
}
  
void NumericalRecipes::ludcmp(Matrix<double> &a, valarray<int> &indx, double &d)
{
  const double TINY=1.0e-20;
  int i,imax,j,k;
  double big,dum,sum,temp;
  int n=a.nRows();
  valarray<double> vv(n); //vv stores the implicit scaling of each row
  d=1.0;
  for (i=0;i<n;i++)  //loop over rows to get the implicit scaling information
    {
      big=0.0;
      for (j=0;j<n;j++)
	{
	  if ((temp=fabs(a[i][j])) > big) 
	    {
	      big=temp;
	    }
	}
      if (big == 0.0)
	{
	  cerr << "Singular matrix in routine ludcmp" << endl;
	  exit(1);
	}
      vv[i]=1.0/big;
    }
  for (j=0;j<n;j++) //loop over columns of crout's method
    {
      for (i=0;i<j;i++) 
	{
	  sum=a[i][j];
	  for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	  a[i][j]=sum;
	}
      big=0.0;
      for (i=j;i<n;i++) 
	{
	  sum=a[i][j];
	  for (k=0;k<j;k++) 
	    {
	      sum -= a[i][k]*a[k][j];
	    }
	  a[i][j]=sum;
	  if ((dum=vv[i]*fabs(sum)) >= big) 
	    {
	      big=dum;
	      imax=i;
	    }
	}
      if (j != imax)
	{
	  for (k=0;k<n;k++) 
	    {
	      dum=a[imax][k];
	      a[imax][k]=a[j][k];
	      a[j][k]=dum;
	    }
	  d = -d;
	  vv[imax]=vv[j];
	}
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n-1) 
	{
	  dum=1.0/(a[j][j]);
	  for (i=j+1;i<n;i++) a[i][j] *= dum;
	}
    }
}


void NumericalRecipes::newt(valarray<double> &x, bool &check,
			    Function &vecfunc)
{
  const int MAXITS=10000;
  const double TOLF=1.0e-12,TOLMIN=1.0e-12,STPMX=100.0;
  const double TOLX=1.0e-13, tolstp=1.e-15;
  //TOLF= set the convergence criterion on function values
  //TOLMIN=sets the criterion for deciding whether spurious convergence 
  //       to a minimum of fmin has occurred
  //TOLX=convergence criterion on delta_x

  int i,j,its;
  double d,den,f,fold,stpmax,sum,temp,test;
  int n=x.size();
  valarray<int> indx(n);
  valarray<double> g(n),p(n),xold(n);
  Matrix<double> fjac(n,n);
  fvec_p=new valarray<double>(n);
  valarray<double> &fvec=*fvec_p;
  f=fmin(x, vecfunc);
  test=0.0; //test for initial guess being a root. use more stringent test than simply TOLF
  for (i=0;i<n;i++)
    {
      if (fabs(fvec[i]) > test) 
	{
	  test=fabs(fvec[i]);
	}
    }
  if (test < 0.01*TOLF) 
    {
      check=false;
      delete fvec_p;
      return;
    } 
  sum=0.0;
  for (i=0;i<n;i++) 
    {
      sum += (x[i]*x[i]);
    }
  stpmax=STPMX*max(sqrt(sum),double(n)); //calculate stpmax for line searches
  stpmax = x[0]*(1.-tolstp);      //stpmax definition 
  double xless = x[0];  
  for(int istp=1; istp<n;istp++)  
    {  
      if(x[istp] < xless)  
	{  
	  xless = x[istp];  
	  stpmax = x[istp]*(1.-tolstp);  
	} 
    }  
  for (its=0;its<MAXITS;its++) //start of iteration loop
    {
      fdjac(x,fvec,fjac,vecfunc);
      for (i=0;i<n;i++) //compute Nabla_f for the line search
	{
	  sum=0.0;
	  for (j=0;j<n;j++) 
	    {
	      sum += fjac[j][i]*fvec[j];
	    }
	  g[i]=sum;
	}
      for (i=0;i<n;i++) 
	{
	  xold[i]=x[i]; //store x
	}
      fold=f; //store f
      for (i=0;i<n;i++) 
	{
	  p[i] = -fvec[i]; //right-hand side for linear equations
	}
      ludcmp(fjac,indx,d); //solve linear equation by LU decomposition
      lubksb(fjac,indx,p); 
      aFunc_ = & NumericalRecipes::fmin;
      lnsrch(xold,fold,g,p,x,f,stpmax,check,vecfunc);
      //returns new x and f.it also calculates fvec at the new x when it calls fmin
      test=0.0;
      for (i=0;i<n;i++)
	{
	  if (fabs(fvec[i]) > test) 
	    {
	      test=fabs(fvec[i]);
	    }
	}
      if (test < TOLF) 
	{
	  check=false;
	  delete fvec_p;
	  return;
	}
      if (check) 
	{
	  test=0.0;
	  den=max(f,0.5*n);
	  for (i=0;i<n;i++) 
	    {
	      temp=fabs(g[i])*max(fabs(x[i]),1.0)/den;
	      if (temp > test) 
		{
		  test=temp;
		}
	    }
	  check=(test < TOLMIN);
	  delete fvec_p;
	  return;
	}
      test=0.0; //test for convergence on delta_x
      for (i=0;i<n;i++) 
	{
	  temp=(fabs(x[i]-xold[i]))/max(fabs(x[i]),1.0);
	  if (temp > test) 
	    {
	      test=temp;
	    }
	}
      if (test < TOLX) 
	{
	  delete fvec_p;
	  return;
	}
    }
  cerr << "MAXITS exceeded in newt" << endl;
  exit(1);
  return;
}
  

int NumericalRecipes::secantMethod( Function &func, double &x1,double &x2,
				    const double &precision,double &root) 
{ 
  const int max_iter = 1000; //che valore devo mettere?
  double yl, yr, dx, swap, xl, rts;
    
  yl = func.function(x1);
  yr = func.function(x2);
  if(x1==x2 || yl==yr) 
    {
      cerr << "error in secantMethod call; x1==x2 or f(x1)==f(x2)" << endl;
      exit(1);
    }
  if(fabs(yl) < fabs(yr)) 
    {
      rts = x1;
      xl = x2;
      swap = yl;
      yl = yr;
      yr = swap;
    }
  else 
    {
      xl = x1;
      rts = x2;
    }
  for(int j = 1; j <= max_iter; j++) 
    {
      dx = (xl-rts)*yr/(yr-yl);
      xl = rts;
      yl = yr;
      rts += dx;
      yr = func.function(rts);
	
      if(fabs(dx) < precision || yr == 0.0) 
	{
	  root = rts;
	  return 1;
	}
    }
  return 0; 
}
  

