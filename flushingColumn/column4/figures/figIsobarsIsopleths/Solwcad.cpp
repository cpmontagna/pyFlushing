#include "Solwcad.hpp"

//#define DEBUG
 
const double Solwcad::R_ = 8.31451;
const double Solwcad::p0_h2o_ = 0.;
const double Solwcad::p0_co2_ = 1.01325e5;
const double Solwcad::whc_ = 0.;
   
Solwcad::Solwcad() :
  wo_(12,12),
  wij_(20,2),
  molecular_weight_(12),
  ghiorso_par_(12),
  fol_(2),
  to_(2), 
  par_(10),
  gas_fraction_(2),
  y_(0.0,2),
  weight_fraction_(0.0,4), //xy(i)
  ctrm_(0.0,2),
  scale_factor_(1.0,2), 
  xtrm_(2),
  rtrm_(2),
  h2o_co2_(0.0,2),

  bm_ ( 0.),
  cm_ ( 0.),
  dm_ ( 0.),
  em_ ( 0.),
  vr_ ( 0.),
  xwo_(0.)
{ 
  molecular_weight_[0] =  18.02 ;
  molecular_weight_[1] =  44.01 ;
  molecular_weight_[2] =  60.085; 
  molecular_weight_[3] =  79.899;
  molecular_weight_[4] = 101.96 ; 
  molecular_weight_[5] = 159.69 ;
  molecular_weight_[6] =  71.846;
  molecular_weight_[7] =  70.937;
  molecular_weight_[8] =  40.304; 
  molecular_weight_[9] = 56.079;
  molecular_weight_[10] = 61.979;
  molecular_weight_[11] = 94.195;
#ifdef DEBUG
  cerr<<"//parametri di ghiorso and Sack 1983"<<endl;
#endif
  ghiorso_par_[0] = 1.;
  ghiorso_par_[1] = 1.;
  ghiorso_par_[2] = 0.25; 
  ghiorso_par_[3] = 0.25;
  ghiorso_par_[4] = 0.375; 
  ghiorso_par_[5] = 0.375;
  ghiorso_par_[6] = 0.25;
  ghiorso_par_[7] = 0.25;
  ghiorso_par_[8] = 0.25;
  ghiorso_par_[9] = 0.25;
  ghiorso_par_[10] = 0.375;
  ghiorso_par_[11] = 0.375;

  fol_[0] = 0.1979347855e+2;
  fol_[1] = 0.2348447513e+2;

  to_[0] = 0.1000000000000000e4;
  to_[1] = 0.1400000000000000e4;

  par_[0] = 0.2418946564e-05;
  par_[1] = 0.1344868775e-07;
  par_[2] = 0.0000000000e-00;
  par_[3] = 0.0000000000e-00;
  par_[4] = 0.2106552365e-14;
  par_[5] = 0.0000000000e-00;
  par_[6] = 0.0000000000e-00;
  par_[7] = 0.0000000000e-00;
  par_[8] = 0.0000000000e-00;
  par_[9] = 0.0000000000e-00;

  for(int i=0; i<12; i++)
    {
      wo_(i,i) = 0;
      for(int j = 0; j <2; j++) 
	{
	  wo_(j,i) = 0;
	}
    }
  wo_(2,3) = -122861.0680;
  wo_(2,4) = -328708.4288;
  wo_(2,5) = 11037.09912;
  wo_(2,6) = -40292.50576;
  wo_(2,7) = 23118.10624;
  wo_(2,8) = -126999.4624;
  wo_(2,9) = -268060.9304;
  wo_(2,10) = -308604.7272;
  wo_(2,11) = -366503.3376;

  wo_(3,4) = -281791.1448;
  wo_(3,5) = -28542.49488;
  wo_(3,6) = -19223.76456;
  wo_(3,7) = -8548.7488;
  wo_(3,8) = 53026.3424;
  wo_(3,9) = -428617.328;
  wo_(3,10) = -422893.616;
  wo_(3,11) = -170291.7288;

  wo_(4,5) = 5189.49888;
  wo_(4,6) = -249067.6624;
  wo_(4,7) = -8023.866;
  wo_(4,8) = -203655.3632;
  wo_(4,9) = -411824.0072;
  wo_(4,10) = -567413.16;
  wo_(4,11) = -733563.984;

  wo_(5,6) = 18930.34064;
  wo_(5,7) = 887.828064;
  wo_(5,8) = -5343.09352;
  wo_(5,9) = 6358.88504;
  wo_(5,10) = -15553.51792;
  wo_(5,11) = 1187.109584;

  wo_(6,7) = -2942.77456;
  wo_(6,8) = -242361.5472;
  wo_(6,9) = -248343.412;
  wo_(6,10) = -154666.5808;
  wo_(6,11) = -353880.628;

  wo_(7,8) = -11757.4584;
  wo_(7,9) = 2925.130632;
  wo_(7,10) = 3264.1476;
  wo_(7,11) = -254.0696344;

  wo_(8,9) = -330220.108;
  wo_(8,10) = -387486.0976;
  wo_(8,11) = -188961.5736;

  wo_(9,10) = -262671.1016;
  wo_(9,11) = -116767.072;

  wo_(10,11) = -75854.6648;

  for(int i=0; i < 12; i++) 
    {
      for(int j = i+1; j <12; j++)
	{
	  wo_(j,i) = wo_(i,j);
	}
    }

  wij_(0,0) = -0.3409307328e5;
  wij_(1,0) =  0.0000000000e0; 
  wij_(2,0) = -0.1891165069e6;
  wij_(3,0) =  0.1359353090e6;
  wij_(4,0) = -0.1957509655e6;
  wij_(5,0) =  0.0000000000e0;
  wij_(6,0) = -0.8641805096e5; 
  wij_(7,0) = -0.2099966494e6;
  wij_(8,0) = -0.3222531515e6;
  wij_(9,0) = -0.3497975003e6;
  for(int i=10; i<20; i++) 
    { 
      wij_(i,0) = 0.;
    }

  wij_(0,1) = -0.5996216205e5;
  wij_(1,1) =  0.0000000000e0;
  wij_(2,1) = -0.5909571939e6;
  wij_(3,1) =  0.4469622832e7;
  wij_(4,1) =  0.2166575511e5;
  wij_(5,1) =  0.0000000000e0;
  wij_(6,1) =  0.5286613028e5;
  wij_(7,1) = -0.3287921030e6;
  wij_(8,1) =  0.1400344548e6;
  wij_(9,1) = 0.3090699699e6;
  wij_(10,1) = 0.6049208301e4;
  wij_(11,1) = 0.0000000000e0;
  wij_(12,1) = 0.4139537240e5;
  wij_(13,1) = -0.5293012049e6;
  wij_(14,1) =  0.1213717202e4;
  wij_(15,1) =  0.0000000000e0;
  wij_(16,1) = -0.1344620202e5;
  wij_(17,1) =  0.1278883196e5;
  wij_(18,1) = -0.3521319385e5;
  wij_(19,1) = -0.5800953852e5;

} 

void Solwcad::exsolution(const double &p, const double &T, const double &h2o,
			 const double &co2, const valarray<double> &magma)
{  
#ifdef DEBUG
  cout << "Solwcad::exsolution(), start " << endl;
#endif
  double gas_tot;
  const double tresho = 5.e-5;
  valarray<double> ctrm2(0.0,2);
  valarray<double> fmol(0.0,12), f(0.0,12);
  valarray<double> oxide(0.0,12), x_fraction(0.0,12);

  valarray<double> gas_oxide(2);
  valarray<double> gas_scaling(2);

  valarray<double> composition(12);

	
  composition[0] = h2o;
  composition[1] = co2;
  for(int i=2; i<12; i++)
    {
      composition[i] = magma[i-2];
    }
     
  p_ = p;
  T_ = T;
  reexpressComposition(composition, oxide, x_fraction, fmol);

  if(x_fraction[0] < tresho)  
    {
      oxide[0] = 0.;
      exsolutionCO2(p,T,oxide);
      return;
    }
  if(x_fraction[1] < tresho)  
    {
      oxide[1] = 0.;
      exsolutionH2O(p,T,oxide);
      return;
    }
  
  for(int i=0; i<2; i++) 
    {
      gas_fraction_[i] = x_fraction[i]; 
      gas_scaling[i] = gas_fraction_[i]; 
      gas_oxide[i] = oxide[i];  
    }
  gas_tot = oxide[0]+oxide[1]; 

  double p_log = log(p_/p0_co2_);

  //set variables to zero
	
  for(int i = 0; i<2; i++) 
    {
      ctrm_[i] = 0;
    }

  xwo_ = 0;
  //
  for(int i = 2; i<11; i++) 
    {
      ctrm_[0] += x_fraction[i]*wij_(i-2,0);
      ctrm_[1] += x_fraction[i]*(wij_(i-2,1)+p_log*wij_(i+8,1));
      ctrm2[1] += x_fraction[i]*wij_(i+8,1);
      for(int j=i+1; j<12; j++) 
	{
	  xwo_ += x_fraction[i]*x_fraction[j]*wo_(i,j);
	}
    }
  ctrm_[0] += x_fraction[11]*wij_(9,0);
  ctrm_[1] += x_fraction[11]*(wij_(9,1)+p_log*wij_(19,1));
  ctrm2[1] += x_fraction[11]*wij_(19,1);
    
  for(int i=0; i<2; i++) 
    {
      xtrm_[i] = (ctrm_[i]-xwo_)/(R_*T_);
    }

  valarray<double> ref_term(0.0, 2), phi(0.0, 2), v_ref(0.0, 2);
  referenceTermsH2O(ref_term, phi);
  rtrm_[0] = fol_[0]+ref_term[0]/(R_*T_);
  v_ref[0] = vr_;

  referenceTermsCO2(ref_term, phi);
  rtrm_[1] = fol_[1]+ref_term[1]/(R_*T_);
  v_ref[1] = vr_;

  //set initial estimate for one gas
  const double sg = 4.1e-6;
  valarray<double> oxide_1g(12), new_oxide_1g(12);
  valarray<double>  x_fraction_1g(12), gas_est(2);
  oxide_1g[0] = sg*pow(p_,0.5); 
  oxide_1g[1] = 0.;

  for(int i = 2; i<12; i++) 
    {
      oxide_1g[i] = oxide[i];
    }
  reexpressComposition(oxide_1g, new_oxide_1g, x_fraction_1g,
		       fmol);
  gas_est[0] = x_fraction_1g[0];
  gas_est[1] = x_fraction_1g[0]*0.1;
  //start iteration
  int errore = 0;
  double ah, bh, yh, xg;
  double eps = 1.e-12; 
  ah = 0.8*gas_est[0];
  bh = 1.2*gas_est[0];

  SolwcadFunction2 solwcad_function2(*this, phi[0], xtrm_[0], rtrm_[0]);
  do
    {
      if(bh > 1.) bh = 1.;
      errore = nr_.secantMethod( solwcad_function2, ah, bh, eps, xg);
      if(errore == 0)
	{
	  ah = 0.5*ah;
	  bh = 1.5*bh;
	}
      if(ah < 1.e-15 && bh > 1.)
	{
	  cout << "error from Solwcad::exsolution" << endl;
	  cout << "at P = " << p_ << " and T = "<< T_<< endl;
	  cout << "the solution cannot be bracketed" << endl;
	  exit(1);
	}
    }while(errore == 0);
  //

  valarray<double> molar_fraction(0.0,4); 
  for(int i=0; i<2; i++) 
    {
      x_fraction[i] = 0.;
    }
     
  // set initial estimate for two gas 
  gas_est[0] = xg;
  gas_est[1] = xg*0.1;

  for(int i=0; i<2; i++) 
    {
      if(gas_est[i]>= gas_fraction_[i]) 
	{
	  gas_est[i] = gas_fraction_[i]-1e-6;
	}
    }

  //calculate for two gas
  //scaling unknowns 
  int kspur = 0;
  int esci = 0; 
  const double tol=1.e-1;
  bool check; 
  valarray<double> xgg(0.0,2), xgg1(0.0,2);

  scalingUnknowns (gas_est,gas_scaling);
  for(int i=0; i<2; i++)
    {
      xgg[i] = gas_est[i];
      xgg1[i] = xgg[i];
    } 
  SolwcadNewt solwcad_newt(*this);
  valarray<double> fvec(0.0,2); 
  nr_.newt(xgg, check, solwcad_newt); 
  vecfunc(xgg, fvec);
  if(check) cout << "convergency problems in newt" << endl;
       
  do
    {
      for(int i = 0; i< 2; i++)
	{
	  if(fabs(fvec[i]) > tol)
	    {
	      kspur += 1;
	      for(int j=0; j<2; j++) 
		{
		  xgg[j] = xgg1[j]*0.8/kspur;
		} 
	      // non abbasso kspur MELI: va bene?
	      if(kspur < 3)
		{
		  scalingUnknowns (gas_est, gas_scaling);
		  nr_.newt(xgg, check, solwcad_newt); 
		  vecfunc(xgg, fvec);
		  if(check) cout << "convergency problems in newt" << endl;
		}
	      if(kspur >= 3) 
		{
		  for(int ispur=0; ispur<2; ispur++)
		    {
		      xgg[ispur] = gas_fraction_[ispur];
		    }
		}
	    }
	  if(fabs(fvec[i]) > tol && kspur < 3)
	    {
	      esci = 0; 
	      break; 
	    }
	  xgg[i] = xgg[i]*scale_factor_[i];
	  x_fraction[i] = xgg[i];
	  molar_fraction[i] = xgg[i];
	  esci = 1; 
	}
    }while( esci == 0); 
  molar_fraction[2] = y_[0];
  molar_fraction[3] = 1 - molar_fraction[2];

  //set properties and weight_fraction_[1->4)
  double xvolo = 1.-(x_fraction[0]+x_fraction[1]);
  double xvolo2 = xvolo*xvolo;
  valarray<double> act(0.0, 2), gam(0.0,2), phig(0.0,2), vexc(0.0,2);
  //H2O
  vexc[0] = x_fraction[1]*(-xvolo*ctrm2[1]+
			   (1.-x_fraction[0])*h2o_co2_[1])/p_;
  gam[0] = (1.-x_fraction[0])*(xvolo*ctrm_[0]+
			       x_fraction[1]*(h2o_co2_[0]+p_log*h2o_co2_[1]))-
    x_fraction[1]*xvolo*ctrm_[1]-xvolo2*xwo_;
  gam[0] = exp(gam[0]/(R_*T_));
  act[0] = gam[0]*x_fraction[0];
  phig[0] = phi[0];
  //CO2 
  vexc[1] = (1.-x_fraction[1])*(xvolo*ctrm2[1]+
				x_fraction[0]*h2o_co2_[1])/p_;
  gam[1] = (1.-x_fraction[1])*(xvolo*ctrm_[1]+
			       x_fraction[0]*(h2o_co2_[0]+p_log*h2o_co2_[1]))-
    x_fraction[0]*xvolo*ctrm_[0]-xvolo2*xwo_;
  gam[1] = exp(gam[1]/(R_*T));
  act[1] = gam[1]*x_fraction[1];
  phig[1] = phi[1];
  //

  for(int i=0; i<2; i++)
    {
      gas_est[i] = x_fraction[i];
    }
  double xyws = molar_fraction[2]*molecular_weight_[0]+
    molar_fraction[3]*molecular_weight_[1];
  for(int i=0; i<12; i++)
    {
      if(i >= 2) x_fraction[i] *= xvolo;
      fmol[i] = x_fraction[i]/ghiorso_par_[i];
    }
  fmol[2]=fmol[2]+0.5*(fmol[6]+fmol[7]+fmol[8]+fmol[9])
    +fmol[10]+fmol[11];

  double ox_sum = 0.;
  for(int i=0; i < 12; i++)
    {
      double oxn = fmol[i]*molecular_weight_[i];
      ox_sum += oxn;
      if(i <= 1) weight_fraction_[i] = oxn;
    }
  weight_fraction_[0] = weight_fraction_[0]/ox_sum;
  weight_fraction_[1] = weight_fraction_[1]/ox_sum;
  weight_fraction_[2] = molar_fraction[2]*molecular_weight_[0]/xyws;
  weight_fraction_[3] =1.-weight_fraction_[2];

  //check mass conservation
  wtfrex_ = (gas_tot-weight_fraction_[0]-weight_fraction_[1])/
    (1.-weight_fraction_[0]-weight_fraction_[1]);

  //undersatured conditions
  if(wtfrex_ <= 1.e-16) 
    {
      wtfrex_ = 0.;

      weight_fraction_[0] = gas_oxide[0];
      weight_fraction_[1] = gas_oxide[1];
      weight_fraction_[2] = 0.;
      weight_fraction_[3] = 0.;
    } 

}

  
void Solwcad::exsolutionCO2(const double &p, const double &T,
			    const valarray<double> &composition) 
{
#ifdef DEBUG
  cout << "Solwcad::exsolutionCO2(), start " << endl;
#endif
  double gas_tot;
  valarray<double> ctrm2(0.0,2);
  valarray<double> fmol(0.0,12), f(0.0,12);
  valarray<double> oxide(0.0,12), x_fraction(0.0, 12);
  valarray<double> gas_oxide(2);
  valarray<double> gas_scaling(2);
    
  p_ = p;
  T_ = T;
    
  reexpressComposition(composition, oxide, x_fraction, fmol);
  for(int i=0; i<2; i++) 
    {
      gas_fraction_[i] = x_fraction[i];  
      gas_scaling[i] = gas_fraction_[i]; 
      gas_oxide[i] = oxide[i]; 
    }
  gas_tot = oxide[2];    
  double p_log = log(p_/p0_co2_);

  //set to zero variables
  for(int i = 0; i<2; i++) 
    {
      ctrm_[i] = 0;
    }
  xwo_ = 0;
  //
  for(int i = 2; i<11; i++) 
    {
      ctrm_[0] += x_fraction[i]*wij_(i-2,0);
      ctrm_[1] += x_fraction[i]*(wij_(i-2,1)+p_log*wij_(i+8,1));
      ctrm2[1] += x_fraction[i]*wij_(i+8,1);
      for(int j=i+1; j<12; j++) 
	{
	  xwo_ += x_fraction[i]*x_fraction[j]*wo_(i,j);
	}
    }
  ctrm_[0] += x_fraction[11]*wij_(9,0);
  ctrm_[1] += x_fraction[11]*(wij_(9,1)+p_log*wij_(19,1));
  ctrm2[1] += x_fraction[11]*wij_(19,1);
    
  for(int i=0; i<2; i++) 
    {
      xtrm_[i] = (ctrm_[i]-xwo_)/(R_*T_);
    }
    
  valarray<double> ref_term(0.0,2), phi(0.0,2), v_ref(0.0,2);

  referenceTermsCO2(ref_term, phi);
  rtrm_[1] = fol_[1]+ref_term[1]/(R_*T_);
  v_ref[1] = vr_;

  //set initial estimate for one gas
  const double sg = 4.1e-6;
  valarray<double> oxide_1g(12), new_oxide_1g(0.0,12);
  valarray<double>  x_fraction_1g(0.0,12), gas_est(2);

  oxide_1g[0] = sg*pow(p_,0.5); 
  oxide_1g[1] = 0.;
  for(int i = 2; i<12; i++) 
    {
      oxide_1g[i] = oxide[i];
    }  

  reexpressComposition(oxide_1g, new_oxide_1g, x_fraction_1g,fmol);
  gas_est[0] = x_fraction_1g[0];
  gas_est[1] = x_fraction_1g[0]*0.1;
 
  //start iteration
  int errore = 0;
  double ah, bh, yh, xg;
  double eps = 1.e-12; 
  ah = 0.8*gas_est[1];
  bh = 1.2*gas_est[1];
 
  SolwcadFunction2 solwcad_function2(*this, phi[1], xtrm_[1], rtrm_[1]);
  do
    { 
      if(bh >= 1.) bh = 1.;
      errore = nr_.secantMethod( solwcad_function2, ah, bh, eps, xg);
      if(errore == 0)
	{
	  ah = 0.5*ah;
	  bh = 1.5*bh;
	}
      if(ah < 1.e-15 && bh > 1.)
	{
	  cout << "error from Solwcad::exsolution" << endl;
	  cout << "at P = " << p_ << " and T = "<< T_<< endl;
	  cout << "the solution cannot be bracketed" << endl;
	  exit(1); 
	} 
    }while(errore == 0);
  //  
  
  valarray<double> molar_fraction(0.0,4); //xym(i)  
  if(xg < gas_fraction_[1])
    {
      x_fraction[1] = xg; 
      molar_fraction[1] = xg;
      molar_fraction[3] = 1;
    }
  else
    {
      x_fraction[1] = gas_fraction_[1];
      molar_fraction[1] = gas_fraction_[1];
    }

  //set properties and weight_fraction_[1->4)
  double xvolo = 1.-(x_fraction[0]+x_fraction[1]);
  double xvolo2=xvolo*xvolo;
  valarray<double> act(0.0,2), gam(0.0,2), phig(0.0,2), vexc(0.0,2);
  //CO2
  vexc[1] = (1.-x_fraction[1])*(xvolo*ctrm2[1]+
				x_fraction[0]*h2o_co2_[1])/p_;
  gam[1] = (1.-x_fraction[1])*(xvolo*ctrm_[1]+
			       x_fraction[0]*(h2o_co2_[0]+p_log*h2o_co2_[1]))-
    x_fraction[0]*xvolo*ctrm_[0]-xvolo2*xwo_;
  gam[1] = exp(gam[1]/(R_*T));
  act[1] = gam[1]*x_fraction[1];
  phig[1] = phi[1];
  //
  for(int i=0; i<2; i++)
    {
      gas_est[i] = x_fraction[i];
    }
  double xyws = molar_fraction[2]*molecular_weight_[0]+
    molar_fraction[3]*molecular_weight_[1];
  for(int i=0; i<12; i++)
    {
      if(i >= 2) x_fraction[i] *= xvolo;
      fmol[i] = x_fraction[i]/ghiorso_par_[i];
    }
  fmol[2]=fmol[2]+0.5*(fmol[6]+fmol[7]+fmol[8]+fmol[9])
    +fmol[10]+fmol[11];

  double ox_sum = 0.;
  for(int i=0; i < 12; i++)
    {
      double oxn = fmol[i]*molecular_weight_[i];
      ox_sum += oxn;
      if(i <= 1) weight_fraction_[i] = oxn;
    }
  weight_fraction_[0] = weight_fraction_[0]/ox_sum;
  weight_fraction_[1] = weight_fraction_[1]/ox_sum;
  weight_fraction_[2] = molar_fraction[2]*molecular_weight_[0]/xyws;
  weight_fraction_[3] =1.-weight_fraction_[2];

  //check mass conservation
  wtfrex_ = (gas_tot-weight_fraction_[0]-weight_fraction_[1])/
    (1.-weight_fraction_[0]-weight_fraction_[1]);
 
  //undersatured conditions
  if(wtfrex_ <= 1.e-16) wtfrex_ = 0.;
    
}
  

void Solwcad::exsolutionH2O(const double &p, const double &T,
			    const valarray<double> &composition)
{
#ifdef DEBUG
  cout << "Solwcad::exsolutionH2O(), start " << endl;
#endif
  double gas_tot;
  valarray<double> ctrm2(0.0,2);
  valarray<double> fmol(0.0,12), f(0.0, 12);
  valarray<double> oxide(0.0, 12), x_fraction(0.0, 12);
  valarray<double> gas_oxide(2);
  valarray<double> gas_scaling(2);
     
  p_ = p;
  T_ = T;

  reexpressComposition(composition, oxide, x_fraction, fmol);
  for(int i=0; i<2; i++) 
    {
      gas_fraction_[i] = x_fraction[i];  
      gas_scaling[i] = gas_fraction_[i]; 
      gas_oxide[i] = oxide[i]; 
    }
  gas_tot = oxide[0];    
  double p_log = log(p_/p0_co2_);

  //set to zero variables
  for(int i = 0; i<2; i++) 
    {
      ctrm_[i] = 0;
    }
  xwo_ = 0;
  //
  for(int i = 2; i<11; i++) 
    {
      ctrm_[0] += x_fraction[i]*wij_(i-2,0);
      ctrm_[1] += x_fraction[i]*(wij_(i-2,1)+p_log*wij_(i+8,1));
      ctrm2[1] += x_fraction[i]*wij_(i+8,1);
      for(int j=i+1; j<12; j++) 
	{
	  xwo_ += x_fraction[i]*x_fraction[j]*wo_(i,j);
	}
    }
  ctrm_[0] += x_fraction[11]*wij_(9,0);
  ctrm_[1] += x_fraction[11]*(wij_(9,1)+p_log*wij_(19,1));
  ctrm2[1] += x_fraction[11]*wij_(19,1);
  for(int i=0; i<2; i++) 
    {
      xtrm_[i] = (ctrm_[i]-xwo_)/(R_*T_);
    }

  valarray<double> ref_term(0.0,2), phi(0.0,2), v_ref(0.0,2);
    
  referenceTermsH2O(ref_term, phi);
  rtrm_[0] = fol_[0]+ref_term[0]/(R_*T_);
  v_ref[0] = vr_;

  //set initial estimate for one gas
  const double sg = 4.1e-6;
  valarray<double> oxide_1g(12), new_oxide_1g(0.0,12);
  valarray<double>  x_fraction_1g(0.0,12), gas_est(2);

  oxide_1g[0] = sg*pow(p_,0.5); 
  oxide_1g[1] = 0.;
  for(int i = 2; i<12; i++) 
    {
      oxide_1g[i] = oxide[i];
    }

  reexpressComposition(oxide_1g, new_oxide_1g, x_fraction_1g,fmol);
  gas_est[0] = x_fraction_1g[0];
  gas_est[1] = x_fraction_1g[0]*0.1;

  //start iteration
  int errore = 0;
  double ah, bh, yh, xg;
  double eps = 1.e-12; 
  ah = 0.8*gas_est[0];
  bh = 1.2*gas_est[0];

  SolwcadFunction2 solwcad_function2(*this, phi[0], xtrm_[0], rtrm_[0]);
  do
    {
      if(bh >= 1.) bh = 1.;
      errore = nr_.secantMethod( solwcad_function2, ah, bh, eps, xg);
      if(errore == 0)
	{
	  ah = 0.5*ah;
	  bh = 1.5*bh;
	}
      if(ah < 1.e-15 && bh > 1.)
	{
	  cout << "error from Solwcad::exsolution" << endl;
	  cout << "at P = " << p_ << " and T = "<< T_<< endl;
	  cout << "the solution cannot be bracketed" << endl;
	  exit(1);
	}
    }while(errore == 0);
  //

  valarray<double> molar_fraction(0.0,4); 
  if(xg < gas_fraction_[0])
    {
      x_fraction[0] = xg;
      molar_fraction[0] = xg;
      molar_fraction[2] = 1;
    }
  else
    {
      x_fraction[0] = gas_fraction_[0];
      molar_fraction[0] = gas_fraction_[0];
    }

  //set properties and weight_fraction_[1->4)
  double xvolo = 1.-(x_fraction[0]+x_fraction[1]);
  double xvolo2=xvolo*xvolo;
  valarray<double> act(0.0, 2 ), gam(0.0, 2), phig(0.0, 2), vexc(0.0,2);
  //H2O
  vexc[0]= x_fraction[1]*(-xvolo*ctrm2[1]+
			  (1.-x_fraction[0])*h2o_co2_[1])/p_;
  gam[0] = (1.-x_fraction[0])*(xvolo*ctrm_[0]+
			       x_fraction[1]*(h2o_co2_[0]+p_log*h2o_co2_[1]))-
    x_fraction[1]*xvolo*ctrm_[1]-xvolo2*xwo_;
  gam[0] = exp(gam[0]/(R_*T));
  act[0] = gam[0]*x_fraction[0];
  phig[0] = phi[0];
  //
  for(int i=0; i<2; i++)
    {
      gas_est[i] = x_fraction[i];
    }
  double xyws = molar_fraction[2]*molecular_weight_[0]+
    molar_fraction[3]*molecular_weight_[1];
  for(int i=0; i<12; i++)
    {
      if(i >= 2) x_fraction[i] *= xvolo;
      fmol[i] = x_fraction[i]/ghiorso_par_[i];
    } 
  fmol[2]=fmol[2]+0.5*(fmol[6]+fmol[7]+fmol[8]+fmol[9])
    +fmol[10]+fmol[11];

  double ox_sum = 0.;  
  for(int i=0; i < 12; i++)
    {
      double oxn = fmol[i]*molecular_weight_[i];
      ox_sum += oxn;
      if(i <= 1) weight_fraction_[i] = oxn;
    }
  weight_fraction_[0] = weight_fraction_[0]/ox_sum;
  weight_fraction_[1] = weight_fraction_[1]/ox_sum;
  weight_fraction_[2] = molar_fraction[2]*molecular_weight_[0]/xyws;
  weight_fraction_[3] =1.-weight_fraction_[2];

  //check mass conservation
  wtfrex_ = (gas_tot-weight_fraction_[0]-weight_fraction_[1])/
    (1.-weight_fraction_[0]-weight_fraction_[1]);
 
  //undersatured conditions
  if(wtfrex_ <= 1.e-16) wtfrex_ = 0.;
    
}
 

double Solwcad::funz(const double &x, const double &bm, const double &cm,
		     const double &dm, const double &em, double &ym)
{
  double a = cm+dm/x+em/(x*x);
  ym = 0.25*bm/x;
  double y = p_-R_*T_*(1.+ym+(ym*ym)-(ym*ym*ym))/
    (x*((1.-ym)*(1.-ym)*(1.-ym)))+a/(sqrt(T_)*x*(x+bm));
  return y;
} 
   
 
double Solwcad::gasCalc(double &x, const double &phi1, 
			const double &xtrm1, const double &rtrm1)
{
  if(x <= 1e-15)
    {
      x = 1e-6;
    }
  double y = log(phi1*p_/x)-(1.-x)*(1-x)*xtrm1-rtrm1;
  return y;
}

void Solwcad::reexpressComposition(const valarray<double> &oxide,
				   valarray<double> &new_oxide, 
				   valarray<double> &x_fraction,
				   valarray<double> &fmol)
{ 
#ifdef DEBUG
  cout << "start reexpresscomposition "<<endl;
#endif
  double fraction_sum = 0.; 
  double oxide_sum = 0.; 

  new_oxide[0] = oxide[0];
  new_oxide[1] = oxide[1]; 
  for(int i=2; i<12; i++) 
    {
      oxide_sum += oxide[i];
    }
  for(int i=0; i<12; i++) 
    {
      if(i>= 2)
	{
	  new_oxide[i] = oxide[i]*(1-oxide[0]-oxide[1])/oxide_sum;
	} 
      fmol[i] = new_oxide[i]/molecular_weight_[i];
    }
  fmol[2] = fmol[2] -0.5*(fmol[6]+fmol[7]+fmol[8]+fmol[9])
    -fmol[10]-fmol[11];
  for(int i=0; i<12; i++) 
    {
      x_fraction[i] = fmol[i]*ghiorso_par_[i];
      fraction_sum += x_fraction[i];
    }
  for(int i=0; i<12; i++) 
    {
      x_fraction[i] = x_fraction[i]/fraction_sum;
      if(i>=2) 
	{
	  x_fraction[i] = x_fraction[i]/(1-x_fraction[0]-x_fraction[1]);
	}
    }
}

void Solwcad::referenceTerms(const double &y_cd, 
			     valarray<double> &ref_term,
			     valarray<double> &phi) 
{

  //ref_term = liquid components reference terms
  //phi = gas components fugacities

#ifdef DEBUG
  cout << "start referenceTerms "<<endl;
#endif
  double bm, cm, dm, em, ym;
  Matrix<double> c(2,2), d(2,2), e(2,2);
  valarray<double> b(2);
 
  double v = vm(y_cd, bm, cm, dm, em, ym, b, c, d, e);
  vr_ = 0. ;

  double vln = log((v+bm)/v);
  double tp = pow(T_,1.5);
  double pb1 = R_*tp*bm;
  double pb2 = R_*tp*(bm*bm);
  double pb3 = R_*tp*(bm*bm*bm);
  double pb4 = R_*tp*(bm*bm*bm*bm);
  double pv1 = R_*tp*v;
  double pv2 = R_*tp*v*v;
  double bv2 = 2.*pb1*v*v;
  double vb = v+bm;
  double ym2=ym*ym;
  double ym3=ym*ym*ym;
  double ym_menus2=(1.-ym)*(1.-ym);
  double ym_menus3=(1.-ym)*(1.-ym)*(1.-ym);

  valarray<double> x(2);
  x[0] = 1.-y_cd;
  x[1] = y_cd;
    
  for(int i=0; i<2; i++) 
    {
      double bc1 = b[i]*cm;
      double bd1 = b[i]*dm;
      double be1 = b[i]*em;
      double c1 = 2.*c(i,i)*x[i]+2.*(1.-x[i])*c(0,1);
      double d1 = 2.*d(i,i)*x[i]+2.*(1.-x[i])*d(0,1)+dm;
      double e1 = 2.*(e(i,i)*x[i]+e(0,1)*(1.-x[i])+em);
	
      double phexp;
      phexp = (4.*ym-3.*ym2)/(ym_menus2)+
	b[i]/bm*(4.*ym-2.*ym2)/(ym_menus3)-
	c1/pb1*vln-bc1/(pb1*vb)+
	bc1/pb2*vln-d1/(pb1*v)+d1/pb2*vln+
	bd1/(pb1*v*vb)+2.*bd1/(pb2*vb)-2.*bd1/pb3*vln-
	e1/bv2+e1/(pb2*v)-e1/pb3*vln+
	be1/(bv2*vb)-3.*be1/(2.*pb2*v*vb)+
	3.*be1/pb4*vln-3.*be1/(pb3*vb)-log(p_*v/(R_*T_));
      phi[i] = exp(phexp);
    };
#ifdef DEBUG
  cout << "exit referenceTerms" << endl;
#endif
} 
  
void Solwcad::referenceTermsCO2(valarray<double> &ref_term,
				valarray<double> &phi) 
{
  //ref_term = liquid components reference terms
  //phi = gas components fugacities

#ifdef DEBUG
  cout << "start referenceTermsCO2 "<<endl;
#endif
  double bm, cm, dm, em, ym;
  Matrix<double> c(2,2), d(2,2), e(2,2);
  valarray<double> b(2);

  double v = vmCO2(bm, cm, dm, em, ym, b, c, d, e);
  double T2=T_*T_;
  double p2=p_*p_;
  double p3=p_*p_*p_;

  vr_ = 0. ;
  vr_ = par_[0] +par_[1] *T_+par_[2] *T2+par_[3] *T_*T_*T_+
    p_*(par_[4] +par_[5] *T_+par_[6] *T2)+p2*
    (par_[7] +par_[8] *T_)+p3*par_[9] ;
    
  double vln = log((v+bm)/v);
    
  double t1 = 2.-T_/to_[1];
  double t2 = 2.*T_-to_[1];
  double t3 = 2.*T2-(to_[1]*to_[1]);
  double tp = pow(T_,1.5);
  double pb1 = R_*tp*bm;
  double pb2 = R_*tp*bm*bm;
  double pb3 = R_*tp*bm*bm*bm;
  double pb4 = R_*tp*bm*bm*bm*bm;
  double pv1 = R_*tp*v;
  double pv2 = R_*tp*v*v;
  double bv2 = 2.*pb1*v*v;
  double vb = v+bm;
  double ym2=ym*ym;
  double ym3=ym*ym*ym;
  double ym_menus3=(1.-ym)*(1.-ym)*(1.-ym);

  double phexp;
  phexp = (8.*ym-9.*ym2+3.*ym3)/ym_menus3-
    cm/(R_*tp*vb)-dm/(pv1*vb)-em/(pv2*vb)-
    cm/pb1*vln-dm/(pv1*bm)+dm/pb2*vln-
    em/bv2+em/(pb2*v)-em/pb3*vln-log(p_*v/(R_*T_));
  phi[1] = exp(phexp);

  ref_term[1] = p_*(par_[0] *t1+par_[1] *T_+par_[2] *T_*t2+par_[3] *T_*t3)+
    p2*(par_[4] *t1+par_[5] *T_+par_[6] *T_*t2)*0.5+
    p3*(par_[7] *t1+par_[8] *T_)/3.+
    p_*p_*p_*p_*par_[9] *t1*0.25;
} 
  
void Solwcad::referenceTermsH2O(valarray<double> &ref_term,
				valarray<double> &phi) 
{
  //ref_term = liquid components reference terms
  //phi = gas components fugacities

#ifdef DEBUG
  cout << "start referenceTermsH2O "<<endl;
#endif
  double bm, cm, dm, em, ym;
  Matrix<double> c(2,2), d(2,2), e(2,2);
  valarray<double> b(2);
  valarray<double> bd(10); 
    
  bd[0] = 9.144e-6;    
  bd[1] = 3.685e-9;     
  bd[2] = 1.168e-11;    
  bd[3] = -1.99e-15;    
  bd[4] = 1.22e-16;    
  bd[5] = -1.945e-17;    
  bd[6] = -1.58e-21;    
  bd[7] = 4.68e-24;    
  bd[8] = 1.144e-26;    
  bd[9] =-3.96e-33 ;    
    
  double T2=T_*T_;
  double p2=p_*p_;
  double p3=p_*p_*p_;
  double v = vmH2O(bm, cm, dm, em, ym, b, c, d, e);
  vr_ = bd[0]+bd[1]*T_+bd[2]*T2+bd[3]*T_*T_*T_+
    p_*(bd[4]+bd[5]*T_+bd[6]*T2)+p2*(bd[7]+bd[8]*T_)+
    p3*bd[9];

  double vln = log((v+bm)/v);

  double t1 = 2.-T_/to_[0];
  double t2 = 2.*T_-to_[0];
  double t3 = 2.*T2-(to_[0]*to_[0]);
  double tp = pow(T_,1.5);
  double pb1 = R_*tp*bm;
  double pb2 = R_*tp*bm*bm;
  double pb3 = R_*tp*bm*bm*bm;
  double pv1 = R_*tp*v;
  double pv2 = R_*tp*v*v;
  double bv2 = 2.*pb1*v*v;
  double vb = v+bm;
  double ym2=ym*ym;
  double ym3=ym*ym*ym;
  double ym_menus3=(1.-ym)*(1.-ym)*(1.-ym);

  double phexp;
  phexp = (8.*ym-9.*ym2+3.*ym3)/ym_menus3-
    cm/(R_*tp*vb)-dm/(pv1*vb)-em/(pv2*vb)-
    cm/pb1*vln-dm/(pv1*bm)+dm/pb2*vln-
    em/bv2+em/(pb2*v)-em/pb3*vln-log(p_*v/(R_*T_));
  phi[0] = exp(phexp);
  ref_term[0] = p_*(bd[0]*t1+bd[1]*T_+bd[2]*T_*t2+bd[3]*T_*t3)+
    p2*(bd[4]*t1+bd[5]*T_+bd[6]*T_*t2)*0.5+
    p3*(bd[7]*t1+bd[8]*T_)/3.+
    p_*p_*p_*p_*bd[9]*t1*0.25;    
} 
  

void Solwcad::scalingUnknowns (valarray<double> &gas_est,
			       valarray<double> &gas_scaling)
{
#ifdef DEBUG
  cout << "start scalingUnknowns "<<endl;
#endif
  for(int i=0; i<2; i++) 
    { 
      const double tolscl = 0.9;
      int inot;
      if(i==0) inot=1;
      if(i==1) inot=0;
      if(gas_est[i]/gas_est[inot] < tolscl) 
	{
	  scale_factor_[i] = gas_est[i]/gas_est[inot];		
	  gas_scaling[i] = gas_scaling[i]/scale_factor_[i] ;
	  gas_est[i] = gas_est[inot];
	} 
    }
}
  
void Solwcad::vecfunc(valarray<double> &x, valarray<double> &fv) 
{
#ifdef DEBUG
  cout << "start vecfunc "<<endl;
#endif
  const double tin = 1.e-16;
  valarray<double> x2(2);
 
  for(int i=0; i<2; i++) 
    { 
      if(x[i] <= 0.) x[i] = tin;
      x2[i] = x[i]*scale_factor_[i];// MELI da controllare
    }

  y_[0] =(gas_fraction_[0]*(1.-x2[1])-x2[0]*(1.-gas_fraction_[1]))/
    (gas_fraction_[1]-x2[1]+gas_fraction_[0]-x2[0]);
  if(y_[0] <= tin) 
    {
      y_[0] = tin;
    }
  if(y_[0] >= (1.-tin)) 
    {
      y_[0] = 1.-tin;
    }
  y_[1] = 1-y_[0];
 
  int inot; 
  for(int i=0; i<2; i++)  
    { 
      if(i == 1) inot = 0;
      if(i == 0) inot = 1;
	
      if(y_[i] <= 0.)   
	{
	  y_[i] = tin;  
	  y_[inot] = 1.-tin;
	}
      if(y_[i] >= 1.) 
	{
	  y_[i] = 1.-tin;
	  y_[inot] = 1.-y_[i];
	}
    }

  valarray<double> f_ref(0.0,2), phi(0.0,2), gtrm(2);

  referenceTerms(y_[1], f_ref, phi);
    
  for(int i=0; i<2; i++) 
    {
      if(i == 1) inot = 0;
      if(i == 0) inot = 1;
      gtrm[i] = log(phi[i]*y_[i]*p_/x2[i]);
      xtrm_[i]= (1.-x2[0]-x2[1])*
	((1.-x2[i])*ctrm_[i]-x2[inot]*ctrm_[inot])
	-(1.-x2[0]-x2[1])*(1.-x2[0]-x2[1])*xwo_+(1.-x2[i])*x2[inot]*whc_;
      xtrm_[i] = xtrm_[i]/(R_*T_);
	
      fv[i] = gtrm[i]-xtrm_[i]-rtrm_[i];
    }
#ifdef DEBUG
  cout << "exit da vecfunc" << endl;
#endif
} 
  
 
double Solwcad::vm(const double &y_cd,double &bm,double &cm,
		   double &dm,double &em,double &ym, 
		   valarray<double> &b, Matrix<double> &c, 
		   Matrix<double> &d, Matrix<double> &e) 
{
#ifdef DEBUG
  cout << "start vm "<<endl;
#endif
  double ah = 0;
  double bh, v;
  double T2=T_*T_;
  valarray<double> x(2);

  x[0] = 1.-y_cd;
  x[1] = y_cd;
  b[0] = 2.9e-5;
  b[1] = 5.8e-5;
  c(0,0) = (290.78-0.30276*T_+1.4774e-4*T2)*1e-1;
  d(0,0) = (-8374.+19.437*T_-8.148e-3*T2)*1e-7;
  e(0,0) = (76600.-133.9*T_+0.1071*T2)*1e-13;
  c(1,1) = (28.31+0.10721*T_-8.81e-6*T2)*1e-1;
  d(1,1) = (9380.-8.53*T_+1.189e-3*T2)*1e-7;
  e(1,1) = (-368654.+715.9*T_+0.1534*T2)*1e-13;
  if(c(0,0)*c(1,1) <= 0.) 
    {
      c(0,1) = 0.;
    }
  else 
    {
      c(0,1) = sqrt(c(0,0)*c(1,1));
    }
  if(d(0,0)*d(1,1) <= 0.) 
    {
      d(0,1) = 0.;
    }
  else 
    {
      d(0,1) = sqrt(d(0,0)*d(1,1));
    }
  if(e(0,0)*e(1,1) <= 0.) 
    {
      e(0,1) = 0.;
    }
  else 
    {
      e(0,1) = sqrt(e(0,0)*e(1,1));
    }
  c(1,0) = c(0,1);
  d(1,0) = d(0,1);
  e(1,0) = e(0,1);
  bm = x[0]*b[0]+x[1]*b[1];
  cm = 0.;
  dm = 0.;
  em = 0.;
  for(int i=0; i<2; i++) 
    {
      for(int j=0; j<2; j++) 
	{
	  cm += x[i]*x[j]*c(i,j);
	  dm += x[i]*x[j]*d(i,j);
	  em += x[i]*x[j]*e(i,j);
	}
    }
  ah = 1.9901e-5-4.6125e-15*p_+5.5523e-25*p_*p_+
    1.3539e-5*y_cd;
  if(T_ <= 1300.) ah = ah*0.8;

  //start iteration
  int errore = 0;
  int counter = 0;
  int counter2 = 0;
  double ah_initial = ah;
  double ahfac = 0.1; 
  double eps = 1.e-16; 
  double fac = 0.2; 
  SolwcadFunction1 solwcad_function1(*this, bm, cm, dm, em, ym);

  do
    { 
      bh = ah*(1.+fac);
      errore = nr_.secantMethod( solwcad_function1, ah, bh, eps, v);
	
      //change boundaries
      if(errore == 0)
	{
	  counter++;
	  if(counter == 20)
	    {
	      counter = 0;
	      counter2++;
	      if(counter2 == 100)
		{
		  counter2 = 0;
		  ah = ah_initial;
		  ahfac = -ahfac;
		}
	      ah = ah*(1+ahfac);
	      fac = 0.1;
	    }
	  fac = 2*fac;
	} 
	
    }    while(errore == 0); 
  ym = solwcad_function1.p5_;
  return v;
}
  
double Solwcad::vmCO2(double &bm,double &cm,
		      double &dm,double &em,double &ym, 
		      valarray<double> &b, Matrix<double> &c, 
		      Matrix<double> &d, Matrix<double> &e) 
{
#ifdef DEBUG
  cout << "start vmCO2 "<<endl;
#endif

  double ah, bh, v;
  double T2=T_*T_;

  b[0] = 0.;
  b[1] = 5.8e-5;
  c(0,0) = 0.;
  d(0,0) = 0.;
  e(0,0) = 0.;
  c(1,1) = (28.31+0.10721*T_-8.81e-6*T2)*1e-1;
  d(1,1) = (9380.-8.53*T_+1.189e-3*T2)*1e-7;
  e(1,1) = (-368654.+715.9*T_+0.1534*T2)*1e-13;
  bm = b[1];
  cm = c(1,1);
  dm = d(1,1);
  em = e(1,1);
  ah = 2.e-5;
    
  //start iteration
  int errore = 0; 
  int counter = 0;
  int counter2 = 0;
  double ah_initial = ah;
  double ahfac = 0.1; 
  double eps = 1.e-16; 
  double fac = 0.2;
  SolwcadFunction1 solwcad_function1(*this, bm, cm, dm, em, ym);
  do
    { 
      bh = ah*(1.+fac);
      errore = nr_.secantMethod( solwcad_function1, ah, bh, eps, v);

      //change boundaries
      if(errore == 0)
	{
	  counter++;
	  if(counter == 20)
	    {
	      counter = 0;
	      counter2++;
	      if(counter2 == 100)
		{
		  counter2 = 0;
		  ah = ah_initial;
		  ahfac = -ahfac;
		}
	      ah = ah*(1+ahfac);
	      fac = 0.1;
	    }
	  fac = 2*fac;
	} 
      //
    }    while(errore == 0); 
  ym = solwcad_function1.p5_;
  return v;
}
double Solwcad::vmH2O(double &bm,double &cm,
		      double &dm,double &em,double &ym, 
		      valarray<double> &b, Matrix<double> &c, 
		      Matrix<double> &d, Matrix<double> &e) 
{
#ifdef DEBUG
  cout << "start vmH2O "<<endl;
#endif
  double ah, bh, v;
  double T2=T_*T_;

  b[0] = 2.9e-5;
  b[1] = 0.; 
  c(0,0) = (290.78-0.30276*T_+1.4774e-4*T2)*1e-1;
  d(0,0) = (-8374.+19.437*T_-8.148e-3*T2)*1e-7;
  e(0,0) = (76600.-133.9*T_+0.1071*T2)*1e-13;
  c(1,1) = 0.;
  d(1,1) = 0.;
  e(1,1) = 0.;
  bm = b[0];
  cm = c(0,0);
  dm = d(0,0);
  em = e(0,0);
  ah = 1.9901e-5-4.6125e-15*p_+5.5523e-25*p_*p_;
 
  //start iteration
  int errore = 0;
  int counter = 0;
  int counter2 = 0;
  double ah_initial = ah;
  double ahfac = 0.1; 
  double eps = 1.e-16; 
  double fac = 0.2;
  SolwcadFunction1 solwcad_function1(*this, bm, cm, dm, em, ym);
  do
    { 
      bh = ah*(1.+fac);
      errore = nr_.secantMethod( solwcad_function1, ah, bh, eps, v);
      //change boundaries
      if(errore == 0)
	{
	  counter++;
	  if(counter == 20)
	    {
	      counter = 0;
	      counter2++;
	      if(counter2 == 100)
		{
		  counter2 = 0;
		  ah = ah_initial;
		  ahfac = -ahfac;
		}
	      ah = ah*(1+ahfac);
	      fac = 0.1;
	    }
	  fac = 2*fac;
	}  
      //
    }    while(errore == 0); 
  ym = solwcad_function1.p5_;
  return v;
}
