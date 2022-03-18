#ifndef __SOLWCAD_HPP
#define __SOLWCAD_HPP

#include <valarray>
            
#include "Matrix.hpp"
#include "Function.hpp"
#include "NumericalRecipes.hpp"

//! class to compute weight fraction of the volatiles exsolved and dissolved
class Solwcad {
    
public: 
  Solwcad();
  ~Solwcad() {cout << "~Solwcad()"<<endl;}
  //! calculate the saturation content of H2O and CO2 in silicate liquids, and the composition of the co-existing gas phase. For any given silicate liquid composition in terms of 10 major oxides, pressure, temperature, and total H2O and CO2 in the sistem.
  void exsolution(const double &p, const double &T, 
		  const double &h2o,const double &co2, 
		  const valarray<double> &magma); 
  //! calculate the saturation content of CO2 in silicate liquids, and the composition of the co-existing gas phase. For any given silicate liquid composition in terms of 10 major oxides, pressure, temperature, and total CO2 in the sistem.
  void exsolutionCO2(const double &p, const double &T,
		     const valarray<double> &composition);
  //! calculate the saturation content of H2O in silicate liquids, and the composition of the co-existing gas phase. For any given silicate liquid composition in terms of 10 major oxides, pressure, temperature, and total H2O in the sistem.
  //! INUPTS
  //! composition: 12-elements vector containing total H2O content (mass fraction), total CO2 content (zero) and 10 major oxides mass fractions in this order
  void exsolutionH2O(const double &p, const double &T,
		     const valarray<double> &composition);
  double funz(const double &x,const double &bm, 
	      const double &cm,const double &dm, 
	      const double &em, double &ym);
  double gasCalc(double &x, 
		 const double &phi1, 
		 const double &xtrm1, 
		 const double &rtrm1);
  void vecfunc(valarray<double> &x, 
	       valarray<double> &fv);  
  //! Inspector
  inline const double wf_h2o_l_onliq()const {return weight_fraction_[0];} 
  //! Inspector
  inline const double wf_co2_l_onliq()const {return weight_fraction_[1];} 
  //! Inspector
  inline const double wf_h2o_g_ongas()const {return weight_fraction_[2];} 
  //! Inspector
  inline const double wf_co2_g_ongas()const {return weight_fraction_[3];}
  //! Inspector
  inline const double wf_gas_ontotal()const {return wtfrex_;}
  //! Inspector
  inline const double wf_h2o_l_ontotal()const {return weight_fraction_[0]*(1-wtfrex_);}
  inline const double wf_co2_l_ontotal()const {return weight_fraction_[1]*(1-wtfrex_);}
  //!Debug inspector  
  inline const double getT()const {double temperature=T_;
    return temperature;} 
  //!Debug inspector  
  inline const double getP()const {double pressure=p_;
    return pressure;} 
  //!Debug inspector  
  inline const double getXwo()const {double xwo=xwo_;
    return xwo;} 
  //!Debug inspector  
  inline const valarray<double> getGasFraction()const {
    valarray<double> gas_fraction(2);gas_fraction=gas_fraction_;
    return gas_fraction;} 
  //!Debug inspector  
  inline const valarray<double> getScaleFactor()const {
    valarray<double> scale_factor(2);scale_factor=scale_factor_;
    return scale_factor;} 
  //!Debug inspector  
  inline const valarray<double> getCtrm()const {
    valarray<double> ctrm(2);ctrm=ctrm_;
    return ctrm;} 
  //!Debug inspector  
  inline const valarray<double> getXtrm()const {
    valarray<double> xtrm(2);xtrm=xtrm_;
    return xtrm;} 
  //!Debug inspector  
  inline const valarray<double> getRtrm()const {
    valarray<double> rtrm(2);rtrm=rtrm_;
    return rtrm;} 

protected:  

  //!weight_fraction_(1),(2) = weight fraction of H2O,CO2 on liquid phase
  //!weight_fraction_(3),(4) = weight fraction of H2O,CO2 on gas phase
  valarray<double> weight_fraction_; 
  //! returns the magma composition adjusted for ghiorso method from the oxides weight fraction.
  void reexpressComposition(const valarray<double> &oxide,
			    valarray<double> &new_oxide,
			    valarray<double> &x_fraction,
			    valarray<double> &fmol);
    
  //! return the liquid components reference terms (ref_term_) and gas components fugacities (fugacity_) for CO2 and H2O.
  void referenceTerms (const double &y_cd,
		       valarray<double> &ref_term,
		       valarray<double> &phi);

  //! return the liquid components reference terms (ref_term_) and gas components fugacities (fugacity_) for CO2 only.
  void referenceTermsCO2 (valarray<double> &ref_term,
			  valarray<double> &phi);

  //! return the liquid components reference terms (ref_term_) and gas components fugacities (fugacity_) for H2O only.
  void referenceTermsH2O (valarray<double> &ref_term,
			  valarray<double> &phi); 

  //! function for scaling the unknowns gas_est[i].
  void scalingUnknowns (valarray<double> &gas_est,
			valarray<double> &gas_scaling);
        
  //! function to calculate molar volumes for H2O and CO2.
  double vm(const double &y_cd,double &bm,
	    double &cm,double &dm,double &em,
	    double &ym,valarray<double> &b,
	    Matrix<double> &c, 
	    Matrix<double> &d, 
	    Matrix<double> &e);
    
  //! function to calculate molar volume for CO2 only.
  double vmCO2(double &bm,double &cm,
	       double &dm,double &em,
	       double &ym,valarray<double> &b, 
	       Matrix<double> &c, 
	       Matrix<double> &d, 
	       Matrix<double> &e);
    
  //! function to calculate molar volume for H2O only.
  double vmH2O(double &bm,double &cm,
	       double &dm,double &em,
	       double &ym,valarray<double> &b, 
	       Matrix<double> &c, 
	       Matrix<double> &d, 
	       Matrix<double> &e);

private:
    
  double bm_, cm_, dm_, em_;
  double p_, T_; 
  double vr_;
  double xwo_;
  //! wtfrex = weight fraction of gas phase
  double wtfrex_; 
  static const double R_; 
  static const double p0_h2o_; 
  static const double p0_co2_; 
  static const double whc_ ;// 
  Matrix<double> wo_, wij_;
  NumericalRecipes nr_;

  
  valarray<double> fol_, to_, par_;
  valarray<double>  ghiorso_par_, h2o_co2_;
  valarray<double> molecular_weight_; 
  valarray<double> gas_fraction_;
  valarray<double> ctrm_, xtrm_, rtrm_;
  valarray<double> scale_factor_;
  valarray<double> y_;
    
};



//! class required to apply secantMethod
class SolwcadFunction1 : public virtual Function {
    
public:
    
  SolwcadFunction1(Solwcad &solwcad,double &p1, double &p2,
		   double &p3, double &p4, double &p5)
  {
#ifdef DEBUG
    cerr << "in solwcadFunction1"<<endl;
#endif
    solwcad_ = &solwcad;
    p1_ = p1;
    p2_ = p2;
    p3_ = p3;
    p4_ = p4;
    p5_ = p5;
  }
  ~SolwcadFunction1() {;}
    
  double function(double &x)
  {
    double solution = solwcad_->funz(x,p1_,p2_,p3_,p4_,p5_);
    return solution;
  }
    
  double p5_;
protected:
    
private:
    
  double p1_;
  double p2_;
  double p3_;
  double p4_;
    
  Solwcad *solwcad_;
    
};
  
  
//! class required to apply secantMethod
class SolwcadFunction2 : public virtual Function {
   
public:

  SolwcadFunction2(Solwcad &solwcad,
		   double &p1, double &p2, double &p3)
  {
#ifdef DEBUG
    cerr << "dentro solwcadFunction2"<<endl;
#endif
    solwcad_ = &solwcad;
    p1_ = p1;
    p2_ = p2;
    p3_ = p3;
  }
  ~SolwcadFunction2() {;}
   
  double function(double &x)
  {
    double solution =  solwcad_->gasCalc(x,p1_, p2_, p3_);
    return solution;
  }
   
protected:
   
private:
   
  double p1_;
  double p2_;
  double p3_;
   
  Solwcad *solwcad_;
   
};

//! class required to apply newt
class SolwcadNewt : public virtual Function {
   
public:

  SolwcadNewt(Solwcad &solwcad)
  {
#ifdef DEBUG
    cerr << "dentro solwcadnewt"<<endl;
#endif
    solwcad_ = &solwcad;
  }
  ~SolwcadNewt() {;}
   
  void function2(valarray<double> &x, valarray<double> &y)
  {
    solwcad_->vecfunc(x,y);
  }
   
protected:
   
private:
   
  Solwcad *solwcad_;
   
};


#endif

