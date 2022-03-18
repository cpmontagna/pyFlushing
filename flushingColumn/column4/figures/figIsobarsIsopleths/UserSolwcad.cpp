// User file

#include <valarray>
#include <iostream>
#include <fstream>
#include <vector>
#include "Solwcad.hpp"
#include "util.h"
#include "Function.hpp"

using namespace std;


int main () {
  Solwcad solwcad;
  valarray<double> composition(10);
  valarray<double> gas_fraction(2);

  ofstream out_file;
  out_file.open("solwcadResultsIsopleths.out");

  double temperature = 1473;
  double initial_pressure = 1.e5;
  double final_pressure = 4.e8;
  double delta_pressure = 1.e7;

  double initial_CO2 = 0.02;
  double final_CO2 = 0.0001;
  double delta_CO2 = 1.e-4;
  double initial_water = 0.08;
  double final_water = 0.005;
  double delta_water = 1.e-3;
  
  cout << "Stromboli Metrich 2010"<<endl;
  composition[0]= 0.4741; //SiO2
  composition[1]= 0.0087; //TiO2
  composition[2]= 0.1496; //Al2O3
  composition[3]= 0.0164; //Fe2O3
  composition[4]= 0.0532; //FeO
  composition[5]= 0.0009; //MnO
  composition[6]= 0.0394; //MgO
  composition[7]= 0.1556; //CaO
  composition[8]= 0.022; //Na2O
  composition[9]= 0.016; //K2O
  //-------------------

  int vector_size = int((final_pressure-initial_pressure)/delta_pressure)+1;
  int co2_size = int((final_CO2 - initial_CO2)/delta_CO2)+1;
  int water_size = int((final_water - initial_water)/delta_water)+1;

  double pressure = 0.;
  double dissolved_h2o_onliq = 0.;
  double dissolved_co2_onliq = 0.;
  double exsolved_h2o_ongas = 0.;
  double exsolved_co2_ongas = 0.;
  double gas_ontotal = 0.;
  double carbon_dioxide = 0.;
  double water = 0.;

  for(int i=0; i<vector_size;i++)
    {
      pressure = final_pressure - (delta_pressure*(i));
      for(int j = 0; j < co2_size; j++)
	{
	  carbon_dioxide = final_CO2 - (delta_CO2 * j);
	  for(int k = 0; k < water_size; k++)
	    {
	      water = final_water - (delta_water * k);
	      solwcad.exsolution(pressure,temperature,water,
				 carbon_dioxide,composition);
	      dissolved_h2o_onliq = solwcad.wf_h2o_l_onliq();
	      dissolved_co2_onliq = solwcad.wf_co2_l_onliq();
	      exsolved_h2o_ongas = solwcad.wf_h2o_g_ongas();
	      exsolved_co2_ongas = solwcad.wf_co2_g_ongas();
	      gas_ontotal = solwcad.wf_gas_ontotal();

	      out_file << pressure << " "<< temperature <<" "<<
		water << " " << carbon_dioxide << " " <<
		dissolved_h2o_onliq<<" "<<dissolved_co2_onliq <<
		" "<<exsolved_h2o_ongas <<" "<<
		exsolved_co2_ongas <<" "<<gas_ontotal <<endl;
	    }
	}
    }
  out_file.close();
}
