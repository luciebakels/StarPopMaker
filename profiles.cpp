#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <random>
#include "proto.h"
#include "allvars.h"
#include "mpi.h"

double hernquistFactor(){
	return 1.0/(*All).CNFW * sqrt(2.0 * (std::log(1.0 + (*All).CNFW) - (*All).CNFW / (1.0 + (*All).CNFW)));
}
double totalHernquistMass(){
	return (*All).TotalMass/pow((1.0-hernquistFactor()), 2);
}
double virialRadius(double mass){ //in kpc
	return pow((3.0/4.0/M_PI*mass/(*All).DeltaVir/constants::rhocrit_Ms_kpci3), 1.0/3.0);
}
double hernquistScalingRadius(){
	return virialRadius((*All).TotalMass)*hernquistFactor()/(1.0-hernquistFactor());
}
double NFWScalingRadius(){
	return virialRadius((*All).TotalMass)/(*All).CNFW;
}
void GasProfile::setProfiles(){
	switch((*All).Profile){
		case PROFILE_HERNQUIST:
			HernquistProfile();
			break;
		case PROFILE_NFW:
			NFWProfile();
			break;
		case PROFILE_CONSTANT:
			ConstantProfile();
			break;
		case PROFILE_ISOTHERMAL:
			IsothermalProfile();
			break;
		default:
			std::cout<< "Assuming constant profile\n";
			IsothermalProfile();
	}
}
void GasProfile::logRadius(double maxRad){
	double logRad{std::log10(maxRad / (*All).LogBins)}; //Gives the radial stepsize in logspace to a max of maxRad
	double logStep{(std::log10(maxRad) - logRad)/(*All).LogBins};
	for(int i = 0; i < (*All).LogBins; i++){
		m_r[i] = pow(10.0, logRad + i*logStep);
	}	
}

// void GasProfile::HernquistProfile(){
// 	using namespace constants;
// 	logRadius(10.0 * virialRadius((*All).TotalMass));

// 	double r_hern{hernquistScalingRadius()}; //kpc
	
// 	double gam{(*All).GasCoeff};
// 	double gamfac{( 2.0 - gam ) * ( 3.0 - 2.0 * gam ) * ( 4.0 - 3.0 * gam ) / 2.0 / pow( kpc_to_cm*r_hern * 
// 		( gam - 1.0 ), 3)};
// 	double numberdensityfac{0.16/0.84*totalHernquistMass()*Msun_g/4.0/M_PI/m_particleMass*gamfac};
// 	double temperaturefac{(gam - 1.0)/gam * 2.0/3.0 * G_cm3_gi_si2 * totalHernquistMass()*Msun_g * 
// 		hydrogen_g / kB_erg_Ki / (r_hern*kpc_to_cm) };
	
// 	for(int i = 0; i < (*All).LogBins; i++){
// 		double temp{ 1.0 / (1.0 + m_r[i] / r_hern)};
// 		m_n[i] = numberdensityfac*pow(temp , 1.0 / (gam - 1.0));
// 		m_T[i] = temperaturefac*temp;
// 	}
// }
void GasProfile::HernquistProfile(){
	using namespace constants;

	double r_vir{virialRadius((*All).TotalMass) * kpc_to_cm}; //cm
	double r_hern{hernquistScalingRadius() * kpc_to_cm}; //cm

	logRadius(10.0 * virialRadius((*All).TotalMass));

	double gam{-0.1637*pow(std::log((*All).CNFW), -0.6614) + 1.312};
	double eta_0{0.3435*pow((*All).CNFW, 0.9020) + 1.075};
	
	//Fit made using python, only correct for Aaron's cNFW values
	double divide_factor_norm{1.4280 * pow((*All).CNFW, 0.04213) - 1.4001};

	double m_c{r_vir*r_vir/r_hern/r_hern /(1.0 + r_vir/r_hern)/(1.0 + r_vir/r_hern)};
	double rho_0{(*All).TotalMass*Msun_g/m_particleMass / (4.0 * M_PI * pow(r_hern, 3) * m_c)};
	double T_0{eta_0 * (*All).TotalMass*Msun_g * hydrogen_g * G_cm3_gi_si2 / ( 3.0 * r_vir * kB_erg_Ki)};
	double fac{3.0/eta_0 * (gam - 1.0)/gam * r_vir / r_hern / m_c};

	for(int i = 0; i < (*All).LogBins; i++){
		double y{1.0 - fac * (1.0 - 1.0 / (1.0 + m_r[i] * kpc_to_cm/r_hern))};
		m_n[i] = rho_0 * pow(y, 1.0/(gam - 1.0))/divide_factor_norm*0.16/0.84;
		m_T[i] = T_0 * y;
	}
}

void GasProfile::NFWProfile(){ //Komatsu and Seljak 2001
	using namespace constants;

	double r_vir{virialRadius((*All).TotalMass) * kpc_to_cm}; //cm
	double r_s{NFWScalingRadius() * kpc_to_cm}; //cm

	logRadius(10.0 * virialRadius((*All).TotalMass));

	double gam{1.15 + 0.01 * ((*All).CNFW - 6.5)};
	double eta_0{0.00676 * pow((*All).CNFW - 6.5, 2) + 0.206 * ((*All).CNFW - 6.5) + 2.48};

	//Fit made using python, only correct for Aaron's cNFW values
	//double divide_factor_norm{6.79 * pow(0.5, 0.6517*std::log10((*All).TotalMass) - 3.718) + 0.2968};
	double divide_factor_norm{5.5224*exp(0.2464*(*All).CNFW - 4.7292) + 0.3003};

	double m_c{std::log(1.0 + (*All).CNFW) - (*All).CNFW/(1.0 + (*All).CNFW)};
	double rho_0{(*All).TotalMass*Msun_g/m_particleMass / (4.0 * M_PI * pow(r_s, 3) * m_c )};
	double T_0{eta_0 * (*All).TotalMass*Msun_g * hydrogen_g * G_cm3_gi_si2 / ( 3.0 * r_vir * kB_erg_Ki)};
	
	double fac{3.0/eta_0 * (gam - 1.0)/gam * (*All).CNFW / m_c};

	for(int i = 0; i < (*All).LogBins; i++){
		double y{1.0 - fac * (1.0 - std::log(1.0 + m_r[i] * kpc_to_cm/r_s)/(m_r[i] * kpc_to_cm/r_s))};
		m_n[i] = rho_0 * pow(y, 1.0/(gam - 1.0))/divide_factor_norm*0.16/0.84;
		m_T[i] = T_0 * y;
	}
}

void GasProfile::ConstantProfile(){
	using namespace constants;
	logRadius((*All).Radius);

	for(int i = 0; i < (*All).LogBins; i++){
		m_n[i] = 0.16/0.84*((*All).TotalMass * Msun_g)/(4.0/3.0*M_PI*pow((*All).Radius*kpc_to_cm, 3));
		m_T[i] = 4.0/3.0 * M_PI * G_cm3_gi_si2 * m_particleMass * m_n[i] * pow(m_r[i]*kpc_to_cm, 2) / 5.0 / kB_erg_Ki; //mo page 374
	}
}

void GasProfile::IsothermalProfile(){ //Binney page 304
	using namespace constants;
	double sigma2{G_kpc_km2_Msi_si2*0.16/0.84*(*All).TotalMass/(*All).Radius/2.0 * 1.0e10}; //cm2/s2
	logRadius((*All).Radius);

	double tempT{sigma2 * m_particleMass / kB_erg_Ki};
	for(int i = 0; i < (*All).LogBins; i++){
		m_n[i] = sigma2 / 2.0 / M_PI / G_cm3_gi_si2 / (m_r[i]*kpc_to_cm) / (m_r[i]*kpc_to_cm) / m_particleMass;
		m_T[i] = tempT;
	}
}

void GasProfile::setCooling(){
	/* find m_TLam closest to m_T, then get m_Lambda in m_cooling */
	using namespace constants;

	for(int j = 0; j < m_T.size(); j++){
		double val1{m_T[j] - m_TLam[0]}, val2{0};
		m_cooling[j] = m_Lambda[0]*m_n[j]*m_n[j]*nH_fraction*nH_fraction; 
		for(int i = 1; i < m_TLam.size(); i++){
			val2 = (m_T[j] - m_TLam[i]);
			if(val2*val2 > val1*val1)
				break;
			else if(val2*val2 < val1*val1){
				m_cooling[j] = m_Lambda[i]*m_n[j]*m_n[j]*nH_fraction*nH_fraction;	 	// Lambda = Cooling/nH^2
				val1 = val2;
			}
			else{
				std::cout << "Error in setCooling: T = " << m_T[j] << "\n";
				throw std::exception();
			}
		}
	}
	/* integrate over logshells for m_coolIntegrated. All profiles are spherically symmetric. */
	double r1;
	double r2{0};		//Middle of sphere, starting point
	for(int j = 0; j < m_r.size(); j++){
		if( j == m_r.size() - 1)
			r1 = m_r[j]*kpc_to_cm; //in cm
		else
			r1 = pow(10.0, std::log10(m_r[j]) + 0.5*std::log10(m_r[j+1]/m_r[j]))*kpc_to_cm; //Logbin + halfway next bin (in logscale)
		double shellVolume{4.0/3.0*M_PI*(r1*r1*r1 - r2*r2*r2)};
		m_coolIntegrated[j] = shellVolume*m_cooling[j];
		r2 = r1;
	}
	/* Cooling time 3kT/(2*Lambda(T)*n) */
	for(int i = 0; i < m_r.size(); i++){ //Mo page 386
		m_coolingTime[i] = 3.0 * m_n[i] * kB_erg_Ki * m_T[i] / 2.0 / m_cooling[i];
	}
	/* Bremsstrahlung cooling time. Peacock page 570 */
	for(int i = 0; i < m_r.size(); i++){
		m_coolingTimeff[i] = 1.8e24 / (1.0/sqrt(m_T[i]/1.0e8) + 0.5*0.5/pow(m_T[i]/1.0e8, (3.0/2.0))) /
			(m_n[i]*m_particleMass/Msun_g*pow(1.0e3*kpc_to_cm, 3));
	}
	// /* Hernquist cooling, test */
	// for(int i = 0; i < m_r.size(); i++){
	// 	m_coolingtest[i] = 8.051e-28*pow((*All).TotalMass, 1.0/3.0)	*
	// 		pow(1.0 + 3.121e2*pow((*All).TotalMass, -1.0/3.0)*m_r[i], -17.167);
	// }
	/* Bremsstrahlung cooling */ //Mo page 368
	for(int i = 0; i < m_r.size(); i++){
		m_coolingff[i] = 1.4e-23*sqrt(m_T[i]/1e8)*m_n[i]*m_n[i];
	}
	/* Cooling time bremsstrahlung test */
	for(int i = 0; i < m_r.size(); i++){
		m_coolingtest[i] = 3.0 * m_n[i] * kB_erg_Ki * m_T[i] / 2.0 / m_coolingff[i];
	}
}

void GasProfile::writeProperties(){
	std::cout << "Writing gas properties...\n";
	std::vector<std::vector<double>> allproperties;
	std::ofstream myfile;

	allproperties.resize(9);
	allproperties[0] = m_r;
	allproperties[1] = m_n;
	allproperties[2] = m_T;
	allproperties[3] = m_cooling;
	allproperties[4] = m_coolIntegrated;
	allproperties[5] = m_coolingTime;
	allproperties[6] = m_coolingTimeff;
	allproperties[7] = m_coolingtest;
	allproperties[8] = m_coolingff;

	std::string filename1;
	switch((*All).Profile){
		case PROFILE_HERNQUIST:
			filename1 = "GasProfilePropertiesH";
			break;
		case PROFILE_NFW:
			filename1 = "GasProfilePropertiesN";
			break;
		case PROFILE_CONSTANT:
			filename1 = "GasProfilePropertiesC";
			break;
		case PROFILE_ISOTHERMAL:
			filename1 = "GasProfilePropertiesI";
			break;
		default:
			filename1 = "GasProfileProperties";
	}
	double massexp{std::log10((*All).TotalMass)};
	myfile.open(filename1 + std::to_string(static_cast<int>(massexp*10)) + "CNFW" + 
		std::to_string(static_cast<int>((*All).CNFW*100)) + ".txt", 
		std::ofstream::out | std::ofstream::trunc);
	myfile << "r\tn\tT\tC\tCintegrated\tCoolingTime\tCoolingTimeff\tCoolingtest\tCoolingff\n";
	myfile.close();
	writeVector(allproperties, filename1 + std::to_string(static_cast<int>(massexp*10)) + "CNFW" + 
		std::to_string(static_cast<int>((*All).CNFW*100)) + ".txt", 'a');
	std::cout << "Done writing gas properties.\n";
}