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

int Particle::s_NumStar = 0;
int Star::s_ActiveStars = 0;

void Star::alphabeta(){ //Power(2013). Only stars between 0.5 and 100
	if ((m_mass>=0.5) && (m_mass<3.0))
		m_alphabeta = std::array<double, 2>{0.8, 4.5};
	else if ((m_mass>=3.0) && (m_mass<10.0))
		m_alphabeta = std::array<double, 2>{1.8, 3.6};
	else if ((m_mass>=10.0) && (m_mass<20.0))
		m_alphabeta = std::array<double, 2>{5.8, 3.1};
	else if ((m_mass>=20.0) && (m_mass<50.0))
		m_alphabeta = std::array<double, 2>{25.8, 2.6};
	else if ((m_mass>=50.0) && (m_mass<=100.0))
		m_alphabeta = std::array<double, 2>{100.0, 2.3};
	else{
		std::cout<< "Error: incorrect mass in alphabeta: " << m_mass << "\n";
		throw std::exception();
	}
}
double Star::lifetime(double mass){//Marigo2001 table 1
	double a, b;
	if (mass < 0.5)
		return (*All).TotalTime;
	else 
		return 1.0e10*(0.29385*pow(mass, -2.50746) + 1.0/(100.0*mass) + 
			0.502346*(pow(exp(mass-0.89872), -5.67346)) + 1.96e-4);
}
// double Star::lifetime(double mass){ //Self made fit from Marigo2001 table 1
// 	double a, b, c;
// 	if (mass <0.1)
// 		return (*All).TotalTime;
// 	else if ((mass>=0.1) && (mass<3.0)){
// 		a = 0.571;
// 		b = 3.8784;
// 		c = 6.7135e-3;
// 	}
// 	else if ((mass>=3.0) && (mass<10.0)){
// 		a = 0.279;
// 		b = 2.7695;
// 		c = 1.4592e-3;
// 	}
// 	else if ((mass>=10.0) && (mass<20.0)){
// 		a = 0.0286;
// 		b = 1.1700;
// 		c = 0.0;
// 	}
// 	else if ((mass>=20.0) && (mass<50.0)){
// 		a = 0.00957;
// 		b = 0.8045;
// 		c = 0.0;
// 	}
// 	else if ((mass>=50.0) && (mass<=100.0)){
// 		a = 0.00280;
// 		b = 0.4902;
// 		c = 0.0;
// 	}
// 	else{
// 		std::cout<< "Error: incorrect mass in lifetime: " << m_mass << "\n";
// 		throw std::exception();
// 	}
// 	return 1.0e10*(a*pow(mass, -1) +c);	
// }
void Star::assignMassRemnant(){ //Heger, Woosley, Fryer, Langer 2003 Fig. 2, solar metallicity
	if(m_mass < 2.7){
		m_massRemnant = 0.096*m_mass + 0.429;
	}
	else if (m_mass < 8){
		m_massRemnant = 0.137*m_mass + 0.318;
	}
	else if(m_mass > 80){
		m_massRemnant = m_mass/10.0;
	}
	else if(m_mass > 46){
		m_massRemnant = 8.0;
	}
	else if(m_mass > 30){
		m_massRemnant = 10.0 - 2.0/16.0*(m_mass - 30.0);
	}
	else if(m_mass > 15){
		m_massRemnant = 0.583*(m_mass - 15.0) + 1.25;
	}
	else if(m_mass > 10){
		m_massRemnant = 1.25;
	}
	else{
		m_massRemnant = 0.10*(m_mass - 8.0) + 1.05;
	}

}
// void Star::assignMassRemnant(){ //Heger 2003 Fig. 5, zero metallicity
// 	if(m_mass < 8.0){
// 		m_massRemnant = 0.0;
// 	}
// 	else if(m_mass > 40.0){
// 		m_massRemnant = m_mass;
// 	}
// 	else if(m_mass > 38.0){
// 		m_massRemnant = m_mass - 23.0;
// 	}
// 	else if(m_mass > 35.0){
// 		m_massRemnant = 0.5*m_mass - 4.0;
// 	}
// 	else if(m_mass > 32.0){
// 		m_massRemnant = m_mass - 21.5;
// 	}
// 	else if(m_mass > 31.0){
// 		m_massRemnant = 0.5*m_mass - 5.5;
// 	}
// 	else if(m_mass > 30.0){
// 		m_massRemnant = m_mass - 21.0;
// 	}
// 	else if(m_mass > 29.0){
// 		m_massRemnant = 8.5*m_mass - 245.0;
// 	}
// 	else if(m_mass > 27.0){
// 		m_massRemnant = 1.5;
// 	}
// 	else if(m_mass > 26.0){
// 		m_massRemnant = -5.5*m_mass + 150.0;
// 	}
// 	else if(m_mass > 24.0){
// 		m_massRemnant = 0.25*m_mass + 0.5;
// 	}
// 	else if(m_mass > 23.0){
// 		m_massRemnant = 5.0*m_mass - 113.5;
// 	}
// 	else if(m_mass > 19.0){
// 		m_massRemnant = 1.5;
// 	}
// 	else if(m_mass > 18.0){
// 		m_massRemnant = -2.5*m_mass + 49.0;
// 	}
// 	else if(m_mass > 17.0){
// 		m_massRemnant = m_mass - 14.0;
// 	}
// 	else if(m_mass > 16.0){
// 		m_massRemnant = 0.6667*m_mass - 9.1667;
// 	}
// 	else if(m_mass >= 8.0){
// 		m_massRemnant = 1.5;
// 	}
// }
void Star::assignRadius(){ //cm
	using namespace constants;
	m_R = pow(m_mass, 0.8)*Rsun_cm;
}
void Star::assignTms(){ //In yr
	using namespace constants;
	m_tms = lifetime(m_mass);
}
void Star::assignLuminosity(){ //In erg s
	using namespace constants;
	if(m_mass < 0.5)
		m_L = 0.0;
	else if(m_mass > 50.0)
		m_L = Lsun_erg_s*(pow(10, (-0.04172 + 4.4954*std::log10(m_mass) - 0.6041*pow((std::log10(m_mass)), 2))));
	else
		m_L = m_alphabeta[0]*pow(m_mass, m_alphabeta[1])*Lsun_erg_s;
}
void Star::assignTemperature(){ //In units of Kelvin
	using namespace constants;
	m_T = pow((m_L/(4.0*M_PI*m_R*m_R*sigmaSB_erg_cm2_si_K4)), 0.25);
}
double Star::blackBody(double E){ //eV cm^-2 * (h_eV_s*h_eV_s*c_cm_si*c_cm_si)
	using namespace constants;
	if(m_T == 0){
		std::cout<< "Error: blackBody(), T = 0\n";
		throw std::exception();
	}
	double below{exp(E/kB_eV_Ki/m_T) - 1.0};
	if(below <= 0)
		return 2.0*E*E/(1.0/kB_eV_Ki/m_T);
	else
		return 2.0*E*E*E/(below);
}
double Star::blackBodyPhoton(double E){ //cm^-2 * (h_eV_s*h_eV_s*c_cm_si*c_cm_si)
	using namespace constants;
	if(m_T == 0){
		std::cout<< "Error: blackBodyPhoton(), T = 0\n";
		throw std::exception();
	}
	double below{exp(E/kB_eV_Ki/m_T) - 1.0};
	if(below <= 0)
		return 2.0*E/(1.0/kB_eV_Ki/m_T);
	else
		return 2.0*E*E/(below);
}
double Star::integrate(double lowlim, double uplim, double (Star::*function)(double)){
	double currentup{uplim}, currentlow{lowlim}, fup, flow, fmid, max;
	if(type() == TYPE_HMXB)
		max = 1e35;
	else
		max = findmax(0, std::max((*All).MaximumEnergy, uplim), function);
	if(max > 1e50){
		std::cout << "Something wrong in findmax()\n";
		throw std::exception();
	}
	double errorfind{max*(uplim - lowlim)};
	double result{0};
	while(currentlow < uplim){
		fup = (*this.*function)(currentup);
		flow = (*this.*function)(currentlow);
		fmid = (*this.*function)(currentup - (currentup - currentlow)/2.0);
		if( (fmid > 0.99999*(std::min(flow, fup) + (std::max(flow, fup) - std::min(flow, fup))/2.0)) &&
			(fmid < 1.00001*(std::min(flow, fup) + (std::max(flow, fup) - std::min(flow, fup))/2.0))){
			double temp{fmid*(currentup - currentlow)};
			if(temp > 10*errorfind){
				std::cout << "Integration failed! Stepsize: " << temp << ", maximum: " << errorfind << 
				" Max curve: " << max << " Type: " << type() << " Temperature HMXB: " << m_THMXB << " Luminosity HMXB: " << m_LHMXB <<"\n";
				temp = 0.0;
				throw std::exception();
			}
			result += fmid*(currentup - currentlow);
			currentlow = currentup;
			currentup = uplim;
		}
		else if( (std::max(flow, fup) < 0.01*max) && ((currentup - currentlow) < 0.1) ){
			result += fmid*(currentup - currentlow);
			currentlow = currentup;
			currentup = uplim;
		}
		else{
			currentup = currentup - (currentup - currentlow)/2.0;
		}
	}
	return result;
}
double Star::findmax(double x1, double x2, double (Star::*function)(double)){
	double xdelta{(x2 - x1)/4.0}, max{(*this.*function)(x1)}, test{0};
	int imax{0};
	if(xdelta < 1.0e-15){
		return x1;
	}
	for(int i = 1; i < 6; i++){
		test = (*this.*function)(x1 + i*xdelta);
		if(test > max){
			imax = i;
			max = test;
		}
	}
	if(imax == 0)
		return x1;
	else if(imax == 5)
		return x2;
	else{
		return findmax(x1+(imax-1)*xdelta, x1 + (imax+1)*xdelta, function);
	}
}
void Star::assignIonPhotonsLum(double lowlimE, double uplimE){
	using namespace constants;
	m_Lion = 0;
	m_Nion = 0;
	if(m_mass > 0.5 && uplimE > 13.6){
		m_Nion = integrate(std::max(13.6, lowlimE), uplimE, &Star::blackBodyPhoton);
		m_Lion = integrate(std::max(13.6, lowlimE), uplimE, &Star::blackBody);
		
		m_Lion = m_Lion/(h_eV_s*h_eV_s*c_cm_si*c_cm_si)*
			m_R*m_R*M_PI*4.0*M_PI/h_eV_s/erg_to_eV; //4*pi*r^2 * pi
		m_Nion = m_Nion/(h_eV_s*h_eV_s*c_cm_si*c_cm_si)*
			m_R*m_R*M_PI*4.0*M_PI/h_eV_s;
	}
}

void Star::accretionLuminosity(GasProperties gasproperties){ //Fender2013
	using namespace constants;
	double mBondi{4.0*M_PI*accretionEfficiency*(G_cm3_gi_si2*m_massRemnant*Msun_g)*
		(G_cm3_gi_si2*m_massRemnant*Msun_g)*constants::hydrogen_g / 
		gasproperties.density*pow((pow(m_Vel.distance(), 2) + gasproperties.c_sound*gasproperties.c_sound), -1.5)};
	m_Lacc = radiativeEfficiency * mBondi*c_cm_si*c_cm_si;
}

void Star::accretionLumNeutral(){ //Fender2013
	using namespace constants;
	double mBondi{4.0*M_PI*(G_cm3_gi_si2*m_massRemnant*Msun_g)*
		(G_cm3_gi_si2*m_massRemnant*Msun_g)*constants::hydrogen_g};
	m_Lacc =  mBondi*c_cm_si*c_cm_si;
}

double Star::HMXBLum(double E){
	using namespace constants;
	if(m_THMXB == 0){
		std::cout<< "Error: HMXBLum(), T = 0\n";
		throw std::exception();
	}
	double below{exp(E/kB_eV_Ki/m_THMXB) - 1.0};
	if(E>5000){
		if(below <= 0)
			return 2.0*E*E/(1.0/kB_eV_Ki/m_THMXB);// + 2.95e14*pow(E, -1.1);
		else{
			return 2.0*E*E*E/(below);// + 2.95e14*pow(E, -1.1);
		}
	}
	else{
		if(below <= 0)
			return 2.0*E*E/(1.0/kB_eV_Ki/m_THMXB);// + 2.148e6*pow(E, 1.1);
		else
			return 2.0*E*E*E/(below);// + 2.148e6*pow(E, 1.1);		
	}
}
double Star::HMXBPhotons(double E){
	using namespace constants;
	if(m_LHMXB < 1.0e37){
		if(E < 2.0e3)
			return 0.0;
		else
			return m_LHMXB/8.6826*pow(E, -1.8)*erg_to_eV;
	}
	else{
		if(m_THMXB == 0){
			std::cout<< "Error: HMXBPhotons(), T = 0\n";
			throw std::exception();
		}
		double below{exp(E/kB_eV_Ki/m_THMXB) - 1.0};
		if(E>5000){
			if(below <= 0)
				return 2.0*E/(1.0/kB_eV_Ki/m_THMXB);// + 2.95e14*pow(E, -2.1);
			else
				return 2.0*E*E/(below);// + 2.95e14*pow(E, -2.1);
		}
		else{
			if(below <= 0)
				return 2.0*E/(1.0/kB_eV_Ki/m_THMXB);// + 0.0859*pow(E, 2.1);
			else
				return 2.0*E*E/(below);// + 0.0859*pow(E, 2.1);			
		}

	}
}
void Star::assignEnergy(){
	m_E.resize((*All).EnergyBins);
	if((*All).EnergyLogOn){
		double x1{std::log10((*All).MinimumEnergy)}, x2{std::log10((*All).MaximumEnergy)};
		double logstep{(x2-x1)/(*All).EnergyBins};
		for(int i = 0; i < (*All).EnergyBins; i++)
			m_E[i] = pow(10.0, x1+i*logstep/(*All).EnergyBins);
	}
	else{
		for(int i = 0; i < (*All).EnergyBins; i++){
			m_E[i] = (*All).MinimumEnergy + i*((*All).MaximumEnergy - (*All).MinimumEnergy)/(*All).EnergyBins;
		}
	}
}
double weibull(double randval){
	double ma{1e38};
	double mi{1e33};
	double k{1.9};
	//1.0e30 + 1.e38/0.337*0.5*pow(std::log(1.0/(1.0 - randomNumber(0.0, 1.0))), 0.5263);
	return pow(10, (std::log10(ma/mi)*pow((-std::log(randval)),(1/k))/pow(((k-1)/k),(1/k)) + std::log10(mi)));
}
void Star::assignHMXBproperties(){
	using namespace constants;
	assignEnergy();
	m_NHMXB.resize((*All).EnergyBins);
	m_LHMXB = weibull(randomNumber(0.0, 1.0));
	m_LIonHMXB = 0.0;
	while(m_LHMXB > 1.26e38*m_massRemnant)
		m_LHMXB = weibull(randomNumber(0.0, 1.0));//1.26e38*m_massRemnant;////1.0e30 + 1.e38/0.337*0.5*pow(std::log(1.0/(1.0 - randomNumber(0.0, 1.0))), 0.5263);
	double energy{m_E[0]}, energynew, norm{h_eV_s};
	m_THMXB = pow(m_LHMXB/lum_cyg_X1, 0.25)*T_cyg_X1;
	double HMXBpowerlawNorm{HMXBLum(500)/0.5384};
	if(m_LHMXB > 1.0e37){
		norm = m_LHMXB/(integrate(0.1, 1e5, &Star::HMXBLum)*1.0004)*erg_to_eV;//+HMXBpowerlawNorm*(4.268))*erg_to_eV;
		m_LIonHMXB = norm/erg_to_eV*integrate(13.6, 1e3, &Star::HMXBLum);
		if(m_LIonHMXB > m_LHMXB){
			std::cout << "HMXB spectrum integration failed: assign average\n";
			m_LIonHMXB = 3.0e36;
		}
	}
	//std::cout << m_LHMXB << "\t" << m_LIonHMXB << "\t" << m_THMXB << "\n";
	for(int i = 0; i < m_E.size(); i++){
		if(i==(m_E.size()-1)){
			if(m_E[i] < 1000){
				m_NHMXB[i] = norm*integrate(energy, m_E[i], &Star::HMXBPhotons);
				if(m_NHMXB[i] > 1e50){
					std::cout << "HMXB photon integration failed: skip energy bin " << m_E[i] << " eV\n";
					m_NHMXB[i] = 0;
				}
			}
			else
				m_NHMXB[i] = norm*(integrate(energy, m_E[i], &Star::HMXBPhotons) + 
				HMXBpowerlawNorm/1.1*(pow(energy, -1.1) - pow(m_E[i], -1.1)));
		}
		else{
			energynew = m_E[i+1] - (m_E[i+1]-energy)/2.0;
			if(m_E[i] < 1000)
				m_NHMXB[i] = norm*integrate(energy, energynew, &Star::HMXBPhotons);
			else
				m_NHMXB[i] = norm*(integrate(energy, energynew, &Star::HMXBPhotons) + 
				HMXBpowerlawNorm/1.1*(pow(energy, -1.1) - pow(energynew, -1.1)));
			energy = energynew;
		}
	}
}
void Star::assignWDproperties(){
	//H.M. van Horn 1971 (Metsel cooling): 1.7e-3*M*(4e23/K0)*mu/mue^2*(Tc/1e7)**7/2
	// - K0 = 4e23 (???)
	// mu = 2 (Helium envelope)
	// mue = 1/(1+X), X=0.7 for solar metallicity
	// Tc = 1.27e7 K

	m_LWD0 = 0.0227 * m_massRemnant;
	m_LWD = m_LWD0;
	m_twd = (*All).TimeStep/2.0;

	//m_A = atomic mass
	//Mass regions from Catalin 2008
	if (m_massRemnant > 1.05)//64% Ne:36% core, Garcia-Burro & Hernanz 1997
		m_A = 17.44;
	else if (m_massRemnant > 0.4)//Xc=40%, Yo=57%, Salaris 1997
		m_A = 14.04;
	else //Helium
		m_A = 4;
	updateWDluminosity();
}

void Star::updateWDluminosity(){
	//Derived from Kepler and Bradley 1995: tcool = 8.8e6 * (A/12)^-1 M^5/7 (mu/2)^-2/7 L^(-5/7)
	// L(t) = L0(1 - 1.14e-7 (A/12) M^-5/7 (mu/2)^-2/7 L0^5/7 t)^7/5
	if(m_twd == -1){
		m_LWD = 0.0;
	}
	else{
		double temp;
		temp = 1.0 - 1.14e-7 * m_A/12 * pow(m_LWD0/m_massRemnant, 0.714285) * m_twd;
		//std::cout << "WD: " << temp << m_twd << m_LWD0 << m_massRemnant << m_A << "\n";
		if(temp <= 0){
			m_LWD = 0.0;
			m_twd = -1;
		}
		else
			m_LWD = m_LWD0 * pow(temp, 1.4);
	}
}

void Star::typechange(){
	if(m_type == TYPE_STAR){
		if((m_mass >= 8.0) && (m_massRemnant + m_massCompanion >= 0.5*m_mass) && (randomNumber(0.0, 1.0) <= (*All).SurvivalFraction) ){
			m_SN = true;
			if(m_massCompanion >= 3.0){
				m_type = TYPE_HMXB;
				m_thmxb = lifetime(m_massCompanion) - lifetime(m_mass);
				assignHMXBproperties();
			}
			else{
				m_thmxb = lifetime(m_massCompanion) - lifetime(m_mass);
				m_type = TYPE_LMXB;
			}
		}
		else if(m_mass >= 8.0 && m_mass <= 20.0){
			m_SN = true;
			m_type = TYPE_NS;
			//m_gastype = getGasType();
			accretionLumNeutral();
		}
		else if(m_mass > 20.0){
			m_SN = true;
			m_type = TYPE_BH;
			//m_gastype = getGasType();
			accretionLumNeutral();
		}
		else
			m_type = TYPE_WD;
			assignWDproperties();
	}
}
void Star::Active(bool active){
	if(m_active == true && active == false){
		s_ActiveStars--;
		m_active = active;
		typechange();
	}
}

void Star::reduceLifeWith(double minus){ 
	switch(m_type){
		case TYPE_STAR:
			m_tms -= minus;
			if(m_tms <= 0)
				Active(false);
			break;
		case TYPE_HMXB:
			m_thmxb -= minus;
			if(m_thmxb <= 0){
				if(m_mass > 20.0){
					m_type = TYPE_BH;
					accretionLumNeutral();
				}
				else{
					m_type = TYPE_NS;
					accretionLumNeutral();
				}
			}
			break;
		case TYPE_LMXB:
			m_thmxb -= minus;
			if(m_thmxb <= 0){
				if(m_mass > 20.0){
					m_type = TYPE_BH;
				}
				else{
					m_type = TYPE_NS;
				}
			}
			break;
		case TYPE_WD:
			m_twd += minus;
			updateWDluminosity();
			break;
		default:
			break;
	}
}
