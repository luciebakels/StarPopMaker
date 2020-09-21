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

int Particle::s_NumGas = 0;

int Gas::s_NumIonized = 0;
Atoms Gas::s_Eth;
Atoms Gas::s_mass; //u


GASTYPE getGasType(){
	GasProperties gasproperties;
	double randval{randomNumber(0.0, 0.955)};
	if(randval <= 0.005){
		return GASTYPE_MOLECULARCLOUD;
	}
	else if(randval > 0.005 && randval <= 0.055){
		return GASTYPE_DIFFUSECLOUD;
	}
	else if(randval > 0.055 && randval <= 0.455){
		return GASTYPE_INTERCLOUDMEDIUM;
	}
	else{
		return GASTYPE_CORONALGAS;
	}
}
GasProperties getGasProperties(GASTYPE gastype){
	GasProperties gasproperties;
	switch(gastype){
		case GASTYPE_MOLECULARCLOUD:
			gasproperties.volume = 0.005;
			gasproperties.massfraction = 0.4;
			gasproperties.density = 1000.0*constants::hydrogen_g;
			gasproperties.temperature = 20.0;
			gasproperties.c_sound = 0.6e5; //cm/s
			break;
		case GASTYPE_DIFFUSECLOUD:
			gasproperties.volume = 0.05;
			gasproperties.massfraction = 0.4;
			gasproperties.density = 100.0*constants::hydrogen_g;
			gasproperties.temperature = 80.0;
			gasproperties.c_sound = 0.9e5; //cm/s
			break;
		case GASTYPE_INTERCLOUDMEDIUM:
			gasproperties.volume = 0.4;
			gasproperties.massfraction = 0.2;
			gasproperties.density = 1.0*constants::hydrogen_g;
			gasproperties.temperature = 8000.0;
			gasproperties.c_sound = 9.0e5; //cm/s
			break;
		case GASTYPE_CORONALGAS:
			gasproperties.volume = 0.545;
			gasproperties.massfraction = 0.001;
			gasproperties.density = 0.001*constants::hydrogen_g;
			gasproperties.temperature = 10.0e6;
			gasproperties.c_sound = 100.0e5; //cm/s
			break;
		default:
			gasproperties.volume = 0.5;
			gasproperties.massfraction = 0.001;
			gasproperties.density = 0.001*constants::hydrogen_g;
			gasproperties.temperature = 10.0e6;
			gasproperties.c_sound = 100.0e5; //cm/s
			break;
	}			
	return gasproperties;
}


void Gas::measureDistanceToStars(std::vector<Star> star){
	m_distanceStar.resize(Particle::s_NumStar);
	for(int i = 0; i < Particle::s_NumStar; i++){
		m_distanceStar[star[i].ID()] = distanceBetween(star[i]);
	}
}
void Gas::changeNum(double H, double He, double HeI, double E){
	m_speciesNum.H -= H;
	m_speciesNum.HI += H;
	m_speciesNum.He -= He;
	m_speciesNum.HeI += (He - HeI);
	m_speciesNum.HeII += HeI;
	m_speciesNum.e = m_speciesNum.HeII*2.0 + m_speciesNum.HI;

	m_speciesFrac.H = m_speciesNum.H/m_speciesNum.X;
	m_speciesFrac.HI = m_speciesNum.HI/m_speciesNum.X;
	m_speciesFrac.He = m_speciesNum.He/m_speciesNum.X;
	m_speciesFrac.HeI = m_speciesNum.HeI/m_speciesNum.X;
	m_speciesFrac.HeII = m_speciesNum.HeII/m_speciesNum.X;

	if(H > 0)
		addEnergy( (E - Gas::s_Eth.H)*H );
	else
		addEnergy( Gas::s_Eth.H*H );
	if(He > 0)
		addEnergy( (E - Gas::s_Eth.He)*He );
	else
		addEnergy( Gas::s_Eth.He*He );
	if(HeI > 0)
		addEnergy( (E - Gas::s_Eth.HeI)*HeI );
	else
		addEnergy( Gas::s_Eth.HeI*HeI );

	m_numberDensity.H = m_speciesNum.H*m_density/m_mass;
	m_numberDensity.HI = m_speciesNum.HI*m_density/m_mass;
	m_numberDensity.He = m_speciesNum.He*m_density/m_mass;
	m_numberDensity.HeI = m_speciesNum.HeI*m_density/m_mass;
	m_numberDensity.HeII = m_speciesNum.HeII*m_density/m_mass;
	m_numberDensity.X = m_speciesNum.X*m_density/m_mass;
	m_numberDensity.e = m_speciesNum.e*m_density/m_mass;	
}
void Gas::assignMassDensity(double mass, double density) {
	m_mass = mass;
	m_density = density;
	m_speciesNum.H = (m_speciesFrac.H*m_mass/s_mass.H*100.0);
	m_speciesNum.HI = (m_speciesFrac.HI*m_mass/s_mass.H);
	m_speciesNum.He = (m_speciesFrac.He*m_mass/s_mass.He*100.0);
	m_speciesNum.HeI = (m_speciesFrac.HeI*m_mass/s_mass.He);
	m_speciesNum.HeII = (m_speciesFrac.HeII*m_mass/s_mass.He);
	m_speciesNum.X = m_speciesNum.H + m_speciesNum.HI + m_speciesNum.He + m_speciesNum.HeI + m_speciesNum.HeII;
	m_speciesNum.e = m_speciesNum.HeII*2.0 + m_speciesNum.HI;

	m_numberDensity.H = m_speciesNum.H*m_density/m_mass;
	m_numberDensity.HI = m_speciesNum.HI*m_density/m_mass;
	m_numberDensity.He = m_speciesNum.He*m_density/m_mass;
	m_numberDensity.HeI = m_speciesNum.HeI*m_density/m_mass;
	m_numberDensity.HeII = m_speciesNum.HeII*m_density/m_mass;
	m_numberDensity.X = m_speciesNum.X*m_density/m_mass;
	m_numberDensity.e = m_speciesNum.e*m_density/m_mass;
}
void Gas::recombinations(Atoms alphaR, double timestep) {
	double H{0}, He{0}, HeI{0};
	H = 1000000*timestep * alphaR.H * m_numberDensity.e;
	HeI = timestep * alphaR.HeI * m_numberDensity.HeII * m_numberDensity.e;
	He = timestep * alphaR.He * (m_numberDensity.HeI + HeI) * m_numberDensity.e;
	//std::cout << H << ", " << He << ", " << HeI << "\n";
	//std::cout << m_numberDensity.HI << ", " << m_numberDensity.e << "\n";
	changeNum(-H, -He, -HeI);
}