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

void eraseStuff(std::string &mystring, std::string erasethis){
	//std::cout << mystring << "\t" << erasethis.size() << "\n";
	mystring.erase(0, erasethis.size());
	mystring.erase(std::remove_if(mystring.begin(), mystring.end(), [](char x){return std::isspace(x);}), mystring.end());
}

void setVolume(){
	double totmass;
	switch((*All).SFRType){
		case SFR_STARBURST || SFR_EXPONENTIAL || SFR_DISK:
			totmass = (*All).TotalMass;
			break;
		case SFR_CONSTANT:
			totmass = (*All).TotalTime*(*All).Sfr;
			break;
		default:
			if(!(*All).TotalMass){
				std::cout << "Give me total mass: ";
				std::cin >> (*All).TotalMass;
			}
			totmass = (*All).TotalMass;
			break;
	}
	GasProperties mol = getGasProperties(GASTYPE_MOLECULARCLOUD);
	GasProperties diff = getGasProperties(GASTYPE_DIFFUSECLOUD);
	GasProperties intcm = getGasProperties(GASTYPE_INTERCLOUDMEDIUM);
	GasProperties cor = getGasProperties(GASTYPE_CORONALGAS);
	(*All).Volume = (*All).GasFraction*totmass*constants::Msun_g*
		(mol.massfraction/mol.density + diff.massfraction/diff.density +
		intcm.massfraction/intcm.massfraction + cor.massfraction/cor.density);
	//std::cout << "Total volume: " << (*All).Volume << "\n";
}

void readParamFile(const std::string& filename){
	std::ifstream param;
	param.open(filename);
	char tijdelijk;
	std::string parameter;

	while (std::getline(param, parameter)){
		if(parameter.find("%") != std::string::npos)
			parameter = parameter.substr(0, parameter.find("%"));
		if(parameter.find("TotalTime") != std::string::npos){
			eraseStuff(parameter, "TotalTime");
			(*All).TotalTime = 1.0e6*std::stod(parameter);
		}
		// else if(parameter.find("StarDataFile") != std::string::npos){
		// 	eraseStuff(parameter, "StarDataFile");
		// 	(*All).StarDataFile = parameter;
		// }
		else if(parameter.find("TimeStep") != std::string::npos){
			eraseStuff(parameter, "TimeStep");
			(*All).TimeStep = 1.0e6*std::stod(parameter);
		}
		else if(parameter.find("TimeRes") != std::string::npos){
			eraseStuff(parameter, "TimeRes");
			(*All).TimeRes = 1.0e6*std::stod(parameter);
		}
		else if(parameter.find("SurvivalFraction") != std::string::npos){
			eraseStuff(parameter, "SurvivalFraction");
			(*All).SurvivalFraction = std::stod(parameter);
		}
		else if(parameter.find("IMF") != std::string::npos){
			eraseStuff(parameter, "IMF");
			tijdelijk = parameter[0];
			if(tijdelijk == 'K')
				(*All).Imf = IMF_KROUPA;
			else
				(*All).Imf = IMF_SALPETER;
		}
		else if(parameter.find("MinimumMass") != std::string::npos){
			eraseStuff(parameter, "MinimumMass");
			(*All).MinimumMass = std::stod(parameter);
		}
		else if(parameter.find("MaximumMass") != std::string::npos){
			eraseStuff(parameter, "MaximumMass");
			(*All).MaximumMass = std::stod(parameter);
		}
		else if(parameter.find("MassBins") != std::string::npos){
			eraseStuff(parameter, "MassBins");
			(*All).MassBins = std::stoi(parameter);
		}
		else if(parameter.find("SFRType") != std::string::npos){
			eraseStuff(parameter, "SFRType");
			tijdelijk = parameter[0];
			if(tijdelijk == 'E')
				(*All).SFRType = SFR_EXPONENTIAL;
			else if(tijdelijk == 'C')
				(*All).SFRType = SFR_CONSTANT;
			else if(tijdelijk == 'D')
				(*All).SFRType = SFR_DISK;
			else if(tijdelijk == 'L')
				(*All).SFRType = SFR_LOGNORM;
			else
				(*All).SFRType = SFR_STARBURST;
		}
		else if(parameter.find("Profile") != std::string::npos && !((parameter.find("GasProfileOn")) != std::string::npos)){
			eraseStuff(parameter, "Profile");
			tijdelijk = parameter[0];
			if(tijdelijk == 'H')
				(*All).Profile = PROFILE_HERNQUIST;
			else if(tijdelijk == 'N')
				(*All).Profile = PROFILE_NFW;
			else if(tijdelijk == 'C')
				(*All).Profile = PROFILE_CONSTANT;
			else if(tijdelijk == 'I')
				(*All).Profile = PROFILE_ISOTHERMAL;
			else{
				std::cout << "Parameter file error: not a valid Profile value: " << tijdelijk << ".\n";
				throw std::exception();
			}
		}
		// else if(parameter.find("CNFW") != std::string::npos){
		// 	eraseStuff(parameter, "CNFW");
		// 	(*All).CNFW = std::stod(parameter);
		// }
		else if(parameter.find("DeltaVir") != std::string::npos){
			eraseStuff(parameter, "DeltaVir");
			(*All).DeltaVir = std::stod(parameter);
		}
		else if(parameter.find("GasCoeff") != std::string::npos){
			eraseStuff(parameter, "GasCoeff");
			(*All).GasCoeff = std::stod(parameter);
		}
		else if(parameter.find("LogBins") != std::string::npos){
			eraseStuff(parameter, "LogBins");
			(*All).LogBins = std::stod(parameter);
		}
		else if(parameter.find("TotalMass") != std::string::npos){
			eraseStuff(parameter, "TotalMass");
			(*All).TotalMass = std::stod(parameter)/NTask;
		}
		else if(parameter.find("Radius") != std::string::npos && !((parameter.find("DiskRadius")) != std::string::npos)){
			eraseStuff(parameter, "Radius");
			(*All).Radius = std::stod(parameter);
		}
		else if(parameter.find("BurstWidth") != std::string::npos){
			eraseStuff(parameter, "BurstWidth");
			(*All).BurstWidth = 1.0e6*std::stod(parameter);
		}
		else if(parameter.find("PeakTime") != std::string::npos){
			eraseStuff(parameter, "PeakTime");
			(*All).PeakTime = 1.0e6*std::stod(parameter);
		}
		else if(parameter.find("SFR") != std::string::npos && !((parameter.find("SFRType")) != std::string::npos)){
			eraseStuff(parameter, "SFR");
			(*All).Sfr = std::stod(parameter);
		}
		// else if(parameter.find("GasType") != std::string::npos){
		// 	eraseStuff(parameter, "GasType");
		// 	tijdelijk = parameter[0];
		// 	if(tijdelijk == 'M')
		// 		(*All).GasType = GASTYPE_MOLECULARCLOUD;
		// 	else if(tijdelijk == 'D')
		// 		(*All).GasType = GASTYPE_DIFFUSECLOUD;
		// 	else if(tijdelijk == 'I')
		// 		(*All).GasType = GASTYPE_INTERCLOUDMEDIUM;
		// 	else
		// 		(*All).GasType = GASTYPE_CORONALGAS;
		// }
		else if(parameter.find("GasProfileOn") != std::string::npos){
			eraseStuff(parameter, "GasProfileOn");
			(*All).GasProfileOn = std::stoi(parameter);
		}
		else if(parameter.find("StromgrenOn") != std::string::npos){
			eraseStuff(parameter, "StromgrenOn");
			(*All).StromgrenOn = std::stoi(parameter);
		}
		else if(parameter.find("PopulationOn") != std::string::npos){
			eraseStuff(parameter, "PopulationOn");
			(*All).PopulationOn = std::stoi(parameter);
		}
		else if(parameter.find("RadiationOn") != std::string::npos){
			eraseStuff(parameter, "RadiationOn");
			(*All).RadiationOn = std::stoi(parameter);
		}
		else if(parameter.find("PhotonsOn") != std::string::npos){
			eraseStuff(parameter, "PhotonsOn");
			(*All).PhotonsOn = std::stoi(parameter);
		}
		else if(parameter.find("EnergyBins") != std::string::npos){
			eraseStuff(parameter, "EnergyBins");
			(*All).EnergyBins = std::stoi(parameter);
		}
		else if(parameter.find("EnergyLogOn") != std::string::npos){
			eraseStuff(parameter, "EnergyLogOn");
			(*All).EnergyLogOn = std::stoi(parameter);
		}
		else if(parameter.find("MinimumEnergy") != std::string::npos){
			eraseStuff(parameter, "MinimumEnergy");
			(*All).MinimumEnergy = std::stod(parameter);
		}
		else if(parameter.find("MaximumEnergy") != std::string::npos){
			eraseStuff(parameter, "MaximumEnergy");
			(*All).MaximumEnergy = std::stod(parameter);
		}
		else if(parameter.find("GasFraction") != std::string::npos){
			eraseStuff(parameter, "GasFraction");
			(*All).GasFraction = std::stod(parameter);
		}
		else if(parameter.find("DiskRadius") != std::string::npos){
			eraseStuff(parameter, "DiskRadius");
			(*All).DiskRadius = std::stod(parameter);
			(*All).DiskRadius *= constants::kpc_to_cm;
		}
		else if(parameter.find("DiskBins") != std::string::npos){
			eraseStuff(parameter, "DiskBins");
			(*All).DiskBins = std::stoi(parameter);
		}
	}
	if(!(*All).GasFraction){
		std::cout << "Enter gas fraction: ";
		std::cin >> (*All).GasFraction;
	}
	if(!(*All).TotalTime){
		std::cout << "Enter total time (Myr): ";
		std::cin >> (*All).TotalTime;
		(*All).TotalTime *= 1.0e6;
	}
	if(!(*All).TimeStep){
		std::cout << "Enter time step (Myr): ";
		std::cin >> (*All).TimeStep;
		(*All).TimeStep *= 1.0e6;
	}
	if(!(*All).TimeRes){
		std::cout << "Enter time resolution (Myr): ";
		std::cin >> (*All).TimeRes;
		(*All).TimeRes *= 1.0e6;
	}
	if(!(*All).SurvivalFraction){
		std::cout << "Enter survival fraction: ";
		std::cin >> (*All).SurvivalFraction;
	}
	if(!(*All).MinimumMass){
		std::cout << "Enter minimum mass: ";
		std::cin >> (*All).MinimumMass;
	}
	if(!(*All).MaximumMass){
		std::cout << "Enter maximum mass: ";
		std::cin >> (*All).MaximumMass;
	}
	if(!(*All).MassBins){
		std::cout << "Enter number of mass bins: ";
		std::cin >> (*All).MassBins;
	}
	if(!(*All).SFRType){
		std::cout << "Enter SFR type.\n Choose from Exponential (E), Constant (C), Star Burst (S), or Disk (D): ";
		std::cin >> tijdelijk;
		if(tijdelijk == 'E')
			(*All).SFRType = SFR_EXPONENTIAL;
		else if(tijdelijk == 'C')
			(*All).SFRType = SFR_CONSTANT;
		else if(tijdelijk == 'D')
			(*All).SFRType = SFR_DISK;
		else
			(*All).SFRType = SFR_STARBURST;
	}
	if(!(*All).Imf){
		(*All).Imf = IMF_KROUPA;
	}

	switch((*All).SFRType){
		case SFR_EXPONENTIAL:
			if(!(*All).TotalMass){
				std::cout << "Enter total mass (Msun): ";
				std::cin >> (*All).TotalMass;
			}
			if(!(*All).BurstWidth){
				std::cout << "Enter burst width (Myr): ";
				std::cin >> (*All).BurstWidth;
				(*All).BurstWidth *= 1.e6;
			}
			break;
		case SFR_CONSTANT:
			if(!(*All).Sfr){
				std::cout << "Enter SFR (yr^-1): ";
				std::cin >> (*All).Sfr;
			}
			break;
		case SFR_STARBURST:
			if(!(*All).TotalMass){
				std::cout << "Enter total mass (Msun): ";
				std::cin >> (*All).TotalMass;
			}
			break;
		case SFR_DISK:
			if(!(*All).DiskRadius){
				std::cout << "Enter disk radius (kpc): ";
				std::cin >> (*All).DiskRadius;
				(*All).DiskRadius *= constants::kpc_to_cm;
			}
			if(!(*All).DiskBins){
				std::cout << "Enter number of disk bins: ";
				std::cin >> (*All).DiskBins;
			}
			break;
		default:
			if(!(*All).TotalMass){
				std::cout << "Enter total mass (Msun): ";
				std::cin >> (*All).TotalMass;
			}
			break;
	}
	// if(!(*All).GasType){
	// 	std::cout << "Enter Gas Type. Molecular cloud (M), Diffuse cloud (D), Intercloud medium (I), Coronal gas (C): ";
	// 	std::cin >> tijdelijk;
	// 	if(tijdelijk == 'M')
	// 		(*All).GasType = GASTYPE_MOLECULARCLOUD;
	// 	else if(tijdelijk == 'D')
	// 		(*All).GasType = GASTYPE_DIFFUSECLOUD;
	// 	else if(tijdelijk == 'I')
	// 		(*All).GasType = GASTYPE_INTERCLOUDMEDIUM;
	// 	else
	// 		(*All).GasType = GASTYPE_CORONALGAS;
	// }

	if((*All).RadiationOn != 0 && (*All).RadiationOn != 1){
		std::cout << "Radiation on (1) or off (0): ";
		std::cin >> (*All).RadiationOn;
	}

	if(!(*All).EnergyBins){
		std::cout << "Number of Energy bins: ";
		std::cin >> (*All).EnergyBins;
	}

	if(!(*All).MinimumEnergy){
		std::cout << "Minimum energy: ";
		std::cin >> (*All).MinimumEnergy;
	}

	if(!(*All).MaximumEnergy){
		std::cout << "Maximum energy: ";
		std::cin >> (*All).MaximumEnergy;
	}
	setVolume();
}

void AllParticles::readRecombRates(){
	std::ifstream infile("InputFiles/recombVerner.txt");
	std::array<double, 3> a, b, T0, T1;
	for(int j = 0; j < 4; j++){
		std::string strInput;
		int Z;
		int N;
		if(j == 0)
			getline(infile, strInput);
		else{
			infile >> strInput;
			infile >> Z;
			infile >> N;
			infile >> a[j-1];
			infile >> b[j-1];
			infile >> T0[j-1];
			infile >> T1[j-1];
		}
	}
	m_H.a = a[0];
	m_H.b = b[0];
	m_H.T0 = T0[0];
	m_H.T1 = T1[0];
	m_He.a = a[1];
	m_He.b = b[1];
	m_He.T0 = T0[1];
	m_He.T1 = T1[1];
	m_HeI.a = a[2];
	m_HeI.b = b[2];
	m_HeI.T0 = T0[2];
	m_HeI.T1 = T1[2];
}

void AllParticles::readCoolingFunction(){
	std::ifstream infile("InputFiles/TemperatureCooling.txt");
	for(int j = 0; j < 701; j++){
		infile >> m_TLam[j];
		infile >> m_Lambda[j];
	}
}
void GasProfile::readCoolingFunction(){
	std::ifstream infile("InputFiles/TemperatureCooling.txt");
	for(int j = 0; j < 701; j++){
		infile >> m_TLam[j];
		infile >> m_Lambda[j];
	}
}

void AllParticles::readElectronCrossSection(){}
void AllParticles::readPhotoCrossSections(){ //Verner 1996
	std::ifstream infile("InputFiles/photocrs.txt");
	std::array<double, 3> Eth, Emax, E0, sigma0, ya, P, yw, y0, y1;
	for(int j = 0; j < 4; j++){
		std::string strInput;
		int Z;
		int N;
		if(j == 0)
			getline(infile, strInput);
		else{
			infile >> strInput;
			infile >> Z;
			infile >> N;
			infile >> Eth[j-1];
			infile >> Emax[j-1];
			infile >> E0[j-1];
			infile >> sigma0[j-1];
			infile >> ya[j-1];
			infile >> P[j-1];
			infile >> yw[j-1];
			infile >> y0[j-1];
			infile >> y1[j-1];
		}
	}
	for(int j = 0; j < m_E.size(); j++){
		std::array<double, 3> F;
		for(int i = 0; i < 3; i++){
			double x{m_E[j]/E0[i] - y0[i]};
			double y{sqrt(x*x + y1[i]*y1[i])};
			if(m_E[j] > E0[i] && m_E[j] < Emax[i]){
				F[i] = ((x - 1.0)*(x - 1.0) + yw[i]*yw[i])*pow(y, (0.5*P[i] - 5.5))*pow((1.0 + sqrt(y/ya[i])),-P[i]);
			}
		}
		m_crossSection[j].H = sigma0[0]*F[0]*1.0e-18; //cm2
		//std::cout << m_crossSection[j].H  << "\n";
		m_crossSection[j].He = sigma0[1]*F[1]*1.0e-18; //cm2
		m_crossSection[j].HeI = sigma0[2]*F[2]*1.0e-18; //cm2
		Gas::s_Eth.H = Eth[0];
		Gas::s_Eth.He = Eth[1];
		Gas::s_Eth.HeI = Eth[2];
	}
}

void GasProfile::readCNFW(double mass){ //From Aaron: Planck_cmz_z0z1z2z3.dat
	/* Obtain righ c_NFW for the given virial mass */
	std::ifstream infile("InputFiles/Planck_cmz_z0z1z2z3.dat");
	double z, m, cnfw, skip1, skip2;
	double m1, m2, cnfw1, cnfw2;
	double massvir{std::log10(mass)};
	for(int j = 0; j < 103; j++){
		infile >> z;
		if(z >= 1.0){
			(*All).CNFW = cnfw;
			std::cout << "C_NFW = " << (*All).CNFW << "\n";
			break;
		}
		infile >> m;
		m += 10; //data is devided by 1e10
		infile >> cnfw;
		infile >> skip1;
		infile >> skip2;
		if(m == massvir){
			(*All).CNFW = cnfw;
			std::cout << "C_NFW = " << (*All).CNFW << "\n";
			break;
		}
		if(m < massvir || j == 0){
			m1 = m;
			cnfw1 = cnfw;
		}
		if(m > massvir){
			m2 = m;
			cnfw2 = cnfw;
			(*All).CNFW = (cnfw2 - cnfw1) * ((mass) - pow(10, m1))/(pow(10, m2) - pow(10, m1)) + cnfw1;
			std::cout << "C_NFW = " << (*All).CNFW << "\n";
			break;
		}
	}	
}