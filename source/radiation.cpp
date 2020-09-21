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
/*-----------------Static Declarations----------------*/

double AllParticles::s_allTime = 0;
//GasProperties *gasProperties = new GasProperties;


/* ----------------Class functions -------------------*/
double AllParticles::equilfunction(double leftfac, double r, double nH, std::vector<double> N_gamma, double x){
	double righthand{0}, lefthand{leftfac*x*x};
	for(int i = 0; i < N_gamma.size(); i++){
		righthand += N_gamma[i]*(1.0 - exp(-m_crossSection[i].H*(nH - x)*r));
	}	
	return lefthand - righthand;
}
double AllParticles::findroot(double leftfac, double r, double nH, std::vector<double> N_gamma, double x1, double x2){
	double verschill{equilfunction(leftfac, r, nH, N_gamma, x1)}, verschilr{equilfunction(leftfac, r, nH, N_gamma, x2)};
	if((verschill < 0) && (verschilr < 0)){
		std::cout << "Findroot error\n";
		throw std::exception();
	}
	if(verschill > 0){
		std::cout << "Findroot error\n";
		throw std::exception();
	}
	if((verschill < 0) && (verschilr > 0)){
		double xdelta{(x2 - x1)/4.0};
		if(xdelta < 1.0e-15){
			return x1;
		}
		double verschil{0};
		for(int i = 1; i < 6; i++){
			verschil = equilfunction(leftfac, r, nH, N_gamma, x1 + i*xdelta);
			if(verschil > 0){
				if(verschil < 0.001)
					return x1 + i*xdelta;
				else
					return findroot(leftfac, r, nH, N_gamma, x1 + (i-1)*xdelta, x1 + i*xdelta);
			}
		}
	}
}

double AllParticles::getEquilibriumHfrac(GasProperties gasproperties, std::vector<double> N_gamma, double stars){
	double r{pow(3.0*(*All).Volume/stars/4.0/M_PI, 0.3333)};
	double alpha{recombinationCoefficient(gasproperties.temperature).H};
	double N_gammacrs{0}, crs{0}, nH{gasproperties.density/constants::hydrogen_g};
	double righthand{0}, lefthand{0};
	for(int i = 0; i < (*All).EnergyBins; i++){
		N_gammacrs += N_gamma[i] * m_crossSection[i].H;
		crs += m_crossSection[i].H;
	}
	if(crs*gasproperties.density/constants::hydrogen_g*r > 0.001){
		return findroot(alpha*4.0/3.0*M_PI*r*r*r, r, gasproperties.density/constants::hydrogen_g, N_gamma, 0, nH);
	}
	else{
		double factor{8.0*M_PI*r*r*alpha/3.0/N_gammacrs};
		if(factor == 0 || N_gammacrs == 0){
			std::cout<< "Zero encountered in getEquilibriumHfrac; alphaR: " << alpha <<", "<< "Temperature: " << gasproperties.temperature<<"\n";
			throw std::exception();
		}
		double equil{(sqrt(1.0+2.0*factor*gasproperties.density/constants::hydrogen_g) - 1.0)/factor};
		if(equil > gasproperties.density/constants::hydrogen_g && equil < gasproperties.density/constants::hydrogen_g*1.0001)
			equil = gasproperties.density/constants::hydrogen_g;
		if(equil > gasproperties.density/constants::hydrogen_g || equil < 0.0){
			std::cout << "Kan nie kloppen (getEquilibriumHfrac): " << equil << ", " << factor << ", " << gasproperties.density/constants::hydrogen_g << "\n";
			throw std::exception();
		}
		return equil;
	}
}
void AllParticles::getEscapeFraction(GasProperties gasproperties){//, std::vector<double> N_gamma){
	double totalPhotons{0}, minphotons{0};
	double equilibrium{getEquilibriumHfrac(gasproperties, m_Nphoton, 1.0/gasproperties.volume)};
	double r{pow(3.0*(*All).Volume*gasproperties.volume/4.0/M_PI, 0.3333)};
	for(int i = 0; i < m_Nphoton.size(); i++){
		totalPhotons += m_Nphoton[i];
	}
	for(int i = 0; i < m_Nphoton.size(); i++){
		minphotons = m_Nphoton[i]*(1.0 - 
			exp(-m_crossSection[i].H*(gasproperties.density/constants::hydrogen_g - equilibrium)*r));
		if(minphotons > m_Nphoton[i]){
			std::cout << "This shouldn't be possible (getEscapeFraction): " << minphotons << ", " << m_Nphoton[i] << "\n";
			throw std::exception();
		}
		m_Nphoton[i] -= minphotons;
		m_escapefraction -= (minphotons/m_NphotonTotal);
		if(m_escapefraction < 0)
			m_escapefraction = 0;
	}
}
std::vector<double> AllParticles::lostPhotons(GasProperties gasprop, std::vector<double> Nphoton, double stars){
	std::vector<double> photons{multiplyDoubleVector(Nphoton, gasprop.volume)};
	double eq{getEquilibriumHfrac(gasprop, photons, stars)};
	double r{pow(3.0*(*All).Volume/stars/4.0/M_PI, 0.3333)};
	for(int i = 0; i < Nphoton.size(); i++){
		photons[i] -= photons[i]*(1.0 - exp(-m_crossSection[i].H*(gasprop.density/constants::hydrogen_g - eq)*r));
		if(photons[i] < 0){
			std::cout << "This shouldn't be possible (lostPhotons): " << photons[i] << ", " << Nphoton[i] << "\n";
			throw std::exception();
		}
	}
	return photons;
}
double AllParticles::equilStromgren(double density, double alpha, double x){
	double deltaR{(*All).DiskRadius/(*All).DiskBins};
	double lefthand{sumDoubleVector(m_Nphoton) * deltaR * 0.5 * (x*x - 0.0833*deltaR*deltaR) / (1.333*x*x*x)};
	double righthand{alpha * density * density * M_PI * m_r * 
		( deltaR * sqrt( x*x - 0.25*deltaR*deltaR ) + 2.0 * x * x * asin( deltaR/2.0/x ) )};
	//std::cout << alpha << "\t" << deltaR << "\t" << x << "\n";
	return lefthand - righthand;	
}
double AllParticles::findStromgrenRad(double density, double alpha){
	double deltaR{(*All).DiskRadius/static_cast<double>((*All).DiskBins)};
	double start{equilStromgren(density, alpha, deltaR)};
	double end{equilStromgren(density, alpha, deltaR*1.0001)};
	double x2{deltaR};
	if(end*end > start*start){
		return sqrt(sumDoubleVector(m_Nphoton)/2.0/alpha/m_r)/M_PI/density;
	}
	else{
		while(start*end > 0){
			x2 *= 2.0;
			end = equilStromgren(density, alpha, x2);
		}
		return findStromgrenRad(density, alpha, 0.5*x2, x2);
	}
}

double AllParticles::findStromgrenRad(double density, double alpha, double x1, double x2){
	double verschill{equilStromgren(density, alpha, x1)}, verschilr{equilStromgren(density, alpha, x2)};
	//std::cout << x1  << "\t" << verschill << "\t" << x2 << "\t" << verschilr << "\t" << density << "\n";
	if(verschilr == 0 && x2 > 1.0)
		return x2;
	else if(verschill == 0 && x1 > 1.0)
		return x1;
	else if(verschill*verschilr > 0)
		throw std::exception();
	else{
		double xdelta{(x2 - x1)/4.0};
		if(xdelta < 1.0e-15*x1){
			return x2;
		}
		double verschil{0};
		for(int i = 1; i < 6; i++){
			verschil = equilStromgren(density, alpha, x1 + i*xdelta);
			if(verschil*verschill < 0){
				if(verschil*verschil < 1.0)
					return x2;
				else
					return findStromgrenRad(density, alpha, x1 + (i - 1.0)*xdelta, x1 + i*xdelta);
			}
		}
	}
}
double AllParticles::getStromgrenRadiusPerGasType(GASTYPE gastype){
	double density{0}, alpha{0};
	switch(gastype){
		case GASTYPE_MOLECULARCLOUD:
			density = getGasProperties(GASTYPE_MOLECULARCLOUD).density/constants::hydrogen_g;
			alpha = recombinationCoefficient(getGasProperties(GASTYPE_MOLECULARCLOUD).temperature).H;
			break;
		case GASTYPE_DIFFUSECLOUD:
			density = getGasProperties(GASTYPE_DIFFUSECLOUD).density/constants::hydrogen_g;
			alpha = recombinationCoefficient(getGasProperties(GASTYPE_DIFFUSECLOUD).temperature).H;
			break;
		case GASTYPE_INTERCLOUDMEDIUM:
			density = getGasProperties(GASTYPE_INTERCLOUDMEDIUM).density/constants::hydrogen_g;
			alpha = recombinationCoefficient(getGasProperties(GASTYPE_INTERCLOUDMEDIUM).temperature).H;
			break;
		case GASTYPE_CORONALGAS:
			density = getGasProperties(GASTYPE_CORONALGAS).density/constants::hydrogen_g;
			alpha = recombinationCoefficient(getGasProperties(GASTYPE_CORONALGAS).temperature).H;
			break;
		default:
			density = getGasProperties(GASTYPE_CORONALGAS).density/constants::hydrogen_g;
			alpha = recombinationCoefficient(getGasProperties(GASTYPE_CORONALGAS).temperature).H;
			break;
	}
	if(m_r < 0.6204){
		return 0.0;
	}
	else{
		double StromgrenRad{findStromgrenRad(density, alpha)};
		//std::cout << StromgrenRad << "\n";
		return StromgrenRad;
	}
}
void AllParticles::getEscFracPerGasBubble(){//, std::vector<double> N_gamma){
	double totalPhotons{0}, minphotons{0};
	for(int i = 0; i < m_Nphoton.size(); i++){
		totalPhotons += m_Nphoton[i];
	}
	std::vector<double> lostphotonsM{lostPhotons(getGasProperties(GASTYPE_MOLECULARCLOUD), 
			m_Nphoton, 1.0/getGasProperties(GASTYPE_MOLECULARCLOUD).volume)};
	std::vector<double> lostphotonsD{lostPhotons(getGasProperties(GASTYPE_DIFFUSECLOUD), 
			m_Nphoton, 1.0/getGasProperties(GASTYPE_DIFFUSECLOUD).volume)};
	std::vector<double> lostphotonsI{lostPhotons(getGasProperties(GASTYPE_INTERCLOUDMEDIUM), 
			m_Nphoton, 1.0/getGasProperties(GASTYPE_INTERCLOUDMEDIUM).volume)};
	std::vector<double> lostphotonsC{lostPhotons(getGasProperties(GASTYPE_CORONALGAS), 
			m_Nphoton, 1.0/getGasProperties(GASTYPE_CORONALGAS).volume)};

	for(int i = 0; i < m_Nphoton.size(); i++){
		m_Nphoton[i] = lostphotonsM[i] + lostphotonsD[i] + lostphotonsI[i] + lostphotonsC[i];
	}
	m_escapefraction = sumDoubleVector(m_Nphoton)/m_NphotonTotal;
}
void AllParticles::getEscFracPerStar(){
    std::vector<double> photM(m_Nphoton.size()), photD(m_Nphoton.size()), photI(m_Nphoton.size()), photC(m_Nphoton.size()), photTot(m_Nphoton.size());
    writeStarPhotonsPerGasType();
    photonsStar(photM, "StarPhotonsM" + std::to_string(ThisTask) + ".txt");
    photonsStar(photD, "StarPhotonsD" + std::to_string(ThisTask) + ".txt");
    photonsStar(photI, "StarPhotonsI" + std::to_string(ThisTask) + ".txt");
    photonsStar(photC, "StarPhotonsC" + std::to_string(ThisTask) + ".txt");
    photM = multiplyDoubleVector(photM, getGasProperties(GASTYPE_MOLECULARCLOUD).volume);
    photD = multiplyDoubleVector(photD, getGasProperties(GASTYPE_DIFFUSECLOUD).volume);
    photI = multiplyDoubleVector(photI, getGasProperties(GASTYPE_INTERCLOUDMEDIUM).volume);
    photC = multiplyDoubleVector(photC, getGasProperties(GASTYPE_CORONALGAS).volume);
    for(int i = 0; i < m_Nphoton.size(); i++){
        m_Nphoton[i] = photM[i] + photD[i] + photI[i] + photC[i];
    }
    m_escapefraction = sumDoubleVector(m_Nphoton)/m_NionStars;
}
//AllParticles

void AllParticles::assignEnergy(){
	m_E.resize((*All).EnergyBins);
	m_Nphoton.resize((*All).EnergyBins);
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

void AllParticles::measureDistanceGasToStars(){
	for(int i = 0; i < Particle::s_NumGas; i++){
		m_gas[i].measureDistanceToStars(m_star);
	}
}
void AllParticles::calculatePhotonNumber(){
	for(int i = 0; i < Particle::s_NumStar; i++){
		m_star[i].photons(0) *= m_timestep;
	}
}
std::vector<int> AllParticles::orderDistanceToStar(int starID){
	std::vector<int> orderedArray(Particle::s_NumGas);
	for(int i = 0; i < Particle::s_NumGas; i++){
		orderedArray[i] = i;
	}
	for(int startIndex = 0; startIndex < Particle::s_NumGas; startIndex++){
		int bestIndex = startIndex;
		for(int currentIndex = startIndex + 1; currentIndex < Particle::s_NumGas; currentIndex++){
			if(m_gas[orderedArray[bestIndex]].distanceStar(starID) > m_gas[orderedArray[currentIndex]].distanceStar(starID)){
				bestIndex = currentIndex;
			}
		}
		std::swap(orderedArray[startIndex], orderedArray[bestIndex]);
	}
	return orderedArray;
}
std::vector<int> AllParticles::orderPhotons(int starID){
	std::vector<int> orderedArray(m_E.size());
	for(int i = 0; i < m_E.size(); i++){
		orderedArray[i] = i;
	}
	for(int startIndex = 0; startIndex < m_E.size(); startIndex++){
		int bestIndex = startIndex;
		for(int currentIndex = startIndex + 1; currentIndex < m_E.size(); currentIndex++){
			if(m_star[starID].photons(orderedArray[bestIndex]) < m_star[starID].photons(orderedArray[currentIndex])){
				bestIndex = currentIndex;
			}
		}
		std::swap(orderedArray[startIndex], orderedArray[bestIndex]);
	}
	return orderedArray;
}
double numberOfIonizations(double photons, double atoms, double crs, double A){
	double ionizations{photons*atoms*crs/A};
	if(10.0*ionizations > atoms){
		ionizations = numberOfIonizations(0.5*photons, atoms, crs, A);
		ionizations += numberOfIonizations(0.5*photons, atoms-ionizations, crs, A);
	}
	return ionizations;
}

double AllParticles::findGoodDividingFactor(Gas gas, Star star, int photonIndex){
	double photons{star.photons(photonIndex)};
	double chance;
	double chanceH{m_crossSection[photonIndex].H/(4.0*M_PI*gas.distanceStar(star.ID())*gas.distanceStar(star.ID()))};
	if(gas.speciesNum().He + gas.speciesNum().HeI + gas.speciesNum().HeII == 0)
		chance = chanceH;
	else{
		double chanceHe{m_crossSection[photonIndex].He/(4.0*M_PI*gas.distanceStar(star.ID())*gas.distanceStar(star.ID()))};
		double chanceHeI{m_crossSection[photonIndex].HeI/(4.0*M_PI*gas.distanceStar(star.ID())*gas.distanceStar(star.ID()))};
		chance = std::max(chanceH, chanceHe);
		chance = std::max(chance, chanceHeI);
	}
	if(photons > 1.0e3 && chance < 1.0e-3)
		return std::min((1.0e-3/(photons*chance)), 1.0);
}
void AllParticles::chanceOfIonizing(Gas& gas, Star& star, int photonIndex){
	int total{ (gas.speciesNum().H > gas.speciesNum().He) ? 
		static_cast<int>(gas.speciesNum().H) : static_cast<int>(gas.speciesNum().He) };
	total = (total > gas.speciesNum().HeI) ? total : static_cast<int>(gas.speciesNum().HeI);
	double A{4.0*M_PI*gas.distanceStar(star.ID())*gas.distanceStar(star.ID())};

	//H = numberOfIonizations(star.photons(photonIndex), gas.speciesNum().H, m_crossSection[photonIndex].H, A);
	//He = numberOfIonizatsions(star.photons(photonIndex), gas.speciesNum().He, m_crossSection[photonIndex].He, A);
	//HeI = numberOfIonizations(star.photons(photonIndex), gas.speciesNum().HeI + He, m_crossSection[photonIndex].HeI, A);
	Atoms alphaR{recombinationCoefficient(gas.Energy()/gas.speciesNum().X*80.8)};
	double ne{(gas.speciesNum().HI + gas.speciesNum().HeI + 2.0*gas.speciesNum().HeII)*gas.denmass()};
	double recombH{0}, ionH{0}, recombHe{0}, ionHe{0}, recombHeI{0}, ionHeI{0};
	double photons{0};
	double dividingFactor{0}; //To make it faster
	dividingFactor = findGoodDividingFactor(gas, star, photonIndex);
	for(int i = 0; i < static_cast<int>(1.0/dividingFactor); i++){

		ionH = dividingFactor*star.photons(photonIndex)*(gas.speciesNum().H)*m_crossSection[photonIndex].H/A;
		ionHe = dividingFactor*star.photons(photonIndex)*(gas.speciesNum().He)*m_crossSection[photonIndex].He/A;
		ionHeI = dividingFactor*star.photons(photonIndex)*(gas.speciesNum().HeI)*m_crossSection[photonIndex].HeI/A;
		photons += ionH + ionHe + ionHeI;

		gas.changeNum(ionH, ionHe, ionHeI, m_E[photonIndex]);

		ne = (gas.speciesNum().HI + gas.speciesNum().HeI + 2.0*gas.speciesNum().HeII)*gas.denmass();

		alphaR = recombinationCoefficient(gas.Energy()/gas.speciesNum().X*80.8);
		recombH = std::min(1.0e10*alphaR.H*ne*ionH, gas.speciesNum().HI);
		recombHe = std::min(alphaR.He*ne*ionHe, gas.speciesNum().HeI);
		recombHeI = std::min(alphaR.HeI*ne*ionHeI, gas.speciesNum().HeII);

		gas.changeNum(-recombH, -recombHe, -recombHeI, m_E[photonIndex]);
		
		if(photons > (star.photons(photonIndex)+1))
			break;
	}
	//std::cout << gas.Energy()/gas.speciesNum().X << ", " << ionH << ", " << recombH << "\n";
	star.photons(photonIndex) -= photons;
	//gas.changeNum(H, He, HeI, m_E[photonIndex]);
	//Atoms alphaR{recombinationCoefficient(gas.Energy()/gas.speciesNum().X*80.8)};
	//gas.recombinations(alphaR, m_timestep);
	if(gas.speciesNum().H == 0 && gas.speciesNum().He == 0 && gas.speciesNum().HeI == 0 && gas.ionized() == false){
		gas.ionized() = true;
		Gas::s_NumIonized += 1;
	}
}
void AllParticles::setIonizationAll(){
	for(int i = 0; i < Particle::s_NumStar; i++){
		std::vector<int> orderedGas = orderDistanceToStar(i);
		for(int j = 0; j < Particle::s_NumGas; j++){
			if (m_gas[orderedGas[j]].ionized() == true)
				continue;
			else if (m_gas[orderedGas[j]].distanceStar(i) > m_star[i].radiationRadius())
				break;
			else{
				std::vector<int> photonIndex{orderPhotons(i)};
				for(int k = 0; k < m_E.size(); k++){
					chanceOfIonizing(m_gas[orderedGas[j]], m_star[i], photonIndex[k]);
				}
			}
		}
	}
	m_totalIonized.H = 0;
	m_totalIonized.He = 0;
	m_totalIonized.HeI = 0;
	for(int j = 0; j < Particle::s_NumGas; j++){
		m_totalIonized.H += m_gas[j].speciesNum().HI/(m_gas[j].speciesNum().H + m_gas[j].speciesNum().HI);
		m_totalIonized.He += m_gas[j].speciesNum().HeI/(m_gas[j].speciesNum().He + m_gas[j].speciesNum().HeI + m_gas[j].speciesNum().HeII);
		m_totalIonized.HeI += m_gas[j].speciesNum().HeII/(m_gas[j].speciesNum().He + m_gas[j].speciesNum().HeI + m_gas[j].speciesNum().HeII);
	}
	m_totalIonized.H /= Particle::s_NumGas;
	m_totalIonized.He /= Particle::s_NumGas;
	m_totalIonized.HeI /= Particle::s_NumGas;
}

void AllParticles::luminosityStarHMXB(){
    m_lumStars = 0;
    m_lumHMXBs = 0;
    m_lumStarsIon = 0;
    m_NionStars = 0;
    m_NionHMXBs = 0;
    m_lumBH = 0;
    m_lumNS = 0;
    for(int i = 0; i < m_star.size(); i++){
        switch(m_star[i].type()){
            case Star::TYPE_STAR:
                m_lumStars += m_star[i].luminosity();
                m_lumStarsIon += m_star[i].ionluminosity();
                m_NionStars += m_star[i].ionN();
                break;
            case Star::TYPE_HMXB:
                m_lumHMXBs += m_star[i].HMXBluminosity();
                m_NionHMXBs += sumDoubleVector(m_star[i].HMXBionN());
                break;
            case Star::TYPE_BH:
                m_lumBH += m_star[i].accluminosity();
                break;
            case Star::TYPE_NS:
            	m_lumNS += m_star[i].accluminosity();
            	break;
            default:
                break;
        }
    }
}

double AllParticles::getTemperature(){
    double averageT{0}, totalMass{0};
    for(int i = 0; i < Particle::s_NumGas; i++){
        averageT += m_gas[i].temperature()*m_gas[i].mass();
        totalMass += m_gas[i].mass();
    }
    return averageT/totalMass;
}

Atoms AllParticles::recombinationCoefficient(double T){ //Verner 1996
	Atoms recombCoeff;
	if(T == 0){
		recombCoeff.H = 0;
		recombCoeff.He = 0;
		recombCoeff.HeI = 0;
		return recombCoeff;
	}
	else{
		recombCoeff.H = m_H.a/(sqrt(T/m_H.T0)*
			pow((1.0 + sqrt(T/m_H.T0)), (1.0 - m_H.b))*
			pow((1.0 + sqrt(T/m_H.T1)), (1.0 + m_H.b)));
		recombCoeff.He = m_He.a/(sqrt(T/m_He.T0)*
			pow((1.0 + sqrt(T/m_He.T0)), (1.0 - m_He.b))*
			pow((1.0 + sqrt(T/m_He.T1)), (1.0 + m_He.b)));
		recombCoeff.HeI = m_HeI.a/(sqrt(T/m_HeI.T0)*
			pow((1.0 + sqrt(T/m_HeI.T0)), (1.0 - m_HeI.b))*
			pow((1.0 + sqrt(T/m_HeI.T1)), (1.0 + m_HeI.b)));
		return recombCoeff;
	}
}
// void AllParticles::assignGasProperties(){
// 	switch((*All).GasType){
// 		case GASTYPE_MOLECULARCLOUD:
// 			(*gasProperties).volume = 0.005;
// 			(*gasProperties).massfraction = 0.4;
// 			(*gasProperties).density = 1000.0*constants::hydrogen_g;
// 			(*gasProperties).temperature = 20.0;
// 			(*gasProperties).c_sound = 0.6e5; //cm/s
// 			break;
// 		case GASTYPE_DIFFUSECLOUD:
// 			(*gasProperties).volume = 0.05;
// 			(*gasProperties).massfraction = 0.4;
// 			(*gasProperties).density = 100.0*constants::hydrogen_g;
// 			(*gasProperties).temperature = 80.0;
// 			(*gasProperties).c_sound = 0.9e5; //cm/s
// 			break;
// 		case GASTYPE_INTERCLOUDMEDIUM:
// 			(*gasProperties).volume = 0.4;
// 			(*gasProperties).massfraction = 0.2;
// 			(*gasProperties).density = 1.0*constants::hydrogen_g;
// 			(*gasProperties).temperature = 8000.0;
// 			(*gasProperties).c_sound = 9.0e5; //cm/s
// 			break;
// 		case GASTYPE_CORONALGAS:
// 			(*gasProperties).volume = 0.5;
// 			(*gasProperties).massfraction = 0.001;
// 			(*gasProperties).density = 0.001*constants::hydrogen_g;
// 			(*gasProperties).temperature = 10.0e6;
// 			(*gasProperties).c_sound = 100.0e5; //cm/s
// 			break;
// 		default:
// 			break;
// 	}
// }

