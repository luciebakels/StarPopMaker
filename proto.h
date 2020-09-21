#ifndef PROTO_H
#define PROTO_H

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <random>
#include "allvars.h"
#include "mpi.h"

//Some simple equations from general.cpp
double hernquistFactor();
double totalHernquistMass();
double virialRadius(double mass);
double hernquistScalingRadius();


//Reading parameterfile and getting the total volume of the halo
void readParamFile(const std::string& filename);
void eraseStuff(std::string &mystring, std::string erasethis);
void setVolume();

//Very slow, but accurate random number generator
double randomNumber(double min, double max);
//Useless
std::array<int, 1000> thousandRandomIndex(double max);
//Don't think I'm using these anymore
double probabilityFunction(int X, double atoms, double photons, double crs, double A);
double combinations(double N, double x);
double stirlingApprox(double n);
double numberOfIonizations(double photons, double atoms, double crs, double A);
//Fast way for multiplying, summing and writing vectors
std::vector<double> multiplyDoubleVector(std::vector<double> vec, double val);
double sumDoubleVector(std::vector<double> vec);
void writeVector(std::vector<std::vector<double>> outputvector, std::string fileName, char x);

//IMFs
double kroupaIMF(double randval);
double salpeterIMF(double randval);

//Integrating a random function (sometimes fails)
double integrate(double lowlim, double uplim, double (*function)(double));

//Stuff that can be set in the parameter file
enum SFR {
	SFR_EXPONENTIAL = 1,
	SFR_STARBURST,
	SFR_CONSTANT,
	SFR_DISK,
	SFR_LOGNORM
};

enum IMF{
    IMF_KROUPA = 1,
    IMF_SALPETER,
    IMF_CHABRIER
};

enum GASTYPE {
	GASTYPE_MOLECULARCLOUD = 1,
	GASTYPE_DIFFUSECLOUD,
	GASTYPE_INTERCLOUDMEDIUM,
	GASTYPE_CORONALGAS
};

enum PROFILE {
	PROFILE_HERNQUIST = 1,
	PROFILE_NFW,
	PROFILE_CONSTANT,
	PROFILE_ISOTHERMAL
};

//Composition gas and properties
struct Atoms {
	double H = 0;
	double He = 0;
	double HeI = 0;
};

struct GasProperties {
	double volume;
	double massfraction;
	double density;
	double temperature;
	double c_sound;
};

//All parameters
struct Parameters {
	//std::string StarDataFile{"StarData.txt"};
	double TotalTime;
	double TimeStep;
	double TimeRes;
	double SurvivalFraction;
	IMF Imf;
	double MinimumMass;
	double MaximumMass;
	int MassBins;
	SFR SFRType;
	PROFILE Profile;
	double CNFW;
	double DeltaVir;
	double GasCoeff;
	double LogBins;
	double TotalMass;
	double Radius;
	double BurstWidth;
	double PeakTime;
	double Sfr;
	//GASTYPE GasType;
	int GasProfileOn;
	int StromgrenOn;
	int PopulationOn;
	int RadiationOn;
	int PhotonsOn;
	int EnergyBins;
	int EnergyLogOn;
	double MinimumEnergy;
	double MaximumEnergy;
	double GasFraction;
	double Volume;
	double DiskRadius;
	int DiskBins;
};

//extern GasProperties *gasProperties;
extern Parameters *All;

GASTYPE getGasType();
GasProperties getGasProperties(GASTYPE gastype);

//Class used by particles to set position and velocity, but not inherited
class Position
{
protected:
	std::array<double, 3> m_Pos;

public:
	Position(double x = 0.0, double y = 0.0, double z = 0.0)
		: m_Pos{{x, y, z}} {}
	Position& operator= (const Position &pos);
	Position& operator= (std::array<double, 3> pos);
	double& operator[](int row);
	friend std::ostream& operator<<(std::ostream& out, const Position &pos);
	Position getPosition(){ return *this; }
	void setPosition(double x = 0.0, double y = 0.0, double z = 0.0){ m_Pos = {{x, y, z}}; }
	double distance(){ return sqrt(m_Pos[0]*m_Pos[0] + m_Pos[1]*m_Pos[1] + m_Pos[2]*m_Pos[2]); }

};

//Properties that all particles have: position, velocity, ID, and distance between particles
class Particle
{
protected:
	Position m_Pos;
	Position m_Vel;
	int m_id;
public:
	static int s_NumStar;
	static int s_NumGas;
	static int s_NumParticle;
	Particle(double x = 0.0, double y = 0.0, double z = 0.0){
		m_Pos.setPosition(x, y, z);
		m_id = s_NumParticle++;
	}
	Position& Pos(){ return m_Pos; }
	void Pos(double x, double y, double z){	m_Pos.setPosition(x, y, z); }
	int ID(){ return m_id; }
	double distanceBetween(Particle x);
	friend double distanceBetween(Particle x, Particle y);
};

//Type of particle: includes stars, HMXBs, LMXBs, NSs, BHs, and WDs
class Star : public Particle
{
public:
	enum TYPE{
		TYPE_STAR,
		TYPE_HMXB,
		TYPE_LMXB,
		TYPE_NS,
		TYPE_BH,
		TYPE_WD
	};
private:
	double m_mass, m_L, m_T, m_tms, m_thmxb, m_twd, m_A, m_R, m_Lion, m_Nion, m_Lacc, m_LWD, m_LWD0, m_LHMXB, m_LIonHMXB, m_THMXB;
	std::vector<double> m_E, m_NHMXB;
	double m_radiationRadius = 1.0;
	bool m_active = true, m_SN = false;
	double m_massCompanion, m_massRemnant;
	std::array<double, 2> m_alphabeta;
	std::array<double, 2> m_photoncount;
	std::vector<double> m_photons{0};
	TYPE m_type;
	GASTYPE m_gastype;

public: //Needed for constructor
	void alphabeta(); //Power(2013). Only stars between 0.5 and 100
	double lifetime(double mass);
	void assignProperties(){
		assignRadius();
		assignLuminosity();
		assignTemperature();
	}
	void assignRadius(); //In units of R_sun
	void assignTms(); //In units of seconds
	void assignLuminosity(); //In erg s
	void assignTemperature(); //In units of Kelvin
	void assignIonPhotonsLum(double lowlimE, double uplimE);
	void assignMassRemnant();
	void assignEnergy();

	double blackBody(double E);
	double blackBodyPhoton(double E);
	GasProperties getGasProperties();
	int HMXBcounted{0}, Compactcounted{0};


public:
	static int s_ActiveStars;
	Star(double x = 0.0, double y = 0.0, double z = 0.0, double mass = 1.0, bool active = true)
		: Particle(x, y, z), m_active{active}, m_mass{mass} { 
			if(m_mass >= 0.5)
				alphabeta();
			if(m_mass >= 8){
				assignProperties();
				assignIonPhotonsLum(13.6, 300);
			}
			//
			assignTms();
			// assignLuminosity();
			assignMassRemnant();
			m_Lacc = 0.0;
			m_type = TYPE_STAR;
			s_ActiveStars++;
			m_gastype = GASTYPE_CORONALGAS;
		}
	void Active(bool active);
	bool SN(){ return m_SN; }
	void resetSN(){ m_SN = false; }
	double getMass(){ return m_mass; }
	double getMassRemnant(){ return m_massRemnant; }
	Star& operator= (bool active){ Active(active); }
	TYPE type(){ return m_type; }
	bool Active(){ return m_active; }
	double radiationRadius(){ return m_radiationRadius; }
	double setRadiationRadius(double alltime){ m_radiationRadius = alltime; }
	double luminosity(){ return m_L; }
	double ionluminosity() { return m_Lion; }
	double accluminosity() { return m_Lacc; }
	double starLifeTime() { return m_tms; }
	void updateWDluminosity();
	double WDluminosity(){ return m_LWD; }
	double HMXBluminosity(){ return m_LHMXB; }
	double HMXBIonluminosity(){ return m_LIonHMXB; }
	std::vector<double> HMXBionN(){ return m_NHMXB; }
	double ionN() { return m_Nion; }
	double temperature() { return m_T; }
	double HMXBtemperature() { return m_THMXB; }
	void assignHMXBproperties();
	void assignWDproperties();
	double HMXBLum(double E);
	double HMXBPhotons(double E);
	double HMXBLifeTime() { return m_thmxb; }
	void typechange();
	void accretionLuminosity(GasProperties gasproperties);
	void accretionLumNeutral();
	void reduceLifeWith(double minus);
	void setCompanion(double massCompanion){ m_massCompanion = massCompanion; }
	double massCompanion(){ return m_massCompanion; }
	double& photons(int i){ return m_photons[i]; }
	void assignPhotons(std::vector<double> E){
		m_photons.resize(E.size());
		for(int j = 0; j < E.size(); j++){
			m_photons[j] = 1.0;
		}
	}
	double integrate(double lowlim, double uplim, double (Star::*function)(double));
	double findmax(double x1, double x2, double (Star::*function)(double));
};

//Particle type, exists out of different elements.
class Gas: public Particle
{
public:
	struct AtomIon {
		double H = 0.7;
		double HI = 0.0;
		double He = 0.3;
		double HeI = 0.0;
		double HeII = 0.0;
		double e = 0.0;
		double X = 1.0;
	};
private:
	double m_U;
	double m_mass;
	double m_density;
	bool m_ionized;

	AtomIon m_numberDensity;
	AtomIon m_speciesFrac;
	AtomIon m_speciesNum;
	std::vector<double> m_distanceStar{0};

public:
	static Atoms s_mass;
	static Atoms s_Eth;
	static int s_NumIonized;
	Gas(double x = 0.0, double y = 0.0, double z = 0.0,  double ionized = 0.0, double U = 0.0)
		: Particle(x, y, z), m_ionized{ionized}, m_U{U} {}
	auto recombinations(Atoms alphaR, double timestep) -> void;
	auto assignMassDensity(double mass, double density) -> void;
	void changeNum(double H, double He, double HeI, double E = 0.0);
	auto measureDistanceToStars(std::vector<Star> star) -> void;
	auto distanceStar(int i = 0) -> double 				{ return m_distanceStar[i]; }
	auto ionized() -> bool&								{ return m_ionized; }
	auto speciesNum() -> AtomIon&						{ return m_speciesNum; }
	auto addEnergy(double E) -> void					{ m_U += E; }
	auto Energy() -> double								{ return m_U; }
	auto denmass() -> double 							{ return m_density/m_mass; }
	auto mass() -> double 								{ return m_mass; }
	auto temperature() -> double						{ return m_U/m_speciesNum.X*80.8; }
	auto speciesFrac() -> AtomIon 						{ return m_speciesFrac; }
	friend std::ostream& operator<<(std::ostream& out, const Gas::AtomIon &atomion);

};

class GasProfile
{
public:
	void readCoolingFunction();
	void setProfiles();
	void writeProperties();
	void readCNFW(double mass);
private:
	std::array<double, 701> m_TLam, m_Lambda;
	std::vector<double> m_r, m_n, m_T, m_cooling, m_coolIntegrated, m_coolingTime, m_coolingTimeff, m_coolingtest, m_coolingff;
	double m_particleMass{constants::hydrogen_g};
public:
	GasProfile(){
		m_r.resize((*All).LogBins);
		m_n.resize((*All).LogBins);
		m_T.resize((*All).LogBins);
		m_cooling.resize((*All).LogBins);
		m_coolIntegrated.resize((*All).LogBins);
		m_coolingTime.resize((*All).LogBins);
		m_coolingTimeff.resize((*All).LogBins);
		m_coolingtest.resize((*All).LogBins);
		m_coolingff.resize((*All).LogBins);
		readCNFW((*All).TotalMass);
		readCoolingFunction();
		setProfiles();
		setCooling();
		writeProperties();
	}
	void HernquistProfile();
	void NFWProfile();
	void ConstantProfile();
	void IsothermalProfile();
	void setCooling();
	void logRadius(double maxRad);
};

//Assembly of all types of particles.
class AllParticles
{
public:
	struct Recomb {
		double a;
		double b;
		double T0;
		double T1;
	};

private:
	std::vector<double> m_maxDistance{0};
	std::vector<int> m_iobjects{0}; //m_starindices of hmxbs, bhs and nss
	std::vector<Gas>  m_gas{0};
	std::vector<Star> m_star{0};
	std::vector<double> m_E, m_Nphoton, m_NphotonHMXB;
	std::vector<Atoms> m_crossSection{0};
	std::array<double, 701> m_TLam, m_Lambda;
	Atoms m_totalIonized;
	Recomb m_H, m_He, m_HeI;
	double m_timestep;
	double m_totalStarMass{0}, m_cumulativeMass{0}, m_newMass{0};
	double m_escapefraction{1};
	double m_NphotonTotal{0};
	double m_velR{0};
	int m_numHMXB, m_numStar, m_numLMXB, m_numBH, m_numNS, m_numWD, m_numSN;
	double m_lumStars, m_lumHMXBs, m_lumStars8, m_lumStarsIon, m_lumHMXBsIon, m_NionStars, m_NionHMXBs, m_lumBH, m_lumNS, m_lumWD;
	double m_r; //In case of disk
	double m_smallStars{0}, m_smallStarsMass{0};
	double m_massStar{0}, m_massHMXB{0}, m_massLMXB{0}, m_massBH{0}, m_massNS{0}, m_massWD{0}, m_massSN{0};

public:
	void readPhotoCrossSections();
	void readElectronCrossSection();
	void readRecombRates();
	void readCoolingFunction();
	//void assignGasProperties();
	void assignEnergy();
	void writeStarSpectra(std::string fileName);
	void writeStarPhotons(std::string fileName);
	void writeStarPhotonsPerGasType();
	void writeMasses();
	static double s_allTime;
public:
	AllParticles(int numGas = 0, int numStar = 0, double mass = 10.0, double density = 2.387)
	{
		assignEnergy();
		m_gas.resize(numGas);
		m_star.resize(numStar);
		m_maxDistance.resize(numStar);
		m_crossSection.resize(m_E.size());
		readPhotoCrossSections();
		readElectronCrossSection();
		//assignGasProperties();
		for(int i = 0; i < numGas; i++){
			m_gas[i].assignMassDensity(mass, density);
		}
		if(!(*All).PopulationOn || ThisTask==0){
			writeStarSpectra("MassFile_" + std::to_string(ThisTask) + ".txt");
			writeStarPhotons("StarPhotons_" + std::to_string(ThisTask) + ".txt");
		}
	}

//Stellar and gas properties
	void giveGasRandomPositions(double rad);
	void giveStarRandomPositions(double rad);
	void measureDistanceGasToStars();
	void computeDensity(Gas &gas, int starID);
	void measureDensityGasToStars();
	void calculatePhotonNumber();
	auto orderDistanceToStar(int starID) -> std::vector<int>;
	auto orderPhotons(int starID) -> std::vector<int>;
	void setRadiationRadius(){
		for(int i = 0; i < Particle::s_NumStar; i++){
			m_star[i].setRadiationRadius(s_allTime);
			m_star[i].assignPhotons(m_E);
		}
	}
	void populationAtTimeT();
	void starsActive();
	void luminosityStarHMXB();
	
//Ionization
	void radiation();
	auto recombinationCoefficient(double T) -> Atoms;
	auto crossSection(int i) -> Atoms 	{ return m_crossSection[i]; }	
	auto chanceOfIonizing(Gas& gas, Star& star, int photonIndex) -> void;
	auto findGoodDividingFactor(Gas gas, Star star, int photonIndex) -> double;
	auto setIonizationAll() -> void;
	auto printStars() -> void;
	auto printGas() -> void;
	auto writePositionsToTxt(int t) -> void;
	auto getTemperature() -> double;

	void makeStars(double mass, double &M);
	void constantSFR();
	void starBurstExp();
	void starBurst();
	void diskSFR();
	void logNorm();
	void writeStarDataToTxt();

	auto countTypePerMassBin(std::vector<double> mass, Star::TYPE type) -> std::vector<double>;
	void luminosityStarHMXB(std::string fileName);
	void photonsStar(std::vector<double>& Nphoton, std::string fileName);
	void photonsHMXB(std::vector<double>& Nphoton);

	auto equilfunction(double leftfac, double r, double nH, std::vector<double> N_gamma, double x) -> double;
	auto findroot(double leftfac, double r, double nH, std::vector<double> N_gamma, double x1, double x2) -> double;
	auto getEquilibriumHfrac(GasProperties gasproperties, std::vector<double> N_gamma, double star) -> double;
	auto lostPhotons(GasProperties gasprop, std::vector<double> Nphotons, double stars) -> std::vector<double>;
	void getEscapeFraction(GasProperties gasproperties);
	void getEscFracPerGasBubble();
	void getEscFracPerStar();
	double equilStromgren(double density, double alpha, double x);
	double findStromgrenRad(double density, double alpha, double x1, double x2);
	double findStromgrenRad(double density, double alpha);
	double getStromgrenRadiusPerGasType(GASTYPE gastype);
	double countTotalStellarMass();
	double countTotalStarMass();
	void countMassPerGroup();

//Running everything
	void runAnnuli(int radbin);
	void runPopulation();
	void addPopFiles();
	auto everythingAtTimeT(double timestep = 1.0) -> void;
};

#endif