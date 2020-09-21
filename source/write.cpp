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

void writeVector(std::vector<std::vector<double>> outputvector, std::string fileName, char x){
	std::ofstream myfile;
	if(x == 'w')
		myfile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	else if(x == 'a')
		myfile.open(fileName, std::ios_base::app);
	for(int j = 0; j < outputvector[0].size(); j++){
		for(int i = 0; i < outputvector.size(); i++){
			myfile << outputvector[i][j];
			if(i == outputvector.size()-1)
				myfile << "\n";
			else
				myfile << "\t";
		}
	}
	myfile.close();
}
void AllParticles::writeStarDataToTxt(){
    std::ofstream myfile;
    myfile.open("Stars.txt", std::ofstream::out | std::ofstream::trunc);
    myfile << "x\ty\tz\tm\n";
    for(int i = 0; i < m_star.size(); i++){
        myfile << m_star[i].Pos() << m_star[i].getMass();
        myfile << "\n";
    }
    myfile.close();
}

void AllParticles::writePositionsToTxt(int t){
	std::ofstream myfile, myfilestar;
	myfile.open("Positions/ParticlePositionsIon_" + std::to_string(t) + ".txt", std::ofstream::out | std::ofstream::trunc);
	myfile << "x\ty\tz\tH\tHI\tHe\tHeI\tHeII\tT\n";
	for(int i = 0; i < m_gas.size(); i++){
		myfile << m_gas[i].Pos() << m_gas[i].speciesFrac() << m_gas[i].temperature();
		myfile << "\n";
	}
	myfile.close();
	myfilestar.open("Positions/StarPositions_" + std::to_string(t) + ".txt", std::ofstream::out | std::ofstream::trunc);
	myfilestar << "\txstar\tystar\tzstar\n";
	for(int i = 0; i < m_star.size(); i++){
		myfilestar << m_star[i].Pos();
		myfilestar << "\n";
	}
	myfilestar.close();
}

void AllParticles::printStars(){
	for(int i = 0; i < m_star.size(); i++){
		std::cout << "StarID: " << m_star[i].ID() << "\t" << m_star[i].Pos() <<  "\t Active: " << m_star[i].Active() << "\n";
	}
}
void AllParticles::printGas(){
	for(int i = 0; i < m_gas.size(); i++){
		std::cout << "GasID: " << m_gas[i].ID() << "\t" << m_gas[i].Pos();
		std::cout << "\t Ionized: " << m_gas[i].ionized() << " ";
		std::cout << "\n";					
	}
}
void AllParticles::writeStarSpectra(std::string fileName){ //Mass, luminosity, lumion, Nion, acclumMol, acclumDiff, acclumInt, acclumCor
	double deltaMass{((*All).MaximumMass - (*All).MinimumMass)/(*All).MassBins};
	double mass{(*All).MinimumMass};
	std::ofstream massfile;
	massfile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	while(mass < (*All).MaximumMass){
		Star startijdelijk(0, 0, 0, mass);

		startijdelijk.assignProperties();

		startijdelijk.assignIonPhotonsLum(13.6, 300);

		massfile << mass << "\t" << startijdelijk.luminosity() << "\t" << startijdelijk.ionluminosity() << "\t" 
		<< startijdelijk.ionN();

		if(mass >= 8.0){
			startijdelijk.accretionLuminosity(getGasProperties(GASTYPE_MOLECULARCLOUD));
			massfile << "\t" << startijdelijk.accluminosity();
			startijdelijk.accretionLuminosity(getGasProperties(GASTYPE_DIFFUSECLOUD));
			massfile << "\t" << startijdelijk.accluminosity();
			startijdelijk.accretionLuminosity(getGasProperties(GASTYPE_INTERCLOUDMEDIUM));
			massfile << "\t" << startijdelijk.accluminosity();
			startijdelijk.accretionLuminosity(getGasProperties(GASTYPE_CORONALGAS));
			massfile << "\t" << startijdelijk.accluminosity();
		}
		else{
			massfile << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
		}
		massfile << "\n";
        mass += deltaMass;
	}
	massfile.close();
}
void AllParticles::writeStarPhotons(std::string fileName){
	double deltaMass{((*All).MaximumMass - (*All).MinimumMass)/(*All).MassBins};
	double mass{(*All).MinimumMass}, energy{m_E[0]}, energynew{0};
	std::ofstream starphotons;
	starphotons.open(fileName, std::ofstream::out | std::ofstream::trunc);
	while(mass < (*All).MaximumMass){
		energy = m_E[0];
		Star startijdelijk(0, 0, 0, mass);
		startijdelijk.assignProperties();
		starphotons << mass << "\t";
		for(int i = 0; i < m_E.size(); i++){
			if(i==(m_E.size()-1)){
				if(mass< 0.5)
					starphotons << 0.0;
				else{
					startijdelijk.assignIonPhotonsLum(energy, m_E[i]);
					starphotons << startijdelijk.ionN();
				}
			}
			else{
				energynew = m_E[i+1] - (m_E[i+1]-energy)/2.0;
				if( mass < 0.5){
					starphotons << 0.0 << "\t";
					energy = energynew;
				}
				else{
					startijdelijk.assignIonPhotonsLum(energy, energynew);
					energy = energynew;
					starphotons << startijdelijk.ionN() << "\t";
				}
			}
		}
		starphotons << "\n";
		mass += deltaMass;
	}
	starphotons.close();
}
void AllParticles::writeMasses(){
	std::ofstream masses;
	masses.open("AllHMXBMasses.txt", std::ofstream::out | std::ofstream::trunc);
	masses << "Mass\tRemnantMass\tCompanionMass\tTemperature\tLuminosity\tIonLuminosity\n";
	for(int i = 0; i < m_star.size(); i++){
		if (m_star[i].type() == Star::TYPE_HMXB){
			masses << m_star[i].getMass() << "\t" << m_star[i].getMassRemnant() << "\t" << m_star[i].massCompanion() << "\t" <<
			m_star[i].temperature() << "\t" << m_star[i].HMXBluminosity() << "\t" << m_star[i].HMXBIonluminosity()<< "\n";
		}
	}
	masses.close();
}

void AllParticles::writeStarPhotonsPerGasType(){
	double deltaMass{((*All).MaximumMass - (*All).MinimumMass)/(*All).MassBins};
	double mass{(*All).MinimumMass}, energy{m_E[0]}, energynew{0};
	std::ofstream starphotM, starphotD, starphotI, starphotC;
	starphotM.open("StarPhotonsM" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
	starphotD.open("StarPhotonsD" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
	starphotI.open("StarPhotonsI" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
	starphotC.open("StarPhotonsC" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);

	while(mass < (*All).MaximumMass){
		energy = m_E[0];
		Star startijdelijk(0, 0, 0, mass);
		startijdelijk.assignProperties();
		std::vector<double> starphot(m_E.size()), starphotlost(m_E.size());
		starphotM << mass;
		starphotD << mass;
		starphotI << mass;
		starphotC << mass;
		for(int i = 0; i < m_E.size(); i++){
			if(i==(m_E.size()-1)){
				if( mass < 0.5)
					starphot[i] = 0.0;
				else{
					startijdelijk.assignIonPhotonsLum(energy, m_E[i]);
					starphot[i] = startijdelijk.ionN();
				}
			}
			else{
				energynew = m_E[i+1] - (m_E[i+1]-energy)/2.0;
				if( mass < 0.5){
					starphot[i] = 0.0;
					energy = energynew;
				}
				else{
					startijdelijk.assignIonPhotonsLum(energy, energynew);
					energy = energynew;
					starphot[i] = startijdelijk.ionN();
				}
			}
		}
		starphotlost = lostPhotons(getGasProperties(GASTYPE_MOLECULARCLOUD), starphot, m_numStar);
		for(int i = 0; i < m_E.size(); i++)
			starphotM << "\t" << starphotlost[i];
		starphotlost = lostPhotons(getGasProperties(GASTYPE_DIFFUSECLOUD), starphot, m_numStar);
		for(int i = 0; i < m_E.size(); i++)
			starphotD << "\t" << starphotlost[i];
		starphotlost = lostPhotons(getGasProperties(GASTYPE_INTERCLOUDMEDIUM), starphot, m_numStar);
		for(int i = 0; i < m_E.size(); i++)
			starphotI << "\t" << starphotlost[i];
		starphotlost = lostPhotons(getGasProperties(GASTYPE_CORONALGAS), starphot, m_numStar);
		for(int i = 0; i < m_E.size(); i++)
			starphotC << "\t" << starphotlost[i];	
		starphotM << "\n";
		starphotD << "\n";
		starphotI << "\n";
		starphotC << "\n";
		mass += deltaMass;
	}
	starphotM.close();
	starphotD.close();
	starphotI.close();
	starphotC.close();
}

void AllParticles::addPopFiles(){
    std::ofstream numfile, lumfile, massfile;
    std::ofstream starmasses;
    std::ofstream hmxb("HMXB.txt");
    std::ofstream bh("BH.txt");
    std::ofstream ns("NS.txt");
    int length{static_cast<int>((*All).TotalTime/(*All).TimeStep)};

    double *Time = new double[length]();
    double *NewMass = new double[length]();
    double *Mass = new double[length]();
    double *StarMass = new double[length]();
    double *HMXBMass = new double[length]();
    double *LMXBMass = new double[length]();
    double *BHMass = new double[length]();
    double *WDMass = new double[length]();
    double *SNMass = new double[length]();
    double *NSMass = new double[length]();
    double *Stars = new double[length]();
    double *HMXBs = new double[length]();
    double *LMXBs = new double[length]();
    double *BH = new double[length]();
    double *NS = new double[length]();
    double *WD = new double[length]();
    double *SN = new double[length]();
    double *LumStars = new double[length]();
    double *LumStars8 = new double[length]();
    double *LumStarsIon = new double[length]();
    double *NionStars = new double[length]();
    double *LumHMXBs = new double[length]();
    double *LumHMXBsIon = new double[length]();
    double *NionHMXBs = new double[length]();
    double *LumBH = new double[length]();
    double *LumNS = new double[length]();
    double *LumWD = new double[length]();

    //double *Masses = new double[length]();
    //double *Number = new double[length]();

    for(int i = 0; i < NTask; i++){
        std::ifstream f2("HMXB_" + std::to_string(i) + ".txt");
        std::ifstream f3("BH_" + std::to_string(i) + ".txt");
        std::ifstream f4("NS_" + std::to_string(i) + ".txt");
        std::string weg;
        if(i > 0){
            getline(f2, weg);
            getline(f3, weg);
            getline(f4, weg);
        }

        hmxb << f2.rdbuf();
        bh << f3.rdbuf();
        ns << f4.rdbuf();
        remove(("HMXB_" + std::to_string(i) + ".txt").c_str());
        remove(("BH_" + std::to_string(i) + ".txt").c_str());
        remove(("NS_" + std::to_string(i) + ".txt").c_str());
    }
    //Adding up the StarData files
    numfile.open("NumData.txt", std::ofstream::out | std::ofstream::trunc);
    lumfile.open("LumData.txt", std::ofstream::out | std::ofstream::trunc);
    massfile.open("MassData.txt", std::ofstream::out | std::ofstream::trunc);
    //starmasses.open("StarMass.txt", std::ofstream::out | std::ofstream::trunc);
    for(int i = 0; i < NTask; i++){
        double tijdelijk;
        std::ifstream infile1("NumData_" + std::to_string(i) + ".txt");
        for(int j = 0; j < length; j++){
            std::string strInput1;
            if(j == 0){
                getline(infile1, strInput1);
                if(i == 0)
                    numfile << strInput1 << "\n";
            }
            else{
                infile1 >> Time[j-1];
                infile1 >> tijdelijk;
                Stars[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                HMXBs[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                LMXBs[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                BH[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                NS[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                WD[j-1] += tijdelijk;
                infile1 >> tijdelijk;
                SN[j-1] += tijdelijk;
            }
        }
        remove(("NumData_" + std::to_string(i) + ".txt").c_str());

        std::ifstream infile2("LumData_" + std::to_string(i) + ".txt");
        for(int j = 0; j < length; j++){
            std::string strInput2;
            if(j == 0){
                getline(infile2, strInput2);
                if(i == 0)
                    lumfile << strInput2 << "\n";
            }
            else{
                infile2 >> Time[j-1];
                infile2 >> tijdelijk;
                LumStars[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumStars8[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumHMXBs[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumBH[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumNS[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumWD[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumStarsIon[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                NionStars[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                LumHMXBsIon[j-1] += tijdelijk;
                infile2 >> tijdelijk;
                NionHMXBs[j-1] += tijdelijk;
            }
        }
        remove(("LumData_" + std::to_string(i) + ".txt").c_str());

        std::ifstream infile3("MassData_" + std::to_string(i) + ".txt");
        for(int j = 0; j < length; j++){
            std::string strInput3;
            if(j == 0){
                getline(infile3, strInput3);
                if(i == 0)
                    massfile << strInput3 << "\n";
            }
            else{
                infile3 >> Time[j-1];
                infile3 >> tijdelijk;
                NewMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                Mass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                StarMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                HMXBMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                LMXBMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                BHMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                NSMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                WDMass[j-1] += tijdelijk;
                infile3 >> tijdelijk;
                SNMass[j-1] += tijdelijk;
            }
        }
        remove(("MassData_" + std::to_string(i) + ".txt").c_str());
    }

    /*
    for(int i = 0; i < NTask; i++){
        double tijdelijk;
        std::ifstream massin("AllStarMasses_" + std::to_string(i) + ".txt");
        for(int j = 0; j < length; j++){
            std::string strInput2;
            if(j == 0){
                getline(massin, strInput2);
                if(i == 0)
                    starmasses << strInput2 << "\n";
            }
            else{
                massin >> Masses[j-1];
                massin >> tijdelijk;
                Number[j-1] += tijdelijk;
            }
        }
        remove(("AllStarMasses_" + std::to_string(i) + ".txt").c_str());
    }
    */

    for(int i = 0; i < length - 1; i++){
    	numfile << Time[i] << "\t" << Stars[i] << "\t" << HMXBs[i] << "\t" << LMXBs[i] << "\t"
    		<< BH[i] << "\t" << NS[i] << "\t" << WD[i] << "\t" << SN[i] << "\n";
    }
    numfile.close();

    for(int i = 0; i < length - 1; i++){
    	lumfile << Time[i] << "\t" << LumStars[i] << "\t" << LumStars8[i] << "\t" << LumHMXBs[i] << "\t" << LumBH[i]  << "\t"
    		<< LumNS[i] << "\t" << LumWD[i] << "\t" << LumStarsIon[i] << "\t" << NionStars[i] << "\t" << LumHMXBsIon[i] << "\t"
    		<< NionHMXBs[i] << "\n";
    }
    lumfile.close();

    for(int i = 0; i < length - 1; i++){
    	massfile << Time[i] << "\t" << NewMass[i] << "\t" << Mass[i] << "\t" << StarMass[i] << "\t" << HMXBMass[i] << "\t"
    		<< LMXBMass[i] << "\t" << BHMass[i] << "\t" << NSMass[i] << "\t" << WDMass[i] << "\t" << SNMass[i] << "\n";
    }
    massfile.close();

    delete[] Time;
    delete[] NewMass;
    delete[] Mass;
    delete[] StarMass;
    delete[] Stars;
    delete[] HMXBs;
    delete[] LMXBs;
    delete[] BH;
    delete[] NS;
    delete[] WD;
    delete[] SN;
    delete[] HMXBMass;
	delete[] LMXBMass;
	delete[] BHMass;
	delete[] NSMass;
	delete[] WDMass;
	delete[] SNMass;
	delete[] NSMass;
    delete[] LumStars;
    delete[] LumStars8;
    delete[] LumStarsIon;
    delete[] NionStars;
    delete[] LumHMXBs;
    delete[] LumHMXBsIon;
    delete[] NionHMXBs;
    delete[] LumBH;
    delete[] LumNS;
    //delete[] Masses;
    //delete[] Number;
}
