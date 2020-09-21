#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <random>
#include <cstdio>
#include <time.h>
#include "proto.h"
#include "allvars.h"
#include "mpi.h"


int Particle::s_NumParticle = 0;

double kroupaIMF(double randval){
    if(randval < 0.435){
        return pow((randval/2.5487), 1.4286);
    }
    else if(randval < 0.8645){
        return pow((0.3*(7.11 - (randval-0.435)/0.14287)), -3.3333);
    }
    else{
        return pow((1.3*(1.894 + (0.8645 - randval)/0.07154)), -0.7692);
    }    
}
double salpeterIMF(double randval){
    if(randval > 0.01){
        return pow(1.995e-3/(0.99 - randval), 0.74074);
    }
    else{
        return 0;
    }
}
double AllParticles::countTotalStellarMass(){
    double mass{m_smallStarsMass};
    for(int i = 0; i < m_star.size(); i++){
        switch(m_star[i].type()){
            case Star::TYPE_WD:
                mass += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_HMXB:
                mass += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_LMXB:
                mass += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_NS:
                mass += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_BH:
                mass += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_STAR:
                mass += m_star[i].getMass();
                break;
        }
    }
    return mass;  
}

void AllParticles::countMassPerGroup(){
    m_massWD = 0;
    m_massHMXB = 0;
    m_massLMXB = 0;
    m_massNS = 0;
    m_massBH = 0;
    m_massStar = 0;

    for(int i = 0; i < m_star.size(); i++){
        switch(m_star[i].type()){
            case Star::TYPE_WD:
                m_massWD += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_HMXB:
                m_massHMXB += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_LMXB:
                m_massLMXB += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_NS:
                m_massNS += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_BH:
                m_massBH += m_star[i].getMassRemnant();
                break;
            case Star::TYPE_STAR:
                m_massStar += m_star[i].getMass();
                break;
        }
    }
}

double AllParticles::countTotalStarMass(){
    double mass{m_smallStarsMass};
    for(int i = 0; i < m_star.size(); i++){
        switch(m_star[i].type()){
            case Star::TYPE_WD:
                break;
            case Star::TYPE_HMXB:
                break;
            case Star::TYPE_LMXB:
                break;
            case Star::TYPE_NS:
                break;
            case Star::TYPE_STAR:
                mass += m_star[i].getMass();
                break;
        }
    }
    return mass;
}

void AllParticles::makeStars(double mass, double &M){
    /*
    Making stars with a total mass of 'mass' or more, depending on the 
    last formed star.
    Distribution depending on the chosen IMF.
    */
    double (*fcnPtr)(double);
    switch((*All).Imf){
        case IMF_KROUPA:
            fcnPtr = kroupaIMF;
            break;
        case IMF_SALPETER:
            fcnPtr = salpeterIMF;
            break;
        default:
            fcnPtr = kroupaIMF;
            break;
    }
    double m;
    double massCompanion{0};
    M = 0;
    int i{0};
    while(M < mass){
        m = fcnPtr(randomNumber(0.0, 0.99999));
        if( m > (*All).MinimumMass && m < (*All).MaximumMass ){
            m_star.push_back(Star(0.0, 0.0, 0.0, m));
            if( m > 8.0 ){
                massCompanion = randomNumber(0.01, m);
                if(massCompanion > m){
                    std::cout<<"Too much mass in makeStars()\n";
                    throw std::exception();
                }
                m_star.back().setCompanion(massCompanion);
            }
            M += m;
        }
    }
}

void AllParticles::constantSFR(){
    m_totalStarMass = countTotalStellarMass();
    if(m_totalStarMass >= s_allTime*(*All).Sfr)
        return;
    if((*All).TotalMass){
        if(m_totalStarMass >= (*All).TotalMass)
            return;
    }
    double M{0};
    makeStars(m_timestep*(*All).Sfr, M);
    m_cumulativeMass += M;
    m_newMass = M;
}

void AllParticles::starBurstExp(){
    /*
    SFR defined by an exponential function
    */
    double stepnew;
    stepnew = exp(-(s_allTime+m_timestep)/(*All).BurstWidth);
    m_totalStarMass = countTotalStellarMass();
    if(m_cumulativeMass >= (*All).TotalMass*(1.0-stepnew))
        return;
    double M{0};
    makeStars((*All).TotalMass*(exp(-(s_allTime)/(*All).BurstWidth) - stepnew), M);
    m_cumulativeMass += M;
    m_newMass = M;
}

void AllParticles::logNorm(){
    /*
    SFR defined by logNorm from Diemer et al. (2017)
    Using sigma and total mass as input to calculate tau, A and T0 from the paper
    */

    double stepnew, stepold;
    double tau{std::asinh((*All).BurstWidth/2.0/(*All).PeakTime)};
    double A{(*All).TotalMass/0.6};
    double T0{std::log((*All).PeakTime) + tau*tau};
    stepold = A/2.0*(1.0 - erf(-(std::log(s_allTime) - T0)/tau/sqrt(2.0)));
    stepnew = A/2.0*(1.0 - erf(-(std::log(s_allTime + m_timestep) - T0)/tau/sqrt(2.0)));
    double M{0};
    m_totalStarMass = countTotalStellarMass();
    if(m_cumulativeMass >= stepnew)
        return;
    makeStars(stepnew - stepold, M);
    m_cumulativeMass += M;
    m_newMass = M;
}

void AllParticles::diskSFR(){
    /*
    SFR that depends on disc properties:
    SFR(deltat) = deltat * M_at_r / r * v
    with M_at_r = Md/Rd * ((Rd + R1) * exp(-R1/Rd) - (Rd + R2) * exp(-R2/Rd))
    with R1 = r - 0.5*deltaR and R2 = r + 0.5*deltaR
    */

    using namespace constants;
    double Rd{(*All).DiskRadius}, deltaR{((*All).DiskRadius/(*All).DiskBins)};
    double Mtot{(*All).TotalMass/0.16};
    double Md{(*All).TotalMass};
    double MatR{Md/Rd*((Rd + m_r - 0.5*deltaR)*exp(-(m_r-0.5*deltaR)/Rd)-(Rd + m_r + 0.5*deltaR)*exp(-(m_r+0.5*deltaR)/Rd))};
    m_totalStarMass = countTotalStellarMass();
    if(m_totalStarMass >= MatR)
        return;
    double M_disk_r{Md*Msun_g/Rd*(Rd - (Rd + m_r)*exp(-m_r/Rd))};
    double r_Hern{virialRadius(Mtot*0.84)*hernquistFactor()/(1.0-hernquistFactor())};
    double M_halo_r{Mtot*0.84*m_r*m_r/r_Hern/r_Hern/pow(1.0 + m_r/r_Hern, 2)};
    m_velR = sqrt(G_cm3_gi_si2 * (M_disk_r+M_halo_r) / m_r);
    double M{0};
    makeStars(m_timestep/s_to_yr*MatR/(m_r/m_velR), M);
    m_totalStarMass += M;
    m_newMass = M;
    //std::cout << "Total mass: " << m_totalStarMass << "\n";
}

void AllParticles::starBurst(){
    double M{0};
    makeStars((*All).TotalMass, M);
    m_totalStarMass += M;
    m_newMass = M;
}

void AllParticles::starsActive(){
    /*
    Aging objects and classifying them appropriately
    */

    m_numHMXB = 0;
    m_numStar = m_smallStars;
    m_numLMXB = 0;
    m_numBH = 0;
    m_numNS = 0;
    m_numWD = 0;
    m_numSN = 0;
    m_massSN = 0;
    m_iobjects.clear();
    for(int i = 0; i < m_star.size(); i++){
        m_star[i].reduceLifeWith(m_timestep); //Aging

        if(m_star[i].SN() == true){ //Counting SN
            m_numSN++;
            m_massSN += m_star[i].getMass();
            m_star[i].resetSN();
        }
        switch(m_star[i].type()){
            case Star::TYPE_STAR:
                m_numStar++;
                break;
            case Star::TYPE_HMXB:
                m_numHMXB++;
                m_iobjects.push_back(i);
                break;
            case Star::TYPE_LMXB:
                m_numLMXB++;
                m_iobjects.push_back(i);
                break;
            case Star::TYPE_BH:
                m_numBH++;
                m_iobjects.push_back(i);
                break;
            case Star::TYPE_NS:
                m_numNS++;
                m_iobjects.push_back(i);
                break;
            case Star::TYPE_WD:
                m_numWD++;
                break;
            default:
                break;
        }
    }
}

std::vector<double> AllParticles::countTypePerMassBin(std::vector<double> mass, Star::TYPE type){
    /*
    Counting objects of certain type per mass bin
    */

    std::vector<double> typePerMassBin;
    typePerMassBin.resize(mass.size());
    for(int j = 0; j < m_star.size(); j++){
        if(m_star[j].type() == type){
            if(m_star[j].getMass() > mass[mass.size() - 1])
                typePerMassBin[mass.size()-1] += 1;
            else{
                for(int i = 0; i < mass.size(); i++){
                    if(m_star[j].getMass() - mass[i] < 0){
                        typePerMassBin[i - 1] += 1;
                        break;
                    }
                }
            }

        }
    }
    return typePerMassBin;
}

void AllParticles::photonsStar(std::vector<double>& Nphoton, std::string fileName){
    /*
    Counting photons from stars using table
    */

    std::ifstream infile(fileName);
    std::vector<double> mass;
    std::vector<std::vector<double>> photons(m_E.size(), std::vector<double>(0));
    double value;
    m_NphotonTotal = 0;
    m_escapefraction = 1.0;
    while( infile >> value ) {
        mass.push_back(value);
        for(int i = 0; i < m_E.size(); i++){
            infile >> value;
            photons[i].push_back(value);
        }
    }
    if(Nphoton.size() != m_E.size())
        Nphoton.resize(m_E.size());
    std::vector<double> starsPerMassBin{countTypePerMassBin(mass, Star::TYPE_STAR)};
    for(int i = 0; i < m_E.size(); i++){
        Nphoton[i] = 0;
    }
    for(int j = 0; j < mass.size(); j++){
        for(int i = 0; i < m_E.size(); i++){
            Nphoton[i] += photons[i][j]*starsPerMassBin[j];
        }
    }
    for(int i = 0; i < m_E.size(); i++){
        m_NphotonTotal += m_Nphoton[i];
    }
}
void AllParticles::photonsHMXB(std::vector<double>& Nphoton){
    /*
    Counting photons from HMXBs for the given range of energies
    */

    if(Nphoton.size() != m_E.size())
        Nphoton.resize(m_E.size());
    for(int i = 0; i < m_E.size(); i++){
        Nphoton[i] = 0;
    }
    for(int i = 0; i < m_E.size(); i++){
        for(int j = 0; j < m_star.size(); j++){
            if(m_star[j].type() == Star::TYPE_HMXB){
                Nphoton[i] += m_star[j].HMXBionN()[i];
            }
        }
    }
}

void AllParticles::luminosityStarHMXB(std::string fileName){
    /*
    Counting luminosities of all stars and HMXBs in the timestep, using the input file
    */

    std::ifstream infile(fileName);
    std::vector<double> mass, lum, lumion, nion, laccM, laccD, laccI, laccC;
    double value, timeFrac;
    while( infile >> value) {
        mass.push_back(value);
        infile >> value;
        lum.push_back(value);
        infile >> value;
        lumion.push_back(value);
        infile >> value;
        nion.push_back(value);
        infile >> value;
        laccM.push_back(value);
        infile >> value;
        laccD.push_back(value);
        infile >> value;
        laccI.push_back(value);
        infile >> value;
        laccC.push_back(value);
    }
    std::vector<double> starsPerMassBin{countTypePerMassBin(mass, Star::TYPE_STAR)};
    std::vector<double> HMXBsPerMassBin{countTypePerMassBin(mass, Star::TYPE_HMXB)};
    std::vector<double> NSPerMassBin{countTypePerMassBin(mass, Star::TYPE_NS)};
    std::vector<double> BHPerMassBin{countTypePerMassBin(mass, Star::TYPE_BH)};
    m_lumStars = 0;
    m_lumStars8 = 0;
    m_lumHMXBs = 0;
    m_lumStarsIon = 0;
    m_lumHMXBsIon = 0;
    m_NionStars = 0;
    m_NionHMXBs = 0;
    m_lumBH = 0;
    m_lumNS = 0;
    m_lumWD = 0;
    for(int i = 0; i < m_star.size(); i++){

        if(m_star[i].type()==Star::TYPE_HMXB){
            if(m_star[i].HMXBLifeTime() < 0)
                timeFrac = 0;
            else if(m_star[i].HMXBLifeTime() < (*All).TimeStep)
                timeFrac = m_star[i].HMXBLifeTime()/(*All).TimeStep;
            else
                timeFrac = 1.0;

            m_lumHMXBs += m_star[i].HMXBluminosity()*timeFrac;
            m_lumHMXBsIon += m_star[i].HMXBIonluminosity()*timeFrac;
            for(int j = 0; j < m_E.size(); j++){
                if( m_E[j] >= 13.6 && m_E[j] <= 1e3){
                    m_NionHMXBs += m_star[i].HMXBionN()[j]*timeFrac;
                }
            }
        }
        else if(m_star[i].type()==Star::TYPE_BH){
            m_lumBH += m_star[i].accluminosity();
        }
        else if(m_star[i].type()==Star::TYPE_NS){
            m_lumNS += m_star[i].accluminosity();
        }
        else if(m_star[i].type()==Star::TYPE_WD){
            m_lumWD += m_star[i].WDluminosity();
        }
        else if((m_star[i].type()==Star::TYPE_STAR)&&(m_star[i].getMass() >= 8)){
            if(m_star[i].starLifeTime() < 0)
                timeFrac = 0;
            else if(m_star[i].starLifeTime() < (*All).TimeStep)
                timeFrac = m_star[i].starLifeTime()/(*All).TimeStep;
            else
                timeFrac = 1.0;

            m_lumStars8 += m_star[i].luminosity()*timeFrac;
            m_lumStars += m_star[i].luminosity()*timeFrac;
            m_lumStarsIon += m_star[i].ionluminosity()*timeFrac;
            m_NionStars += m_star[i].ionN()*timeFrac;
        }
    }
    for(int i = 0; i < mass.size(); i++){
        if(mass[i] >= 8)
            continue; //Already added up
        m_lumStars += starsPerMassBin[i]*lum[i];
        m_lumStarsIon += starsPerMassBin[i]*lumion[i];
        m_NionStars += starsPerMassBin[i]*nion[i];
        // m_lumAcc += getGasProperties(GASTYPE_MOLECULARCLOUD).volume*(NSPerMassBin[i]+BHPerMassBin[i])*laccM[i];
        // m_lumAcc += getGasProperties(GASTYPE_DIFFUSECLOUD).volume*(NSPerMassBin[i]+BHPerMassBin[i])*laccD[i];
        // m_lumAcc += getGasProperties(GASTYPE_INTERCLOUDMEDIUM).volume*(NSPerMassBin[i]+BHPerMassBin[i])*laccI[i];
        // m_lumAcc += getGasProperties(GASTYPE_CORONALGAS).volume*(NSPerMassBin[i]+BHPerMassBin[i])*laccC[i];
    }
}

void AllParticles::populationAtTimeT(){
    /* 
    Calling right SFR function depending on parameter file input 
    */

    char sfr;
    void (AllParticles::*fcnPtrM)();

    if(!(*All).SFRType){
        std::cout<< "No SFR in populationAtTimeT()\n";
        throw std::exception();
    }
    switch((*All).SFRType){
        case SFR_EXPONENTIAL:
            fcnPtrM = &AllParticles::starBurstExp;
            break;
        case SFR_CONSTANT:
            fcnPtrM = &AllParticles::constantSFR;
            break;
        case SFR_STARBURST:
            fcnPtrM = &AllParticles::starBurst;
            break;
        case SFR_DISK:
            fcnPtrM = &AllParticles::diskSFR;
            break;
        case SFR_LOGNORM:
            fcnPtrM = &AllParticles::logNorm;
            break;
        default:
            fcnPtrM = &AllParticles::starBurst;
            break;
    }
    if((*All).SFRType == SFR_STARBURST && s_allTime > m_timestep);
    else
        (*this.*fcnPtrM)();
}

void AllParticles::runAnnuli(int radbin){
    /*
    Calculating stromgren radius for N annuli in a exponential disc, for different types of gas.
    */

    m_r = (radbin + 0.5)*(*All).DiskRadius/(*All).DiskBins;
    std::ofstream stromfile;
    stromfile.open("StromgrenRadius_" + std::to_string(radbin) + ".txt", std::ofstream::out | std::ofstream::trunc);
    stromfile << "Time\tAnnulusRad\tVelRad\tStars\tNion\tStromgrenRadM\tStromgrenRadD\tStromgrenRadI\tStromgrenRadC\n";
    m_timestep = (*All).TimeRes;
    readRecombRates();
    while(s_allTime < (*All).TotalTime){
        s_allTime += m_timestep;
        populationAtTimeT();
        starsActive();
        if(fmod(s_allTime,((*All).TimeStep))==0){
            luminosityStarHMXB("MassFile_" + std::to_string(ThisTask) + ".txt");
            if((*All).PhotonsOn)
                photonsStar(m_Nphoton, "StarPhotons_" + std::to_string(ThisTask) + ".txt");
            stromfile << s_allTime << "\t" << m_r << "\t" << m_velR << "\t" << m_numStar << "\t";
            if((*All).PhotonsOn)
                stromfile << sumDoubleVector(m_Nphoton) << "\t"; 
            stromfile << getStromgrenRadiusPerGasType(GASTYPE_MOLECULARCLOUD) << "\t" << getStromgrenRadiusPerGasType(GASTYPE_DIFFUSECLOUD) << 
            "\t" << getStromgrenRadiusPerGasType(GASTYPE_INTERCLOUDMEDIUM) << "\t" << getStromgrenRadiusPerGasType(GASTYPE_CORONALGAS) << "\n";
        }   
    }
    stromfile.close();
}

void AllParticles::runPopulation(){
    /*
    Depending on your choice of parameters for IMF and SFR, a population of stars, HMXB, and compact objects
    is calculated over time.
    There is an option for radiative transfer (RadiationOn), but it is not worked out completely yet.
    */

    /* Making file for stardata: mass, stellar mass, total objects, stars, HMXBs, LMXBs, etc. */
    time_t start = time(0);
    std::ofstream escapefracs;

    std::ofstream numfile, lumfile, massfile;
    numfile.open("NumData_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
    lumfile.open("LumData_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
    massfile.open("MassData_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);

    numfile << "Time\tStar\tHMXB\tLMXB\tBH\tNS\tWD\tSN\n";
    lumfile << "Time\tLumStar\tLumStar8\tLumHMXB\tLumBH\tLumNS\tLumWD\tLumStarIon\tNionStar\tLumHMXBIon\tNionHMXB\n";
    massfile << "Time\tNewMass\tMass\tStarMass\tHMXBMass\tLMXBMass\tBHMass\tNSMass\tWDMass\tSNMass\n";

    /* Making file for every hmxb, bh and ns that is born and their properties */
    std::ofstream hmxb, bh, ns;
    hmxb.open("HMXB_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
    hmxb << "Mass\tRemnantMass\tCompanionMass\tTemperature\tLuminosity\tIonLuminosity\tLifeTime\tTime\n";
    bh.open("BH_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
    bh << "Mass\tRemnantMass\tLuminosity\tTime\n";
    ns.open("NS_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
    ns << "Mass\tRemnantMass\tLuminosity\tTime\n";

    /*
    std::ofstream starmasses;
    if((*All).SFRType == SFR_STARBURST){
        starmasses.open("AllStarMasses_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
        starmasses << "Mass\tNumber\n";
    }
    */

    /* If you chose the option "PhotonsOn", a file will be written with the amount of photons for a range of energies */
    std::vector<std::ofstream> photons(m_E.size());
    if((*All).PhotonsOn){
        for(int i = 0; i < m_E.size(); i++){
            photons[i].open("PhotonData_" + std::to_string(i) + "_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
            photons[i] << "Time\tE\tStars\tHMXBs\n";
        }
    }
    /* If you chose the option "RadiationOn", a file will be created of the number of photons escaping the system */
    if((*All).RadiationOn){
        escapefracs.open("EscapeFractions_" + std::to_string(ThisTask) + ".txt", std::ofstream::out | std::ofstream::trunc);
        escapefracs << "Time\tNionStars\tNionStarsEscape\tStarsEscapeFraction\n";
    }

    /* Initializing values and reading recombination rates */
    m_timestep = (*All).TimeRes;
    double escpfr1{0}, escpfr2{0}, escpfr3{0}, escpfr4{0};
    if((*All).PhotonsOn)
        readRecombRates();
    double timepop{0}, timestar{0}, timelum{0}, timewrite{0}, timeobjectcount{0}, tijdelijk, allmass{0};
    while(s_allTime < (*All).TotalTime){
        s_allTime += m_timestep;

        /* Visual bar */
        if(ThisTask == 0 && fmod(s_allTime,((*All).TotalTime/100.0)) == 0){
            std::cout<< "[";
            int pos = 70* s_allTime/(*All).TotalTime;
            for(int i = 0; i < 70; ++i){
                if(i <= pos) std::cout << "|";
                else std::cout << " ";
            }
            std::cout << "] " << int(s_allTime/(*All).TotalTime*100.0) << " %\r";
            std::cout.flush();
        }

        /* Calculating distance gas to stars and ionization + writing to file */
        if((*All).RadiationOn){
            radiation();
        }

        /* Making stars */
        tijdelijk = difftime( time(0), start );
        populationAtTimeT();
        timepop += difftime( time(0), start ) - tijdelijk;

        /* Setting age and type of objects */
        tijdelijk = difftime( time(0), start );
        starsActive();
        timestar += difftime( time(0), start ) - tijdelijk;


        countMassPerGroup();
        allmass = m_massWD + m_massHMXB + m_massLMXB + m_massNS + m_massBH + m_massStar;

        /* Writing step in file, depending on the TimeStep from parameter file */
        if(fmod(s_allTime,((*All).TimeStep))==0){
            tijdelijk = difftime( time(0), start );
            luminosityStarHMXB("MassFile_0.txt"); //Only on node 0 this file was written
            timelum += difftime(time(0), start ) - tijdelijk;

            tijdelijk = difftime( time(0), start );

            numfile << s_allTime << "\t" << m_numStar << "\t" << m_numHMXB << "\t" << m_numLMXB << "\t" 
                << m_numBH << "\t" << m_numNS << "\t" << m_numWD << "\t" << m_numSN << "\n";
            lumfile << s_allTime << "\t" << m_lumStars << "\t" << m_lumStars8 << "\t" << m_lumHMXBs  << "\t"
                << m_lumBH << "\t" << m_lumNS << "\t" << m_lumWD << "\t" << m_lumStarsIon << "\t"
                << m_NionStars << "\t" << m_lumHMXBsIon << "\t" << m_NionHMXBs << "\n";
            massfile << s_allTime << "\t" << m_newMass << "\t" << m_totalStarMass << "\t" << allmass << "\t"
                << m_massHMXB << "\t" << m_massLMXB << "\t" << m_massBH << "\t" << m_massNS << "\t" << m_massWD << "\t"
                << m_massSN << "\n";

            timewrite += difftime(time(0), start ) - tijdelijk;
        }

        /* Counting stars and put it in starmasses file */
        /*
        if(s_allTime == m_timestep && (*All).SFRType == SFR_STARBURST){
            double deltaMass{((*All).MaximumMass - (*All).MinimumMass)/(*All).MassBins};
            double masstemp{(*All).MinimumMass};
            std::vector<double> starsmass;
            while(masstemp < (*All).MaximumMass){
                starsmass.push_back(masstemp);
                masstemp += deltaMass;
            }
            std::vector<double> starsnumber{countTypePerMassBin(starsmass, Star::TYPE_STAR)};
            for(int i = 0; i < starsmass.size(); i++){
                starmasses << starsmass[i] << "\t" << starsnumber[i] << "\n";
            }
            
        }
        */
        /* Put HMXB, BH and NS info in files */
        tijdelijk = difftime(time(0), start);
        for(int i = 0; i < m_iobjects.size(); i++){
            switch(m_star[m_iobjects[i]].type()){
                case Star::TYPE_STAR:
                    break;
                case Star::TYPE_WD:
                    break;
                case Star::TYPE_LMXB:
                    break;
                case Star::TYPE_HMXB:
                    if(m_star[m_iobjects[i]].HMXBcounted == 0){
                        hmxb << m_star[m_iobjects[i]].getMass() << "\t" << 
                        m_star[m_iobjects[i]].getMassRemnant() << "\t" << 
                        m_star[m_iobjects[i]].massCompanion() << "\t" << 
                        m_star[m_iobjects[i]].HMXBtemperature() << "\t" << 
                        m_star[m_iobjects[i]].HMXBluminosity() << "\t" << 
                        m_star[m_iobjects[i]].HMXBIonluminosity() << "\t" <<
                        m_star[m_iobjects[i]].HMXBLifeTime() << "\t" <<
                        s_allTime << "\n";
                        m_star[m_iobjects[i]].HMXBcounted = 1;
                    }
                    break;
                case Star::TYPE_BH:
                    if(m_star[m_iobjects[i]].Compactcounted == 0){
                        bh << m_star[m_iobjects[i]].getMass() << "\t" << 
                        m_star[m_iobjects[i]].getMassRemnant() << "\t" << 
                        m_star[m_iobjects[i]].accluminosity() << "\t" << 
                        s_allTime <<"\n";
                        m_star[m_iobjects[i]].Compactcounted = 1;
                    }
                    break;
                case Star::TYPE_NS:
                    if(m_star[m_iobjects[i]].Compactcounted == 0){
                        ns << m_star[m_iobjects[i]].getMass() << "\t" << m_star[m_iobjects[i]].getMassRemnant() << "\t" << 
                        m_star[m_iobjects[i]].accluminosity() << "\t" <<
                        s_allTime << "\n";
                        m_star[m_iobjects[i]].Compactcounted = 1;
                    }
                    break;
                default:
                    break;
            }
        }
        timeobjectcount += difftime(time(0), start) - tijdelijk;

        /* Write photons per energy and escapefraction per time */
        if((*All).PhotonsOn){
            photonsStar(m_Nphoton, "StarPhotons_0.txt"); //Only on node 0 this file was written
            photonsHMXB(m_NphotonHMXB);
            for(int i = 0; i < m_E.size(); i++)
                photons[i] << s_allTime << "\t" << m_E[i] << "\t" << m_Nphoton[i] << "\t" << m_NphotonHMXB[i] << "\n";
            if((*All).RadiationOn){
                //getEscapeFraction(getGasProperties(GASTYPE_MOLECULARCLOUD));
                //getEscapeFraction(getGasProperties(GASTYPE_DIFFUSECLOUD));
                getEscFracPerStar();
                escapefracs << s_allTime << "\t" << m_NionStars << "\t" << sumDoubleVector(m_Nphoton) << "\t" << m_escapefraction << "\n";
                //getEscapeFraction(getGasProperties(GASTYPE_DIFFUSECLOUD));
                //getEscapeFraction(getGasProperties(GASTYPE_MOLECULARCLOUD));
                //escapefracs << m_escapefraction << "\n";
            }
        }
    }

    /* Close files */
    /*
    if((*All).SFRType == SFR_STARBURST)
        starmasses.close();
    */
    hmxb.close();
    bh.close();
    ns.close();
    numfile.close();
    lumfile.close();
    massfile.close();
    if((*All).PhotonsOn){
        for(int i = 0; i < m_E.size(); i++){
            photons[i].close();
        }
    }
    if((*All).RadiationOn)
        escapefracs.close();
    if(ThisTask==0){
        std::cout << std::endl;
        std::cout << "Total mass: " << NTask*m_totalStarMass << "; Number of stars (approx): " << NTask*m_star.size() << '\n';
        std::cout << "Time it took to compute:\npopulationAtTimeT:\t" << timepop << "\nstarsActive:\t" << 
            timestar << "\nluminosityStarHMXB:\t" << timelum << "\nWriting file:\t" << timewrite << 
            "\nCounting and writing objects:\t" << timeobjectcount << "\nTotal time:\t" << difftime( time(0), start ) << "\n";
    }
}

void AllParticles::radiation()
{
    /* Gives the distance between every gas particle and star and stores it in a vector within Gas */
    measureDistanceGasToStars();

    /* Calculates the radius for which the star is able to reach as a function of time (speed of light) */
    setRadiationRadius();

    /* Is not right anymore! Update! */
    calculatePhotonNumber();

    /* Calculate ionization fraction per particle */ 
    setIonizationAll();

    /* Write properties of gas and stars to files on snapshot time */
    if(fmod(s_allTime,((*All).TimeStep))==0){
        writePositionsToTxt(static_cast<int>(s_allTime/m_timestep));
    }
}

void AllParticles::everythingAtTimeT(double timestep)
{
    s_allTime += timestep;
    m_timestep = timestep;
    int t{0};
    if(timestep > 0)
        t = static_cast<int>(s_allTime/timestep);

    readRecombRates();
    if((*All).RadiationOn){
        radiation();
    }
    starsActive();
    std::cout << "Photons left: " << m_star[0].photons(0) << "\n";
    std::cout << "Percentage HI: " << m_totalIonized.H*100 << "\n";
    std::cout << "Percentage HeI: " << m_totalIonized.He*100 << "\n";
    std::cout << "Percentage HeII: " << m_totalIonized.HeI*100 << "\n";
    std::cout << "Temperature gas: " << getTemperature() << "\n";
}