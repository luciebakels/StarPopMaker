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


Parameters *All = new Parameters;


int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	if(argc < 2 && ThisTask == 0){
		std::cout << "Where is the parameter file?\n";
		throw std::exception();
	}

	//read parameter file
	std::string paramfile;
	std::stringstream convert;
	convert << argv[1];
	convert >> paramfile;

	readParamFile(paramfile);
	if((*All).StromgrenOn){
		if(!(*All).DiskBins || !(*All).DiskRadius || (*All).SFRType != SFR_DISK){
			std::cout<< "You should check your parameter file. I need a value for DiskBins and " <<
				"DiskRadius, and SFRType has to be D. Otherwise, turn StromgrenOn off.\n";
			throw std::exception();
		}
		if(!(*All).PhotonsOn){
			std::cout<< "You should probably turn on PhotonsOn\n";
		}
		int diskBinsPerTask{static_cast<int>((*All).DiskBins/NTask) + 1};
		for(int i = 0; i < (*All).DiskBins; i++){
			int rad{i + ThisTask*diskBinsPerTask};
			if((i < diskBinsPerTask) && 
				( rad <= (*All).DiskBins)){
				AllParticles::s_allTime = 0;
				AllParticles allpart;
				allpart.runAnnuli(rad);
			}
		}
	}
	else if(((*All).PopulationOn || (*All).RadiationOn) && !(*All).StromgrenOn){
		AllParticles::s_allTime = 0;
		AllParticles allpart;
		MPI_Barrier( MPI_COMM_WORLD );  //The common files have to be written before proceeding
		if((*All).PopulationOn){
			allpart.runPopulation();
			MPI_Barrier( MPI_COMM_WORLD ); //All Files have to be completed before you can add them together
			if(ThisTask == 0)
				allpart.addPopFiles();
		}
		else if((*All).RadiationOn)
			allpart.everythingAtTimeT();
	}


	if((*All).GasProfileOn){
		std::array<double, 19> masses{1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10, 1e11, 5e11, 1e12, 5e12, 
			1e13, 5e13, 1e14, 5e14, 1e15};		
		for(int i = 0; i < masses.size(); i++){
			std::cout << "Calculating gas profile...\n";
			if(ThisTask == 0){
				(*All).TotalMass = masses[i];
				GasProfile gasprofile;
				std::cout << "Done!\n";
			}
		}
	}

	// Star star(0, 0, 0, 60);
	// star.assignHMXBproperties();
	// std::cout << star.HMXBluminosity() << "\t" << star.HMXBIonluminosity() << "\t" << star.HMXBtemperature() << "\n";
	delete[] All;
	//Star star(0, 0, 0, 1);
	//std::cout << star.temperature() << "\t" <<  star.luminosity() << "\t" << 
	//	star.ionluminosity() << "\n"; //constants::erg_to_eV << "\n";
	//writeStarDataToTxt(makeStars(totalMass, IMF_SALPETER, minMass, maxMass));
	/*
	std::cout << "Enter number of gas particles: ";
	std::cin >> numGas;
	std::cout << "Enter number of stars: ";
	std::cin >> numStar;
	std::cout << "Enter timestep: ";
	std::cin >> timestep;
	std::cout << "Enter number of elapses: ";
	std::cin >> loops;

	
	all.giveStarRandomPositions(rad);
	all.giveGasRandomPositions(rad);
	for(int i = 0; i < loops; i++){
		all.everythingAtTimeT(timestep);//+i*timestep);
	}*/

	MPI_Finalize();
	return 0;
}