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

double randomNumber(double min, double max){
	std::random_device rd;
	std::mt19937 mersenne(rd());
	std::uniform_real_distribution<double> uni(min, max);
	return uni(mersenne);
}
std::array<int, 1000> thousandRandomIndex(double max){
	std::array<int, 1000> partIndex;
	// for(int i = 0; i < partIndex.size(); i++){
	// 	int randvalue{static_cast<int>(randomNumber(0, max))};
	// 	bool samevalue{false};
	// 	do{
	// 		samevalue = false;
	// 		for(int j = 0; j < i; j++){
	// 			if(randvalue == partIndex[j]){
	// 				samevalue = true;
	// 			}
	// 		}
	// 		randvalue = static_cast<int>(randomNumber(0, max));
	// 	}
	// 	while( samevalue );
	// 	partIndex[i] = randvalue;
	// }
	for(int i = 0; i < partIndex.size(); i++){
		int randvalue{static_cast<int>(randomNumber(0, max))};
		partIndex[i] = randvalue;
	}
	return partIndex;
}
std::vector<double> multiplyDoubleVector(std::vector<double> vec, double val){
	std::vector<double> vecout(vec.size());
	for(int i = 0; i < vec.size(); i++){
		vecout[i] = vec[i]*val;
	}
	return vecout;
}
double sumDoubleVector(std::vector<double> vec){
	double sum;
	for(int i = 0; i < vec.size(); i++){
		sum += vec[i];
	}
	return sum;
}
//Position
Position& Position::operator= (const Position &pos){
	if (this == &pos)
		return *this;
	m_Pos = pos.m_Pos;
	return *this;
}
double& Position::operator[](int row){
	if (row <= 2)
		return m_Pos[row];
	else
		return m_Pos[2];
}
std::ostream& operator<<(std::ostream& out, const Position &pos){
	out << pos.m_Pos[0] << "\t" << pos.m_Pos[1] << "\t" << pos.m_Pos[2] << "\t";
	return out;
}
std::ostream& operator<<(std::ostream& out, const Gas::AtomIon &atomion){
	out << atomion.H << "\t" << atomion.HI << "\t" << atomion.He << 
	"\t" << atomion.HeI << "\t" << atomion.HeII << "\t";
	return out;
}
Position& Position::operator= (std::array<double, 3> pos){
	m_Pos = pos;
}

//Particle
double Particle::distanceBetween(Particle x){
	return sqrt((x.m_Pos[0] - m_Pos[0])*(x.m_Pos[0] - m_Pos[0])
		+ (x.m_Pos[1] - m_Pos[1])*(x.m_Pos[1] - m_Pos[1])
		+ (x.m_Pos[2] - m_Pos[2])*(x.m_Pos[2] - m_Pos[2]));
}
double distanceBetween(Particle x, Particle y){
	return sqrt((x.m_Pos[0] - y.m_Pos[0])*(x.m_Pos[0] - y.m_Pos[0])
		+ (x.m_Pos[1] - y.m_Pos[1])*(x.m_Pos[1] - y.m_Pos[1])
		+ (x.m_Pos[2] - y.m_Pos[2])*(x.m_Pos[2] - y.m_Pos[2]));
}

void AllParticles::giveGasRandomPositions(double rad){
	Particle tijdelijk(0, 0, 0);
	for(int i = 0; i < Particle::s_NumGas; i++){
		if(i == 0){
			m_gas[i].Pos() = std::array<double, 3>{rad, 0, 0};
			continue;
		}
		do{
			m_gas[i].Pos() = std::array<double, 3>{randomNumber(-rad, rad), randomNumber(-rad, rad), randomNumber(-rad, rad)};
		}
		while( distanceBetween(tijdelijk, m_gas[i]) > rad );
	}
}
void AllParticles::giveStarRandomPositions(double rad){
	for(int i = 0; i < Particle::s_NumStar; i++){
		if(i == 0)
			m_star[i].Pos() = std::array<double, 3>{0, 0, 0};
		else{
			do{
				m_star[i].Pos() = std::array<double, 3>{randomNumber(-rad, rad), randomNumber(-rad, rad), randomNumber(-rad, rad)};
				std::cout << m_star[i].Pos() << "\n";
			}
			while( distanceBetween(m_star[0], m_star[i]) > rad );
		}
	}
}