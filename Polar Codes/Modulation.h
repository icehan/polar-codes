#ifndef MODULATION_H
#define MODULATION_H

#include"ice_type.h"
#include"InfoSource.h"
#include<vector>

class Modulation
{
public:
	Modulation(const std::vector<uint> &info);	// uint -> int
	Modulation(const InfoSource& src);
	Modulation();
	~Modulation();
	void bpsk();
public:
	std::vector<double> infoOut;	//modulated
	std::vector<int> infoIn;
	std::size_t infoLen;
};

#endif