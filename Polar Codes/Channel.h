#ifndef CHANNEL_H
#define CHANNEL_H

#include"Modulation.h"
#include<vector>

class Channel
{
public:
	Channel(const std::vector<double>& info, double _u=0, double _sigma2 = 1);
	Channel(const Modulation& modu, double _u=0, double _sigma2 = 1);
	Channel();
	~Channel();

	void add_gauss();

public:
	std::vector<double> infoIn;
	std::vector<double> infoOut;
	std::size_t infoLen;
	double u;
	double sigma;
	double sigma2;
};

#endif