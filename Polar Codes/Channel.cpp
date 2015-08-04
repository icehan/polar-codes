#include"Channel.h"
#include<random>

void
Channel::add_gauss()
{
	std::random_device rd;
	std::mt19937 mt_gen(rd());
	std::normal_distribution<double> norm_dist(u, sigma);

	for (std::size_t i = 0; i < infoLen; ++i)
		infoOut[i] = infoIn[i] + norm_dist(mt_gen);
}

Channel::Channel(const std::vector<double>& info, double _u, double _sigma2)
{
	infoIn = info;
	infoLen = info.size();
	infoOut.resize(infoLen);
	u = _u;
	sigma = sqrt(_sigma2);
	sigma2 = _sigma2;
}

Channel::Channel(const Modulation& modu, double _u, double _sigma2)
{
	infoIn = modu.infoOut;
	infoLen = modu.infoLen;
	infoOut.resize(infoLen);
	u = _u;
	sigma = sqrt(_sigma2);
	sigma2 = _sigma2;
}

Channel::Channel()
{
	infoLen = 0;
	u = 0.0;
	sigma = 1.0;
	sigma2 = 1.0;
}

Channel::~Channel()
{
}
