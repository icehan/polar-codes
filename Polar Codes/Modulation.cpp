#include"Modulation.h"

void
Modulation::bpsk()
{
	//0 --> 1, 1 --> -1
	for(std::size_t i = 0; i < infoLen; ++i)
		infoOut[i] = 1 - 2*infoIn[i];
}

Modulation::Modulation(const std::vector<uint> &info)
{
	infoLen = info.size();
	infoIn.resize(infoLen);
	infoOut.resize(infoLen);

	for(std::size_t i = 0; i < infoLen; ++i)
		infoIn[i] = info[i];
}

Modulation::Modulation(const InfoSource& src)
{
	infoLen = src.infoLen;
	infoIn.resize(infoLen);
	infoOut.resize(infoLen);

	for(std::size_t i = 0; i < infoLen; ++i)
		infoIn[i] = src.infoBits[i];
}

Modulation::Modulation()
{
	infoLen = 0;
}

Modulation::~Modulation()
{
}
