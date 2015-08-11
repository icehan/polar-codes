#include"InfoSource.h"
#include<random>

void
InfoSource::ShowInfo()
{
	for (std::size_t i = 0; i < infoBits.size(); ++i)
		std::cout << infoBits[i] << " ";
}

InfoSource::InfoSource(int a[], int len)
{
	infoLen = len;
	for (int i = 0; i < len; ++i)
		infoBits.push_back(a[i]);
}

InfoSource::InfoSource(uint infoLength)
{
	infoLen = infoLength;
	infoBits.resize(infoLength);

	std::random_device rd;
	std::mt19937 mt_gen(rd());
	std::uniform_int_distribution<uint> unif_dist(0, 1);

	for (std::size_t i = 0; i < infoLength; ++i)
		infoBits[i] = unif_dist(mt_gen);

#ifdef DEBUG	
	for (size_t i = 0; i < 10; i++)
		infoBits[i] = 0;
	infoBits[11] = infoBits[12] = infoBits[13] = 1;
#endif
}
