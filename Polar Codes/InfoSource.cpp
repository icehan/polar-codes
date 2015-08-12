#include"InfoSource.h"
#include<random>

InfoSource::InfoSource(int len)
{
	infoLen = len;
	infoBits.resize(len);

	std::random_device rd;
	std::mt19937 mt_gen(rd());
	std::uniform_int_distribution<uint> unif_dist(0, 1);

	for (int i = 0; i < len; ++i)
		infoBits[i] = unif_dist(mt_gen);

#ifdef DEBUG	
	for (size_t i = 0; i < 10; i++)
		infoBits[i] = 0;
	infoBits[11] = infoBits[12] = infoBits[13] = 1;
#endif
}

InfoSource::InfoSource(uint a[], int len)
{
	infoLen = len;
	infoBits.resize(len);
	for (int i = 0; i < len; ++i)
		infoBits[i] = a[i];
}

InfoSource::InfoSource()
{
}

void
InfoSource::ShowInfo()
{
	for (int i = 0; i < infoLen; ++i)
		std::cout << infoBits[i] << " ";
}