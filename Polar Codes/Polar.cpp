#include"Polar.h"
#include<iostream>

Polar::Polar(){}

Polar::~Polar(){}

uint* 
Polar::universal_encode(const uint* data, uint length, uint *dst)
{
	if(NULL == data || NULL == dst || length == 0)
		return NULL;

	for (int i = 0; i < length; i++)
		dst[i] = data[i];
	uint M = log(length*1.0)/log(2.0);
	for (uint i = 0; i < M; i++)
	{
		uint block_len = pow(2.0, i+1.0);
		uint jMax = pow(2.0, M - i - 1.0);
		for (uint j = 0; j < jMax; j++)
		{
			uint kLen = pow(2.0, i*1.0), sep = kLen;
			uint kBeg = j*block_len, kEnd = kBeg + kLen;
			for (uint k = kBeg; k < kEnd; k++)
				dst[k] = dst[k] ^ dst[k + sep];
		}
	}
	return dst;
}

void
Polar::set_pattern(std::size_t _K)
{
	const uint *pt = get_point_of_pattern(_K*2);
	if (pt == NULL)
	{
		std::cout << "Only pattern 4/128/512 served!\n";
		exit(EXIT_FAILURE);
	}

	//set_pattern
	pattern.resize(_K*2);
	for (std::size_t i = 0; i < _K*2; ++i)
		pattern[i] = pt[i];
}

const uint*
Polar::get_point_of_pattern(std::size_t _N)
{
	const uint *pt = NULL;

	if (_N == 1024)
		pt = PATTERN1024;
	else if (_N == 256)
		pt = PATTERN256;
	else if (_N == 8)
		pt = PATTERN8;
	else if (_N == 4)
		pt = PATTERN4;

	return pt;
}