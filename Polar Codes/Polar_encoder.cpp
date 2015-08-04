#include"Polar_encoder.h"


void
Polar_encoder::Encode()
{
	universal_encode(infoBitsAddPattern.data(), N, codeword.data());
/*
	codeword = infoBitsAddPattern;
	for (std::size_t i = 0; i < log2N; i++)
	{
		uint block_len = pow(2.0, i+1.0);
		uint jMax = pow(2.0, log2N - i - 1.0);
		for (uint j = 0; j < jMax; j++)
		{
			uint kLen = pow(2.0, i*1.0), sep = kLen;
			uint kBeg = j*block_len, kEnd = kBeg + kLen;
			for (uint k = kBeg; k < kEnd; k++)
				codeword[k] = codeword[k] ^ codeword[k + sep];	//CODEWORD must be continuing updating
		}
	}
*/
}

Polar_encoder::Polar_encoder(const std::vector<uint> &bits, double code_rate)
{
	K = bits.size();
	codeRate = code_rate;
	N = K/code_rate;
	R = N - K;
	log2N = log(N*1.0)/log(2.0);
	infoBits = bits;
	codeword.resize(N);
	infoBitsAddPattern = set_info_using_pattern(bits);
	set_pattern(bits.size());
}

Polar_encoder::Polar_encoder(const InfoSource &src, double code_rate)
{
	K = src.infoLen;
	codeRate = code_rate;
	N = K/code_rate;
	R = N - K;
	log2N = log(N*1.0)/log(2.0);
	infoBits = src.infoBits;
	codeword.resize(N);
	infoBitsAddPattern = set_info_using_pattern(infoBits);
	set_pattern(src.infoLen);
}

Polar_encoder::~Polar_encoder(){}

//private
std::vector<uint> 
Polar_encoder::set_info_using_pattern(const std::vector<uint>& info)
{
	const uint *pt = NULL;
	if (info.size() == 512)
		pt = PATTERN1024;
	else if (info.size() == 128)
		pt = PATTERN256;
	else if (info.size() == 4)
		pt = PATTERN8;
	else if (info.size() == 2)
		pt = PATTERN4;
	else 
	{
		std::cout << "Only pattern 4/128/512 served!\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<uint> ret(info.size()*2, 0);
	for (std::size_t i = 0, k = 0; (i<ret.size()) && (k<info.size()); i++)
		if (pt[i]) ret[i] = info[k++];

	return ret;
}

/*
uint
Polar_encoder::inv_permutation(uint pos)
{
	for (int i = 0; i < log2N/2; ++i)
	 	pos = ( pos & ~(1U<<i) & ~(1U<<(log2N-1-i)) )
			  | ( ((pos & (1U<<i))) << log2N-1-2*i )
			  | ( (pos & (1U<<(log2N-1-i))) >> log2N-1-2*i );

	return pos;
}
*/