#include"CRC.h"

void
CRC::set_gen()
{
	//default generate polynomial
	uint coef_one[] = {24, 23, 18, 17, 14, 11, 10, 7, 6, 5, 4, 3, 1, 0};
	for (int i = 0; i < sizeof(coef_one)/sizeof(uint); ++i)
		posCoefOne.push_back(coef_one[i]);

	//real coefficient of gene: 1 1 0 ...
	coefOfGen.resize(R+1, 0);	// all 0s
	for(std::size_t i = 0; i < sizeof(coef_one)/sizeof(uint); ++i)
	{
		coefOfGen[posCoefOne[i]] = 1;
	}
}

std::vector<uint>
CRC::really_encode(const std::vector<uint>& info)
{
	//R registers set to 0
	std::vector<uint> regist(R, 0);

	//m0 m1 ... mk-1 0 0 ...0
	std::vector<uint> code = info;
	for (std::size_t i = 0; i < regist.size(); ++i)
		code.push_back(0);

	//encode
	/*
	Attention: std::size_t  unsigned
	*/
	for (std::size_t i = 0; i < info.size(); ++i)
	{
		uint tmp = regist[regist.size()-1] ^ code[i];
		for (std::size_t j = regist.size()-1; j >= 1; --j)
			regist[j] = regist[j-1] ^ (tmp & coefOfGen[j]);
		regist[0] = tmp;
	}
	
	for (std::size_t i = 0; i < regist.size(); ++i)
		code[info.size() + i] = regist[regist.size()-1-i];
	
	return code;
}

void
CRC::Encode(void)
{
	codeword = really_encode(infoBits);
}

bool
CRC::IsCorrectCodeword(void)
{
	std::vector<uint> info(codeword.size()-R);

	for (std::size_t i = 0; i < info.size(); i++)
		info[i] = codeword[i];
	std::vector<uint> new_codeword = really_encode(info);

	for (std::size_t i = info.size(); i < info.size()+R; i++)
		if (new_codeword[i] != codeword[i])
			return false;
	return true;
}

CRC::CRC(const std::vector<uint> &bits, uint r, std::string gen_ploy, bool is_codeword)
{
	R = r;
	if(is_codeword == true)
	{
		N = bits.size();
		K = N - R;
		codeRate = 1.0*K/N;
		codeword = bits;
	}else
	{
		K = bits.size();
		N = K + R;
		codeRate = 1.0*K/N;
		infoBits = bits;
	}

	genPloy = gen_ploy;
	set_gen();
}

CRC::CRC(const InfoSource &src,  std::string gen_ploy, uint r)
{
	R = r;
	K = src.infoLen;
	N = K + R;
	codeRate = 1.0*K/N;
	infoBits = src.infoBits;
	genPloy = gen_ploy;
	set_gen();
}

CRC::CRC(uint k, std::string gen_ploy, uint r)
{
	R = r;
	K = k;
	N = k + r;
	codeRate = 1.0*K/N;
	genPloy = gen_ploy;
	set_gen();
}

CRC::~CRC() {}