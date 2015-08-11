#ifndef CRC_H
#define CRC_H

#include"ice_type.h"
#include"InfoSource.h"
#include<vector>
#include<string>

class CRC
{
public:
	CRC(const std::vector<uint> &bits, uint r=24, std::string gen_ploy = "", bool is_codeword = false);
	CRC(const InfoSource &src, std::string gen_ploy = "", uint r=24);
	CRC(uint k, std::string gen_ploy = "", uint r=24);
	CRC();
	~CRC();
	void Encode();
	bool IsCorrectCodeword(const std::vector<uint> &_cw, uint _r=24);

public:
	std::size_t K;
	std::size_t N;
	std::size_t R;
	double codeRate;
	std::vector<uint> infoBits;
	std::vector<uint> codeword;

private:
	std::vector<uint> posCoefOne;	//{24, 23, 18, 17, 14, 11, 10, 7, 6, 5, 4, 3, 1, 0}µÄÄæÐò
	std::vector<uint> coefOfGen; //1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1µÄÄæÐò
	std::string genPloy; 

	void set_gen();
	std::vector<uint> really_encode(const std::vector<uint>& info);
//	vector<uint> string2num(const string strnums);
};

#endif