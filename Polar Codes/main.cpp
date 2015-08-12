#include"allheader.h"
#include<omp.h>
using namespace std;


int main()
{
	//L=1
	uint Counts_Of_Each_Ebno[] =  { 1000,  1000,  1000,  5000,  5000, 10000, 30000, 50000,  500000 };
	const double Eb_N0[] =        { 1.5,   1.75,  2,     2.25,  2.5,  2.75,  3,     3.25,   3.5};
/*
	//L=2
	uint Counts_Of_Each_Ebno[] =  { 1000,  1000,  5000,  10000, 50000, 100000, 500000};
	const double Eb_N0[] =        { 1.5,   1.75,  2,     2.25,  2.5,   2.75,   3};
	//L=4
	uint Counts_Of_Each_Ebno[] =  { 1000,  1000,  5000,  5000,  10000, 50000,  200000 };
	const double Eb_N0[] =		  { 1.0,   1.25,  1.5,   1.75,  2.0,   2.25,   2.5 };
	//L=8
	uint Counts_Of_Each_Ebno[] =  { 1000,  1000,  5000,  10000, 50000, 200000 };
	const double Eb_N0[] =		  { 1.0,   1.25,  1.5,   1.75,  2.0,   2.25};
	//L=16
	uint Counts_Of_Each_Ebno[] =  { 1000,  2000,  5000,  20000, 100000, 500000 };
	const double Eb_N0[] =		  { 1.0,   1.25,  1.5,   1.75,  2.0,    2.25 };
	//L=32
	uint Counts_Of_Each_Ebno[] =  { 1000,  1000,  2000,  10000, 50000, 300000 };
	const double Eb_N0[] =        { 0.75,  1.0,   1.25,  1.5,   1.75,  2.0};
*/
	long double FER[sizeof(Eb_N0)/sizeof(double)] = {0};

	St_CodeInfo st_code_info;
	st_code_info.method_of_decoding = "CA-SCL";
	st_code_info.search_width = 1;
	st_code_info.code_len = 1024;

	//simulation
	time_t t_start = time(NULL);
	for (int i = 0; i < sizeof(Eb_N0) / sizeof(double); i++)	// loop for every Eb/N0
		run_in_a_ebno(st_code_info, Counts_Of_Each_Ebno[i], Eb_N0[i], FER[i]);
	time_t t_end = time(NULL);

	//write results into file
	save_simulation_result(st_code_info, sizeof(Eb_N0) / sizeof(double), Counts_Of_Each_Ebno, Eb_N0, FER, t_start, t_end);

#ifndef LJQ_COM
	//plot curve using matlab
	plot_curve(st_code_info, sizeof(Eb_N0) / sizeof(double), Eb_N0, FER);
#endif

	return 0;
}

/*
int main()
{
	#pragma omp parallel for
	for (int i = 0; i < 1; i++)
	{
		cout << i << "\n";
		InfoSource src(128 - 24); // 2 4 128 512

		CRC crc(src);
		crc.Encode();

		Polar_encoder polenc(crc.codeword);
		polenc.Encode(); 

		Modulation modu(polenc.codeword);
		modu.bpsk();

//		Channel chan(modu.infoOut, 0, 0.05);
//		chan.add_gauss();

		Polar_decoder poldec(modu.infoOut);
		poldec.Decode("CA-SCL", 4);

#ifdef DEBUG
		cout << "info + crc:\n";
		for (size_t i = 0; i < crc.N; i++)
			cout << crc.codeword[i] << " ";
		cout << endl;
#endif	
	}
	return 0;
}
*/