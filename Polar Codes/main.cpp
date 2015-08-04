#include"allheader.h"
#include<omp.h>
using namespace std;

/*
int main()
{
	//256 SC
	uint Counts_Of_Each_Ebno[] =  {50,   50,   50, 50,  50, 100,  500, 500, 5000, 20000, 50000, 200000, 1000000};
	const double Eb_N0[] =        {-1.0, -0.5, 0,  0.5, 1,  1.5,  2,   2.5, 3,    3.5,   4,     4.5,    5};

	//1024 CA-SCL L = 1
//	uint Counts_Of_Each_Ebno[] =  {500, 500,  500, 1000, 1000, 10000, 10000, 50000, 50000};
//	const double Eb_N0[] =        {1.5, 1.75, 2.0, 2.25, 2.5,  2.75,  3.0,   3.25,  3.5};

	long double BER[sizeof(Eb_N0)/sizeof(double)] = {0};
	long double FER[sizeof(Eb_N0)/sizeof(double)] = {0};

	string Method_Of_Decoding = "SC";
	uint Code_Len = 256;

	//simulation
	time_t t_start = time(NULL);
	for (int i = 0; i < sizeof(Eb_N0)/sizeof(double); i++)	// loop for every Eb/N0
		run_in_a_ebno(Method_Of_Decoding, Code_Len, Eb_N0[i], Counts_Of_Each_Ebno[i], FER[i], BER[i]);
	time_t t_end = time(NULL);

	//write results into file
	save_simulation_result(Method_Of_Decoding, Eb_N0, BER, FER, Counts_Of_Each_Ebno, Code_Len, sizeof(Eb_N0)/sizeof(double), t_start, t_end);

#ifdef LJQ_COM

#else
	//plot curve using matlab
	plot_curve(Method_Of_Decoding, Eb_N0, BER, FER, Code_Len, sizeof(Eb_N0)/sizeof(double));
#endif

	return 0;
}
*/

int main()
{
//	#pragma omp parallel for
	for (int i = 0; i < 1; i++)
	{
		InfoSource src(4); // 2 4 128 512

		Polar_encoder polenc(src);
		polenc.Encode();

		Modulation modu(polenc.codeword);
		modu.bpsk();

		Channel chan(modu.infoOut, 0, 0.5);
		chan.add_gauss();

		Polar_decoder poldec(chan);//modu.infoOut);		//	//cout << i+1 << "th decoding...\n";
		poldec.Decode("SCL", 3);
	}
	return 0;
}

/* test data
		double y[] = {-1.1, -0.5, -1.6, 1.2};
		for (int i = 0; i < 4; i++)
			chan.infoOut[i] = y[i];

		int j = 0;
		for (int i = 0; i < poldec.N; i++)
		{
			if (i < src.infoLen)
			{
				cout << src.infoBits[i] << "\t" << polenc.infoBitsAddPattern[i] << "\t" << polenc.codeword[i] << "\t" 
					<< poldec.recCodeword[i] << "\t" << poldec.deCodeword[i] << "\n";
			}else
			{
				cout << "*" << "\t" << polenc.infoBitsAddPattern[i] << "\t"  << polenc.codeword[i] << "\t" 
					<< poldec.recCodeword[i] << "\t" << poldec.deCodeword[i] << "\n";
			}
			if(polenc.infoBitsAddPattern[i] != poldec.deCodeword[i])
				j++;
		}
		cout << j << endl;
*/
