#include"InfoSource.h"
#include"Polar_encoder.h"
#include"CRC.h"
#include"Modulation.h"
#include"Channel.h"
#include"Polar_decoder.h"
#include"ice_process.h"
#include"ice_type.h"
#include<fstream>
#include<string>
#include<cmath>
#include<omp.h>
#include<iomanip>

#ifdef LJQ_COM
#else

#include<engine.h>
void plot_curve(const std::string &_method, const double _ebno[], const long double _BER[], const long double _FER[], 
							const uint _code_len, const uint _snr_cnt)
{
	/*
	*	title('图形名称');
	*	xlabel('x轴说明');
	*	ylabel('y轴说明');
	*	x = [1,2,3];
	*	y = [1,2,3];
	*	semilogy（x, y1, r-);
	*   semilogy（x, y2, b-);
	*	legend('图例1', '图例2');
	*/
	//construct matlab statements
	std::string xlabel = "xlabel('Eb/N0(dB)');";
	std::string ylabel = "ylabel('FER');";
	std::string plot = "semilogy(x, y, 'r-');";
	std::string legend = "legend('Polar "; 
	legend += (_method + ", N = " + int2str(_code_len) + "');");

	std::string x_data = "x = [";
	std::string y_data = "y = [";
	for (uint i = 0; i < _snr_cnt; i++)
	{
		x_data += (double2str(_ebno[i]) + ",");
		y_data += (double2str(_FER[i]) + ",");
	}
	x_data[x_data.size()-1] = ']'; x_data += ";";
	y_data[y_data.size()-1] = ']'; y_data += ";";

	//invoke matlab engine
	Engine* EnginePt;
	EnginePt = engOpen(NULL);
	if( EnginePt==NULL ) { 
		fprintf(stderr,"engine not found!\n");
		exit(-1);
	}
	engEvalString(EnginePt, xlabel.c_str());
	engEvalString(EnginePt, ylabel.c_str());
	engEvalString(EnginePt, x_data.c_str());
	engEvalString(EnginePt, y_data.c_str());
	engEvalString(EnginePt, plot.c_str());
	engEvalString(EnginePt, legend.c_str());
	engEvalString(EnginePt, "grid on;");
	system("pause");
	engClose(EnginePt);
}
#endif

void save_simulation_result(const std::string &_method, const double _ebno[], const long double _BER[], const long double _FER[], 
							const uint _Counts_Of_Each_Ebno[], const uint _code_len, const uint _snr_cnt, const time_t &_start, const time_t &_end)
{
	std::string cur_time = ice_cur_time();
	std::string name_of_output_file = "D:/Simulation Result/Polar_Res__" + cur_time + "__" +  _method + "_" + int2str(_code_len) + ".txt";
	std::ofstream out_file(name_of_output_file);

	save_run_time(out_file, _start, _end);
	out_file << "Eb/N0\t\tBER\t\tFER\t\tCounts\n";
	for (int i = 0; i < _snr_cnt; i++)
		out_file << _ebno[i] << "\t\t" 
		<< std::setprecision(6)  << _BER[i] << "\t\t" 
		<< std::setprecision(6) << _FER[i] << "\t\t" 
		<< _Counts_Of_Each_Ebno[i] << "\n";

	out_file.close();
}

void save_run_time(std::ofstream& fname, const time_t &_start, const time_t &_end)
{
	long secondsPerDay = 24*60*60;
	long secondsPerHour = 60*60;
	long secondsPerMin = 60;

	long days, hours, mins, seconds;
	long during_seconds = _end - _start;
	days = during_seconds / secondsPerDay;
	hours = (during_seconds - days * secondsPerDay)/secondsPerHour;
	mins = (during_seconds - days * secondsPerDay - hours * secondsPerHour) / secondsPerMin;
	seconds = during_seconds - days * secondsPerDay - hours * secondsPerHour - mins * secondsPerMin;

	fname << "Running " << days << "d " << hours << "h " << mins << "min " << seconds << "sec\n";
	fname << "----------------------------------------\n";
}

void run_in_a_ebno(const std::string &_method, uint _code_len, double _ebno, uint _counts, long double &__FER, long double &__BER)
{

#ifdef ADD_CRC
	double sigma2 = 1.0 * pow(10.0, -_ebno/10) * (_code_len-48)/_code_len;
#else
	double sigma2 = 1.0 * pow(10.0, -_ebno/10);
#endif

	long long counts_of_bit_err = 0;
	long long counts_of_word_err = 0;

	//loop for specified Eb/N0
	//using openMP multi-core program
	#pragma omp parallel for
	for(int j = 0; j < _counts; ++j)
	{
		long long cur_round_bit_err = run_a_codeword(_method, _code_len, sigma2);
		if( cur_round_bit_err ){
			counts_of_bit_err  += cur_round_bit_err;
			++counts_of_word_err;
		}

		if(!(j%100))
		{
			std::cout << "Eb/N0 = " << _ebno
					  << ", the " << j+1
					  << "th.\t | BER = " << 1.0*counts_of_bit_err/(_code_len*_counts) 
					  << " FER = " << 1.0*counts_of_word_err/_counts
					  << ".\n";
		}
	}
	__BER = 1.0*counts_of_bit_err/(_code_len*_counts);
	__FER = 1.0*counts_of_word_err/_counts;

	std::cout << "Eb/N0 = " << _ebno << "\tBER = " << __BER << "\tFER = "  << __FER << std::endl;
}

uint run_a_codeword(const std::string &_method, uint _code_len, double _sigma2)
{
#ifdef ADD_CRC
	InfoSource src(_code_len/2 - 24); 
	CRC crc(src);
	crc.Encode();
	Polar_encoder polenc(crc.codeword);
	polenc.Encode();
#else 
	InfoSource src(_code_len/2); 
	Polar_encoder polenc(src);
	polenc.Encode();
#endif

	Modulation modu(polenc.codeword);
	modu.bpsk();

	Channel chan(modu.infoOut, 0, _sigma2);
	chan.add_gauss();

	Polar_decoder poldec(chan);
	poldec.Decode(_method);

	uint bit_err_num = 0;
#ifdef ADD_CRC
	uint LEN = poldec.K - 24;
#else
	uint LEN = poldec.K;
#endif
	for (int i = 0; i < LEN; i++)
	{
		if (polenc.infoBits[i] != poldec.deInfoBit[i])
			++bit_err_num;
	}

	return bit_err_num;
}
