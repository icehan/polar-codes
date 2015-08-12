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
void plot_curve(const St_CodeInfo& _st_code_info, const uint _snr_cnt,
				const double _ebno[], const long double _FER[])
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
	legend += (_st_code_info.method_of_decoding + ", N = " + int2str(_st_code_info.code_len) + "');");

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
		fprintf(stderr,"Engine not found!\n");
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

void save_simulation_result(const St_CodeInfo& _st_code_info, const uint _snr_cnt,
							const uint _Counts_Of_Each_Ebno[],
							const double _ebno[], const long double _FER[],
							const time_t &_start, const time_t &_end)
{
	std::string cur_time = ice_cur_time();
	std::string name_of_output_file = "D:/Simulation Result/Polar_Res__" + cur_time + "__" 
										+  _st_code_info.method_of_decoding + "_" 
										+ int2str(_st_code_info.code_len) + ".txt";
	std::ofstream out_file(name_of_output_file);

	save_run_time(out_file, _start, _end);
	out_file << "Eb/N0\t\tFER\t\tCounts\n";
	for (int i = 0; i < _snr_cnt; i++)
		out_file << _ebno[i] << "\t\t" 
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

void run_in_a_ebno(const St_CodeInfo& _st_code_info, const uint _counts,
					const double _ebno, long double &__FER)
{

#ifdef ADD_CRC
	double sigma2 = 1.0 * pow(10.0, -_ebno/10) * (_st_code_info.code_len-48)/ _st_code_info.code_len;
#else
	double sigma2 = 1.0 * pow(10.0, -_ebno/10);
#endif

	long long counts_of_word_err = 0;

	//loop for specified Eb/N0
	//using openMP multi-core program
	#pragma omp parallel for
	for(int j = 0; j < _counts; ++j)
	{
		if ( FALSE == run_a_codeword(_st_code_info, sigma2) )
			++counts_of_word_err;

		if(!(j%100))
		{
			std::cout << "Eb/N0 = " << _ebno
					  << ", the " << j+1
					  << " FER = " << 1.0*counts_of_word_err/_counts
					  << ".\n";
		}
	}
	__FER = 1.0*counts_of_word_err/_counts;

	std::cout << "Eb/N0 = " << _ebno << "\tFER = "  << __FER << std::endl;
}

BOOL run_a_codeword(const St_CodeInfo& _st_code_info, const double _sigma2)
{
#ifdef ADD_CRC
	InfoSource src(_st_code_info.code_len/2 - 24); 
	
	CRC crc(src);
	crc.Encode();

	Polar_encoder polenc(crc.codeword);
#else 
	InfoSource src(_st_code_info.code_len /2);
	Polar_encoder polenc(src.infoBits);
#endif
	polenc.Encode();

	Modulation modu(polenc.codeword);
	modu.bpsk();

	Channel chan(modu.infoOut, 0, _sigma2);
	chan.add_gauss();

	Polar_decoder poldec(chan);
//	poldec.Decode(_st_code_info.method_of_decoding, _st_code_info.search_width);

	return poldec.decode_correct;
}
