#ifndef ICE_PROCESS_H
#define ICE_PROCESS_H

#include"ice_type.h"
#include<time.h>

//multi-thread version
void run_in_a_ebno(const std::string &_method, uint _code_len, double _ebno, uint _counts, long double &__FER, long double &__BER);
uint run_a_codeword(const std::string &_method, uint _code_len, double _sigma2);	//return val: number of error bits

#ifdef LJQ_COM
#else
void plot_curve(const std::string &_method, const double _ebno[], 
							const long double _BER[], const long double _FER[], 
							const uint _code_len, const uint _snr_cnt);
#endif

void save_simulation_result(const std::string &_method, const double _ebno[], 
							const long double _BER[], const long double _FER[], 
							const uint _Counts_Of_Each_Ebno[],
							const uint _code_len, const uint _snr_cnt, 
							const time_t &_start, const time_t &_end);

void save_run_time(std::ofstream& fname, const time_t &_start, const time_t &_end);

#endif