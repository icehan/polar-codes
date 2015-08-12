#ifndef ICE_PROCESS_H
#define ICE_PROCESS_H

#include"ice_type.h"
#include<time.h>

//multi-thread version
void run_in_a_ebno(const St_CodeInfo& _st_code_info, const uint _counts,
				   const double _ebno, long double &__FER);
BOOL run_a_codeword(const St_CodeInfo& _st_code_info, const double _sigma2);


#ifndef LJQ_COM
void plot_curve(const St_CodeInfo& _st_code_info, const uint _snr_cnt,
				const double _ebno[], const long double _FER[]);
#endif

void save_simulation_result(const St_CodeInfo& _st_code_info,
							const uint _snr_cnt,  const uint _Counts_Of_Each_Ebno[],
							const double _ebno[], const long double _FER[], 
							const time_t &_start, const time_t &_end);

void save_run_time(std::ofstream& fname, const time_t &_start, const time_t &_end);

#endif
