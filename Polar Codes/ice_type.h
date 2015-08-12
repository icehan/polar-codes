#ifndef ICE_TYPE_H

#define ICE_TYPE_H

#include<string>

#define ADD_CRC
//#define LJQ_COM
//#define DEBUG

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned int BOOL;

typedef struct st_path_metric_fork {
	uint path;
	BOOL u_bit;
	double metric;
}St_PMF, *St_PMF_Ptr;

typedef struct st_path_metric {
	uint path;
	double metric;
}St_PM, *St_PM_Ptr;

typedef struct st_code_info {
	std::string method_of_decoding;
	uint search_width;
	uint code_len;
}St_CodeInfo, *St_CodeInfo_Ptr;

#define TRUE 1
#define FALSE 0
#define INF -999
#define SPACE_WARNING(value) \
	 {\
		 if (value == NULL){\
			 std::cout << "space not enough!\n";\
			 std::exit(EXIT_FAILURE);\
		 }\
	 }
#define FREE_SPACE(value)\
	{\
		if(NULL != value){\
			free(value);\
		}else{\
			std::cout << "can't free null space!\n";\
			std::exit(EXIT_FAILURE);\
		}\
	}

std::string int2str(const uint &double_temp);
std::string double2str(const long double &int_temp);
std::string ice_cur_time();


#endif