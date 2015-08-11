#ifndef ICE_TYPE_H
#define ICE_TYPE_H

#include<string>

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned int BOOL;
typedef struct st_path_metric_fork {
	uint path;
	BOOL u_bit;
	double metric;
}St_PMF;
typedef struct st_path_metric {
	uint path;
	double metric;
}St_PM;

#define TRUE 1
#define FALSE 0
#define INF -999
#define SPACE_WARNING(value) \
	 {\
		 if (value == NULL){\
			 std::cout << "space not enough";\
			 std::exit(0);\
		 }\
	 }
//#define ADD_CRC
//#define LJQ_COM
#define DEBUG


std::string int2str(const uint &double_temp);
std::string double2str(const long double &int_temp);
std::string ice_cur_time();


#endif