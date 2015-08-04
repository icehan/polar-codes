#include"ice_type.h"
#include<sstream>
#include<string>
#include<ctime>

std::string ice_cur_time()
{
	time_t t = time(NULL); 
	tm cur_t;
	localtime_s(&cur_t, &t);
	std::stringstream stream;
	stream << 1900 + cur_t.tm_year << "_" 
		   << 1 + cur_t.tm_mon << "_" 
		   << cur_t.tm_mday << "_" 
		   << cur_t.tm_hour << "_" 
		   << cur_t.tm_min;
	return stream.str();
}


std::string double2str(const long double &double_temp)
{ 
        std::stringstream stream; 
        stream << double_temp;  
        return stream.str();   //此处也可以用 stream>>string_temp  
}

std::string int2str(const uint &int_temp)
{ 
        std::stringstream stream; 
        stream << int_temp;  
        return stream.str();   //此处也可以用 stream>>string_temp  
}