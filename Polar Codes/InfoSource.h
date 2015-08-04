#ifndef INFO_SOURCE_H
#define INFO_SOURCE_H

#include<iostream>
#include<vector>
#include"ice_type.h"

class InfoSource
{
public:
	void ShowInfo();
	InfoSource(uint infoLength);
	InfoSource(int a[], int len);
	~InfoSource() {};

public:
	std::vector<uint> infoBits;
	std::size_t infoLen;
};

#endif