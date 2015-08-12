#ifndef INFO_SOURCE_H
#define INFO_SOURCE_H

#include<iostream>
#include<vector>
#include"ice_type.h"

class InfoSource
{
public:
	InfoSource(int len);
	InfoSource(uint a[], int len);
	InfoSource();
	~InfoSource() {};
	void ShowInfo();

public:
	std::vector<uint> infoBits;
	int infoLen;
};

#endif