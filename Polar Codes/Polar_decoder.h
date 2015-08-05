#ifndef POLAR_DECODER_H
#define POLAR_DECODER_H

#include"Polar.h"
#include"Channel.h"
#include<vector>
#include<string>
#include<stack>

class Polar_decoder : public Polar
{
public:
	Polar_decoder(const std::vector<double> &_rec_codeword, double _code_rate = 0.5, double _sigma2 = 1.0);
	Polar_decoder(const Channel &_chan, double _code_rate = 0.5);
	Polar_decoder();
	~Polar_decoder();
	void Decode(std::string _method, uint _search_width = 1, uint _CRC_para_r = 24); // _method = SC, SCL, CA-SCL
private:
	//----------------Interface of several decoding algorithms -------------
	void decode_CASCL(uint _search_width = 1, uint _CRC_para_r = 24);

	//---------------- functions of SCL ------------------------------------
	void decode_SCL(uint _search_width);
	void extendPath_UnfrozenBit(uint global_phase, uint _L);
	void extendPath_FrozenBit(uint global_phase, uint _L);
	uint clonePath(uint path_ind_cloned);
	void killPath(uint path_ind_killed);

	void update_BAT_List(uint global_phase, uint _L);
	void update_LLR_List(uint global_phase, uint _L);
	uint assignInitPath();

	//---------------- data struct of SCL ----------------------------------
	void init_data_struct_of_SCL(uint _L);
	void set_data_struct_of_SCL(uint _L);
	void free_data_struct_of_SCL(uint _L);
	void show_data_struct_of_SCL(uint _L, uint _cur_phase);
	void show_fixed_struct_of_SCL(uint _L);
	void show_array_struct_of_SCL(uint _L);
	void show_code_struct_of_SCL(uint _L, uint _cur_phase);

	//---------------- SC ---------------------------------------------------
//	void decode_SC();
//	void update_BAT(uint global_phase, uint path);
//	void update_LLR(uint global_phase, uint path);
	void compute_channel_llr(double *_y_in, uint _length, double *_z_out);	// layer = 0, _length = N
	void compute_inner_llr(double *_llr_in, uint _length, double *_llr_out, char _node_type, BOOL *_bat_arr);	// layer > 0
	double f_blue(double L1, double L2);
	double g_red(double L1, double L2, BOOL u);

public:
	double sigma2;
	std::vector<double> recCodeword;	//after gaussin channel, received code
	std::vector<uint> deCodeword;		//codeword decoded by algorithm using recCodeword
	std::vector<uint> deInfoBit;
private:
	//**************  fixed  ***************************************
	char **NodeType;				//N*(M+1) mat of type of node, red/blue
	uint *LayerBlockLen;			//M+1 array of each layer
	uint *LayerRenewed;				//M+1 array, when goes to i, which layer need to be renewed

	//**************  array  ***************************************
	double ***LLR;					//L*(M+1) mat of pointer of LLR
	BOOL ***BAT;					//L*(M+1) mat of pointer of u, bit array of assisting
	int **PathIndexToArrayIndex;	// L*(M+1)
	uint **ArrayReferenceCount;		// L*(M+1)
	std::vector<std::stack<uint>> ArrayInactive; //(M+1)*(0:L)

	//***************  code  ***************************************
	uint **EstimatedWord;			//L*N mat of estimated codeword
	double **PathMetricsOfForks;	//2*L array of extended PM
	double *PathMetrics;			//L array of PM
	BOOL **ForkActiveOrNot;			//2*L array of bool, keep or not
	BOOL *PathActiveOrNot;			//L array
	std::stack<uint> PathInactive;	//0:L
};

#endif		