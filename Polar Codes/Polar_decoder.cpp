#include"Polar_decoder.h"
#include<iostream>
#include<cmath>
#include<algorithm> 

//------------------------ Constructor & Destructor ------------------------------------
Polar_decoder::Polar_decoder(const Channel &_chan, double _code_rate)
{
	codeRate = _code_rate;
	N = _chan.infoLen;
	K = N*_code_rate;
	R = N-K;
	log2N = log(N*1.0)/log(2.0);
	set_pattern(K);
	
	sigma2 = _chan.sigma2;
	recCodeword = _chan.infoOut;
	deCodeword.resize(N);
	deInfoBit.resize(K);
}
Polar_decoder::Polar_decoder(const std::vector<double> &_rec_codeword, double _code_rate, double _sigma2)
{
	codeRate = _code_rate;
	N = _rec_codeword.size();
	K = N*_code_rate;
	R = N-K;
	log2N = log(N*1.0)/log(2.0);
	set_pattern(K);

	sigma2 = _sigma2;
	recCodeword.resize(_rec_codeword.size());
	recCodeword = _rec_codeword;
	deCodeword.resize(_rec_codeword.size());
	deInfoBit.resize(K);
}
Polar_decoder::Polar_decoder() {}
Polar_decoder::~Polar_decoder() {}

//-------------------------Interface of several decoding algorithms ------------------------------
void
Polar_decoder::Decode(std::string _method, uint _L, uint _CRC_r)
{
	if (_method == "SC")
		decode_SC();
	else if (_method == "SCL")
		decode_SCL(_L);
	else if (_method == "CASCL")
		decode_CASCL(_L, _CRC_r);
	else
	{
		std::cout << "Polar decode method wrong!\n";
		exit(EXIT_FAILURE);
	}
}
void
Polar_decoder::decode_CASCL(uint _L, uint _CRC_r){}


//------------------------- I'm working on SC -----------------------------------
void
Polar_decoder::decode_SCL(uint _L)
{
std::cout << "##############  INIT       ################\n";
	init_data_struct_of_SCL(_L, log2N, N);
	set_data_struct_of_SCL(_L, log2N, N);
show_fixed_struct_of_SCL(_L, log2N, N);
std::cout << "##############  PROCESSING ################\n";
	for (int cur_phase = 0; cur_phase < N; cur_phase++)
	{

//show_data_struct_of_SCL(_L, log2N, N, cur_phase);
	}

	free_data_struct_of_SCL(_L, log2N, N);
std::cout << "##############  END        ################\n";
}
void
Polar_decoder::decode_SC()
{
	uint PATH = 0;
	uint L = 1;

	init_data_struct_of_SCL(L, log2N, N);

//std::cout << "##############  INIT       ################\n";
//show_fixed_struct_of_SCL(L, log2N, N);
//std::cout << "##############  PROCESSING ################\n";

	for (int cur_phase = 0; cur_phase < N; cur_phase++)
	{
		//update B LLR
		update_BAT(cur_phase, PATH);
		update_LLR(cur_phase, PATH);

		//update PM && Path extend
		double cur_llr = LLR[PATH][log2N][0];
		if (0 == pattern[cur_phase]){ // frozen bit
			EstimatedWord[PATH][cur_phase] = 0;
			if (cur_llr < 0)
				PathMetricsOfForks[0][PATH] += abs(cur_llr);
		}else {						 // info bit
			if (cur_llr < 0)
				EstimatedWord[PATH][cur_phase] = 0;
			else
				EstimatedWord[PATH][cur_phase] = 1;
		}
//std::cout << "##############  PHASE = " << cur_phase << " ################\n";
//show_varied_struct_of_SCL(L, log2N, N);
	}

//std::cout << "##############  END    ################\n";
	for (int i = 0, j = 0; i < N; i++)
	{
		deCodeword[i] = EstimatedWord[PATH][i];
		if (1 == pattern[i])
			deInfoBit[j++] = deCodeword[i];
	}
	free_data_struct_of_SCL(L, log2N, N);
}

//----------------------------- compute LLR & BAT --------------------------------------
void
Polar_decoder::update_BAT(uint cur_phase, uint path)
{
	//[path=0][LayerRenewed[cur_phase]][*]
	if (cur_phase <= 0)
		return;

	uint layer_renewed = LayerRenewed[cur_phase];
	uint block_len = LayerBlockLen[layer_renewed];
	uint pos_start = cur_phase - block_len;

	universal_encode(
		&EstimatedWord[path][pos_start], 
		block_len, 
	    BAT [PathIndexToArrayIndex[path][layer_renewed]] [layer_renewed]
	);
}
void
Polar_decoder::update_LLR(uint cur_phase, uint path)
{
	//update LLR[path=0][LayerRenewed[cur_phase]:M][*]
	uint layer_start_renewed = LayerRenewed[cur_phase];

	for (int layer = layer_start_renewed; layer < log2N+1; layer++){
		for (int phase = 0; phase < LayerBlockLen[layer]; phase++){
			switch (layer){
			case 0:
				compute_channel_llr(
					recCodeword.data(), N, 
					LLR[PathIndexToArrayIndex[path][layer]][layer]
				);
				break;
			default:
				compute_inner_llr(LLR, path, layer, cur_phase);
				break;
			}
		}
	}
}
void
Polar_decoder::compute_inner_llr(double ***LLR, uint _path, uint _layer, uint _global_phase)
{
	if (_layer == 0)
		return;

	uint block_len = LayerBlockLen[_layer];
	for (int phase = 0; phase < block_len; phase++)
	{
		double L1 = LLR[PathIndexToArrayIndex[_path][_layer-1]][_layer-1][phase];
		double L2 = LLR[PathIndexToArrayIndex[_path][_layer-1]][_layer-1][phase + block_len];
		uint u = BAT[PathIndexToArrayIndex[_path][_layer]][_layer][phase];

		switch (NodeType[_global_phase][_layer])
		{
		case 'f':
			LLR[PathIndexToArrayIndex[_path][_layer]][_layer][phase] = f_blue(L1, L2);
			break;
		case 'g':
			LLR[PathIndexToArrayIndex[_path][_layer]][_layer][phase] = g_red(L1, L2, u);
			break;
		default:
			std::cout << "Node Type WRONG! in layer:" << _layer << ", phase:" << phase;
			std::exit(0);
			break;
		}
	}
}
void
Polar_decoder::compute_channel_llr(double *_y_in, uint _length, double *_z_out)
{
	if (0 == _length)
		return;
	for (int i = 0; i < _length; i++)
		_z_out[i] = 2*_y_in[i]/sigma2;
}


//------------------------- Basic compute function  ----------------------------
double
Polar_decoder::f_blue(double L1, double L2)
{
	if (L1<0 && L2<0){ 
		if (L1 < L2) //-2, -1
			return -L2;
		else         //-1, -2  
			return -L1;
	}else if (L1<0 && L2>0){
		if (-L1 < L2) //-1 2
			return L1;
		else          //-2 1
			return -L2;
	}else if (L1>0 && L2>0){
		if (L1 < L2)
			return L1;
		else
			return L2;
	}else{ // L1>0 L2<0
		if (L1 < -L2) 
			return -L1;
		else 
			return L2;
	}
}
double
Polar_decoder::g_red(double L1, double L2, BOOL u)
{
	return  u ? L2 - L1 : L2 + L1;
}


//--------------------------  data_struct  ----------------------------------
void
Polar_decoder::init_data_struct_of_SCL(uint L, uint M, uint N)
{
	//***********************  fixed  ***************************************
	// NodeType -- N*(M+1)
	NodeType = (char **) malloc ( N * sizeof(char*) ); SPACE_WARNING(NodeType);
	for (int i = 0; i < N; i++){
		NodeType[i] = (char *) malloc ( (M+1) * sizeof(char) ); SPACE_WARNING(NodeType[i]);
	}
	for (int phase = 0; phase < N; phase++)	{
		NodeType[phase][0] = 'f'; // function f, blue
		for (int layer = M; layer > 0; layer--)	{
			if ( phase & (0x1 << (M-layer)))
				NodeType[phase][layer] = 'g';
			else
				NodeType[phase][layer] = 'f';
		}
	}
	// LayerBlockLen -- M+1
	LayerBlockLen = (uint *) malloc ( (M+1) * sizeof(uint) ); SPACE_WARNING(LayerBlockLen);
	for (int layer = 0; layer < M+1; layer++)
		LayerBlockLen[layer] = pow(2.0, M - layer*1.0);
	// LayerRenewed -- N
	LayerRenewed = (uint *) malloc ( N * sizeof(uint) ); SPACE_WARNING(LayerRenewed);
	for (int phase = 1; phase < N; phase++){
		int tmp = 0, ii = phase;
		while( ii%2 == 0 ){
			tmp += 1;
			ii /= 2;
		}
		LayerRenewed[phase] = M - tmp;
	} LayerRenewed[0] = 0;

	//**************  array  ***************************************
	//LLR -- L*(M+1)*(LayerBlockLen[])
	LLR = (double ***) malloc (L * sizeof(double**)); SPACE_WARNING(LLR);
	for (int i = 0; i < L; i++)	{
		LLR[i] = (double **) malloc ((M+1) * sizeof(double*)); SPACE_WARNING(LLR[i]);
		for (int j = 0; j < M+1; j++) {
			LLR[i][j] = (double *) malloc (LayerBlockLen[j] * sizeof(double)); SPACE_WARNING(LLR[i][j]);
		}
	}
	//BAT -- L*(M+1)*(LayerBlockLen[])
	BAT = (BOOL ***) malloc (L*sizeof(BOOL**)); SPACE_WARNING(BAT);
	for (int i = 0; i < L; i++){
		BAT[i] = (BOOL **) malloc ((M+1) * sizeof(BOOL*)); SPACE_WARNING(BAT[i]);
		for (int j = 0; j < M+1; j++){
			BAT[i][j] = (BOOL *) malloc (LayerBlockLen[j] * sizeof(BOOL)); SPACE_WARNING(BAT[i][j]);
		}
	}
	//PathIndexToArrayIndex -- L*(M+1)
	PathIndexToArrayIndex = (int **) malloc (L*sizeof(int*)); SPACE_WARNING(PathIndexToArrayIndex);
	for (int i = 0; i < L; i++){
		PathIndexToArrayIndex[i] = (int *) malloc ((M+1)*sizeof(int)); SPACE_WARNING(PathIndexToArrayIndex[i]);
	}
	//ArrayReferenceCount -- L*(M+1)
	ArrayReferenceCount = (uint **) malloc (L*sizeof(uint*)); SPACE_WARNING(ArrayReferenceCount);
	for (int i = 0; i < L; i++){
		ArrayReferenceCount[i] = (uint *) malloc ((M+1)*sizeof(uint)); SPACE_WARNING(ArrayReferenceCount[i]);
	}
	//ArrayInactive -- (M+1)*L
	ArrayInactive.resize(M+1);

	//***************  code  ***************************************
	//EstimatedWord -- L*N 
	EstimatedWord = (uint **) malloc ( L*sizeof(uint*)); SPACE_WARNING(EstimatedWord);
	for (int i = 0; i < L; i++)	{
		EstimatedWord[i] = (uint*) malloc ( N*sizeof(uint)); SPACE_WARNING(EstimatedWord[i]);
	}
	//PathMetrics -- L
	PathMetrics = (double *) malloc ( L*sizeof(double) ); SPACE_WARNING(PathMetrics);

	//PathMetricsOfForks -- 2*L
	PathMetricsOfForks = (double **) malloc ( 2 * sizeof(double*) ); SPACE_WARNING(PathMetricsOfForks);
	for (int i = 0; i < 2; i++){
		PathMetricsOfForks[i] = (double *) malloc ( L * sizeof(double) ); SPACE_WARNING(PathMetricsOfForks[i]);
	}

	//KeptDescendant -- 2*L
	ForkActiveOrNot = (BOOL **) malloc ( 2 * sizeof(BOOL*) ); SPACE_WARNING(ForkActiveOrNot);
	for (int i = 0; i < 2; i++){
		ForkActiveOrNot[i] = (BOOL *) malloc ( L * sizeof(BOOL) ); SPACE_WARNING(ForkActiveOrNot[i]);
	}

	//ActivePath -- L
	PathActiveOrNot = (BOOL *) malloc ( L*sizeof(BOOL) ); SPACE_WARNING(PathActiveOrNot);
}
void
Polar_decoder::free_data_struct_of_SCL(uint L, uint M, uint N)
{
	//**************  fixed  ***************************************
	for (int i = 0; i < N; i++)
		free(NodeType[i]);
	free(NodeType);
	free(LayerBlockLen);
	free(LayerRenewed);

	//**************  array  ***************************************
	for (int i = 0; i < L; i++)	{
		for (int j = 0; j < M+1; j++)
			free(LLR[i][j]);
		free(LLR[i]);
	}free(LLR);
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < M+1; j++)
			free(BAT[i][j]);
		free(BAT[i]);
	}free(BAT);

	for (int i = 0; i < L; i++){
		free(PathIndexToArrayIndex[i]);
	}free(PathIndexToArrayIndex);
	for (int i = 0; i < L; i++){
		free(ArrayReferenceCount[i]);
	}free(ArrayReferenceCount);

	//***************  code  ***************************************
	for (int i = 0; i < L; i++)
		free(EstimatedWord[i]);
	free(EstimatedWord);
	free(PathMetrics);
	for (int i = 0; i < 2; i++){
		free(PathMetricsOfForks[i]);		
	} free(PathMetricsOfForks);
	for (int i = 0; i < 2; i++)	{
		free(ForkActiveOrNot[i]);
	} free(ForkActiveOrNot);

	free(PathActiveOrNot);
}

void
Polar_decoder::set_data_struct_of_SCL(uint L, uint M, uint N)
{
	//**************  array  ***************************************
	//BAT & LLR -- L*(M+1)*(LayerBlockLen[])
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			for (int phase = 0; phase < LayerBlockLen[layer]; phase++){
				LLR[path][layer][phase] = INF;
				BAT[path][layer][phase] = -INF;
			}
		}
	}
	//PathIndexToArrayIndex -- L*(M+1)
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			PathIndexToArrayIndex[path][layer] = INF;
		}
	}
	//ArrayReferenceCount -- L*(M+1)
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			ArrayReferenceCount[path][layer] = 0;
		}
	}
	//ArrayInactive -- (M+1)*L
	for (int layer = 0; layer < M+1; layer++){
		for (uint path = 0; path < L; path++) {
			ArrayInactive[layer].push(path);
		}
	}


	//***************  code  ***************************************
	//PathMetrics
	for (int path = 0; path < L; path++){
		PathMetrics[path] = 0;
	}
	//PathMetricsOfForks -- 2*L
	for (int i = 0; i < 2; i++){
		for (int path = 0; path < L; path++){
			PathMetricsOfForks[i][path] = 0;
		}
	}
	//KeptDescendant -- 2*L
	for (int i = 0; i < 2; i++){
		for (int path = 0; path < L; path++){
			ForkActiveOrNot[i][path] = FALSE;
		}
	}
	//ActivePath -- L
	for (uint path = 0; path < L; path++){
		PathActiveOrNot[path] = FALSE;
		PathInactive.push(path);
	}
}

void
Polar_decoder::show_fixed_struct_of_SCL(uint L, uint M, uint N)
{
	// lambda & NodeType
	std::cout << "net struct\n";
	for (int i = 0; i < N; i++){
		std::cout << i << "\t" << LayerRenewed[i] << "\t";
		for (int j = 0; j < M+1; j++)
			std::cout << NodeType[i][j];
		std::cout << "\n";
	}

	//block len
	std::cout << "layer:\t";
	for (int layer = 0; layer < M+1; layer++)
		std::cout << layer << "\t";
	std::cout << "\nlength:\t";
	for (int layer = 0; layer < M+1; layer++)
		std::cout << LayerBlockLen[layer] << "\t";
	std::cout << std::endl;
}
void
Polar_decoder::show_array_struct_of_SCL(uint L, uint M, uint N)
{	
	//BAT
	std::cout << "\n###BAT\n";
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			for (int phase = 0; phase < LayerBlockLen[layer]; phase++){
				std::cout << BAT[path][layer][phase] << " ";
			} std::cout << "\n";
		}std::cout << "\n";
	}
	//LLR
	std::cout << "\n###LLR\n";
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			for (int phase = 0; phase < LayerBlockLen[layer]; phase++){
				std::cout << LLR[path][layer][phase] << " ";
			} std::cout << "\n";
		}std::cout << "\n";
	}	
	//PathIndexToArrayIndex
	std::cout << "\n###To\n";
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			std::cout << PathIndexToArrayIndex[path][layer] << " ";
		} std::cout << "\n";
	}
	//ArrayReferenceCount
	std::cout << "\n###Cnt\n";
	for (int path = 0; path < L; path++){
		for (int layer = 0; layer < M+1; layer++){
			std::cout << ArrayReferenceCount[path][layer] << " ";
		} std::cout << "\n";
	}
	//ArrayInactive
	std::cout << "\n###Array Inactive\n";
	std::stack<uint> s_tmp;
	for (int layer = 0; layer < M+1; layer++){
		while (!ArrayInactive[layer].empty()){
			uint top = ArrayInactive[layer].top();
			ArrayInactive[layer].pop();
			s_tmp.push(top);
			std::cout << top << " ";
		}std::cout << "\n";
		while (!s_tmp.empty()){
			uint top = s_tmp.top();
			s_tmp.pop();
			ArrayInactive[layer].push(top);
		}
	}
}
void
Polar_decoder::show_code_struct_of_SCL(uint L, uint M, uint N, uint _cur_phase)
{
	std::cout << "\n###EstWord\n";
	for (int path = 0; path < L; path++){
		std::cout << path << "\t";
	} std::cout << "\n";
	for (int phase = 0; phase <= _cur_phase; phase++){
		for (int path = 0; path < L; path++){
			std::cout << EstimatedWord[path][phase] << "\t";
		} std::cout << "\n";
	}
	std::cout << "\n###PM\n";
	for (int path = 0; path < L; path++){
		std::cout << PathMetrics[path] << "\t";
	}
	std::cout << "\n###PM_Forks\n";
	for (int i = 0; i < 2; i++){
		for (int path = 0; path < L; path++){
			std::cout << PathMetricsOfForks[i][path] << "\t";
		} std::cout << "\n";
	}
	std::cout << "\n###Fork Active OrNot\n";
	for (int i = 0; i < 2; i++){
		for (int path = 0; path < L; path++){
			std::cout << ForkActiveOrNot[i][path] << "\t";
		} std::cout << "\n";
	}
	std::cout << "\n###Path Active OrNot\n";
	for (int path = 0; path < L; path++){
		std::cout << PathActiveOrNot[path] << "\t";
	}
	std::cout << "\n###Path Inactive\n";
	std::stack<uint> s_tmp;
	while (!PathInactive.empty()){
		uint top = PathInactive.top(); PathInactive.pop();
		s_tmp.push(top);
		std::cout << top << " ";
	} std::cout << "\n";
	while (!s_tmp.empty()){
		uint top = s_tmp.top(); s_tmp.pop();
		PathInactive.push(top);
	}
}
void
Polar_decoder::show_data_struct_of_SCL(uint L, uint M, uint N, uint _phase)
{
	show_fixed_struct_of_SCL(L,M,N);
	show_array_struct_of_SCL(L,M,N);
	show_code_struct_of_SCL(L, M,N,_phase);
}

