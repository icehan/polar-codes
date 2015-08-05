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
		decode_SCL(1);
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


//------------------------- I'm working on SCL -----------------------------------
void
Polar_decoder::decode_SCL(uint _L)
{
	std::cout << "##############  INIT       ################\n";

	init_data_struct_of_SCL(_L);
	uint init_path = assignInitPath();

	show_fixed_struct_of_SCL(_L);
	std::cout << "##############  PROCESSING ################\n";

	for (int global_phase = 0; global_phase < N; global_phase++)
	{
		//compute llr of current global_phase
		update_BAT_List(global_phase, _L);
		update_LLR_List(global_phase, _L);

		// extend  paths according to llr
		if (0 == pattern[global_phase]) //frozen bit
			extendPath_FrozenBit(global_phase, _L);
		else							//info bit
			extendPath_UnfrozenBit(global_phase, _L);

		// find best codeword which passed CRC in the list


		//show_data_struct_of_SCL(_L, log2N, N, cur_phase);
	}

	free_data_struct_of_SCL(_L);
	std::cout << "##############  END        ################\n";
}
void
Polar_decoder::extendPath_UnfrozenBit(uint global_phase, uint _L)
{
	uint active_path_cnt = 0;

	//populate PathMetricsOfForks
	for (uint path = 0; path < _L; path++)
	{
		if (TRUE == PathActiveOrNot[path])
		{
			//get llr of bit in cur phase of this path
			double llr = LLR[PathIndexToArrayIndex[path][log2N]][log2N][0];
			//PathMetricsOfForks
			if (llr < 0) {
				PathMetricsOfForks[path].u_bit = 0;
				PathMetricsOfForks[path].path = path;
				PathMetricsOfForks[path].metric = PathMetrics[path] + abs(llr);
				PathMetricsOfForks[path + _L].u_bit = 1;
				PathMetricsOfForks[path + _L].path = path;
				PathMetricsOfForks[path + _L].metric = PathMetrics[path];
			}
			else {
				PathMetricsOfForks[path].u_bit = 0;
				PathMetricsOfForks[path].path = path;
				PathMetricsOfForks[path].metric = PathMetrics[path];
				PathMetricsOfForks[path + _L].u_bit = 1;
				PathMetricsOfForks[path + _L].path = path;
				PathMetricsOfForks[path + _L].metric = PathMetrics[path] + abs(llr);
			}
			active_path_cnt++;
		}
		else {
			PathMetricsOfForks[path].u_bit = 0;
			PathMetricsOfForks[path + _L].u_bit = 1;
			PathMetricsOfForks[path].path = path;
			PathMetricsOfForks[path + _L].path = path;
			PathMetricsOfForks[path].metric = INF;
			PathMetricsOfForks[path + _L].metric = INF;
		}
	}

	//copy PathMetricsOfForks
	St_PMF *copy_of_PMF = (St_PMF *)malloc(2 * _L*sizeof(St_PMF));
	for (size_t i = 0; i < 2*_L; i++){
		copy_of_PMF[i] = PathMetricsOfForks[i];
	}

	//select p largest forks, p = min(active_path_cnt*2, _L)
	qsort(PathMetricsOfForks, _L, sizeof(St_PMF), comparePM);
	uint p = std::min(_L, 2 * active_path_cnt);

	//populate ForkActiveOrNot, max p forks is active
	for (uint i = 0; i < 2*_L; i++)
	{
		uint u = PathMetricsOfForks[i].u_bit;
		uint path = PathMetricsOfForks[i].path;
		ForkActiveOrNot[u][path] = (i < p ? TRUE : FALSE);
	}

	//populate PathActiveOrNot
	//kill path whose forks both inactive
	for (uint path = 0; path < _L; path++){
		PathActiveOrNot[path] = ForkActiveOrNot[0][path] | ForkActiveOrNot[1][path];
		FALSE == PathActiveOrNot[path] ? killPath(path) : 0;
	}
	
	//continue other paths and dup if necessary
	for (uint path = 0; path < _L; path++)
	{
		if (FALSE == PathActiveOrNot[path])//killed
			continue; 	
		else if (TRUE == ForkActiveOrNot[0][path] & ForkActiveOrNot[1][path]) {	//both active
			//copy path struct
			uint path_copy = clonePath(path, global_phase);
			//copy word 
			for (uint i = 0; i < global_phase; i++){
				EstimatedWord[path_copy][i] = EstimatedWord[path][i];
			}
			//extend word & populate path metrics
			EstimatedWord[path][global_phase] = 0;
			PathMetrics[path] = copy_of_PMF[path].metric;
			EstimatedWord[path_copy][global_phase] = 1;
			PathMetrics[path_copy] = copy_of_PMF[path + _L].metric;
		}else if (TRUE == ForkActiveOrNot[0][path]) {	//fork 0 active
			EstimatedWord[path][global_phase] = 0;
			PathMetrics[path] = copy_of_PMF[path].metric;
		}else{	//fork 1 acitve
			EstimatedWord[path][global_phase] = 1;
			PathMetrics[path] = copy_of_PMF[path + _L].metric;
		}
	}
}
void
Polar_decoder::extendPath_FrozenBit(uint global_phase, uint _L)
{
	double llr = 0;
	for (uint path = 0; path < _L; path++){
		if (TRUE == PathActiveOrNot[path]) {
			//get llr of bit in cur phase of this path
			llr = LLR[PathIndexToArrayIndex[path][log2N]][log2N][0];
			//update Path Metric
			llr < 0 ? PathMetrics[path] += abs(llr) : 0;
			//extend Path
			EstimatedWord[path][global_phase] = 0;
		}
	}
}
void
Polar_decoder::update_BAT_List(uint cur_phase, uint _L)
{
	uint layer_renewed = LayerRenewed[cur_phase];
	uint block_len = LayerBlockLen[layer_renewed];
	uint pos_start = cur_phase - block_len;
	uint array_ind_of_path_layer = 0;

	for (uint path = 0; path < _L; path++)
	{
		if (TRUE == PathActiveOrNot[path])
		{
			array_ind_of_path_layer = PathIndexToArrayIndex[path][layer_renewed];
			universal_encode(
				&EstimatedWord[path][pos_start],
				block_len,
				BAT[array_ind_of_path_layer][layer_renewed]
				);
		}
	}
}
void
Polar_decoder::update_LLR_List(uint global_phase, uint _L)
{
	uint layer_start_renewed = LayerRenewed[global_phase];
	uint arr_ind_of_path_layer = 0;
	uint arr_ind_of_path_prelayer = 0;

	for (uint path = 0; path < _L; path++)
	{
		if (TRUE == PathActiveOrNot[path])		//for each active path
		{
			for (uint cur_layer = layer_start_renewed;
					    cur_layer < log2N + 1; 
						cur_layer++) {			//for each layer of the active path
				arr_ind_of_path_layer = PathIndexToArrayIndex[path][cur_layer];
				for (uint cur_phase = 0;
							cur_phase < LayerBlockLen[cur_layer]; 
							cur_phase++) {		//for each phase of the current layer of current active path
					switch (cur_layer)
					{
					case 0:						//when in layer 0
						compute_channel_llr(
							recCodeword.data(),
							N,
							LLR[arr_ind_of_path_layer][cur_layer]
							);
						break;
					default:					//when in layer>0
						arr_ind_of_path_prelayer = PathIndexToArrayIndex[path][cur_layer-1];
						compute_inner_llr(
							LLR[arr_ind_of_path_prelayer][cur_layer - 1],
							LayerBlockLen[cur_layer],
							LLR[arr_ind_of_path_layer][cur_layer],
							NodeType[global_phase][cur_layer],
							BAT[arr_ind_of_path_layer][cur_layer]
							);
						break;
					}
				}
			}
		}
	}

}
uint
Polar_decoder::clonePath(uint _path_cloned, uint global_phase)
{
	//first, get a inactive path index from PathInactive.
	//Then, set this path ind to TRUE. reps. convert to active
	uint _path_copy = PathInactive.top(); PathInactive.pop();
	PathActiveOrNot[_path_copy] = TRUE;

	//as for arr struct.
	//first, get the renewed layer p in next globla phase.
	//then, for 0£ºp - 1, two path simple share the array, 
	//we just point the same array & modify cnt arr +1.
	uint _arr_copy_from_layer = LayerRenewed[global_phase + 1];
	for (uint layer = 0; layer < _arr_copy_from_layer; layer++) {
		uint arr_ind_of_cloned = PathIndexToArrayIndex[_path_cloned][layer];
		PathIndexToArrayIndex[_path_copy][layer] = arr_ind_of_cloned;
		ArrayReferenceCount[arr_ind_of_cloned][layer]++;
	}
	//but,  for p:log2N+1, each path have their own array.
	//so, we should point diff array & set cnt to 1 
	for (uint layer = _arr_copy_from_layer; layer < log2N + 1; layer++){
		uint arr_ind = ArrayInactive[layer].top(); ArrayInactive[layer].pop();
		PathIndexToArrayIndex[_path_copy][layer] = arr_ind;
		ArrayReferenceCount[arr_ind][layer] = 1;
	}

	//in the end, return the index of the copy of that path being cloned.
	return _path_copy;
	//so, we can safely update the arrary in the next global phase, and don't
	//worry about conflict. the most important is that we avoid copy of arr of
	//diff paths which shared the same arr. this is the so called 'lazy copy'.
}
void
Polar_decoder::killPath(uint path_ind_killed)
{
	PathActiveOrNot[path_ind_killed] = FALSE;
	PathInactive.push(path_ind_killed);
	for (uint layer = 0; layer < log2N + 1; layer++) {
		uint arr_ind = PathIndexToArrayIndex[path_ind_killed][layer];
		uint cnt = --ArrayReferenceCount[path_ind_killed][layer];
		if (0 == cnt) {
			ArrayInactive[layer].push(arr_ind);
		}
	}
}
uint
Polar_decoder::assignInitPath()
{
//	if (PathInactive.empty()){
//		std::cout << "Stack of inactive path is empty!\n";
//		std::exit(EXIT_FAILURE);
//	}
	uint init_path = PathInactive.top();
					 PathInactive.pop();

	PathActiveOrNot[init_path] = TRUE;
	for (uint layer = 0; layer < log2N+1; layer++)
	{
		uint s = ArrayInactive[layer].top();
		PathIndexToArrayIndex[init_path][layer] = s;
		ArrayReferenceCount[init_path][s] = 1;
	}

	return init_path;
}


//----------------------------- SC --------------------------------------
/*void
Polar_decoder::decode_SC()
{
	uint PATH = 0;
	uint L = 1;

	init_data_struct_of_SCL(L);

	//std::cout << "##############  INIT       ################\n";
	//show_fixed_struct_of_SCL(L);
	//std::cout << "##############  PROCESSING ################\n";

	for (int cur_phase = 0; cur_phase < N; cur_phase++)
	{
		//update B LLR
		update_BAT(cur_phase, PATH);
		update_LLR(cur_phase, PATH);

		//update PM && Path extend
		double cur_llr = LLR[PATH][log2N][0];
		if (0 == pattern[cur_phase]) { // frozen bit
			EstimatedWord[PATH][cur_phase] = 0;
			if (cur_llr < 0)
				PathMetricsOfForks[0][PATH] += abs(cur_llr);
		}
		else {						 // info bit
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
	free_data_struct_of_SCL(L);
}
void
Polar_decoder::update_BAT(uint global_phase, uint path)
{
	//[path=0][LayerRenewed[cur_phase]][*]
//	if (cur_phase <= 0)
//		return;
	uint layer_renewed = LayerRenewed[global_phase];
	uint block_len = LayerBlockLen[layer_renewed];
	uint pos_start = global_phase - block_len;

	universal_encode(
		&EstimatedWord[path][pos_start], 
		block_len, 
	    BAT [PathIndexToArrayIndex[path][layer_renewed]] [layer_renewed]
	);
}
void
Polar_decoder::update_LLR(uint global_phase, uint path)
{
	//update LLR[path=0][LayerRenewed[cur_phase]:M][*]
	uint layer_start_renewed = LayerRenewed[global_phase];

	for (int cur_layer = layer_start_renewed; cur_layer < log2N+1; cur_layer++){
			switch (cur_layer)
			{
			case 0:
				compute_channel_llr(
					recCodeword.data(),
					N,
					LLR[PathIndexToArrayIndex[path][cur_layer]][cur_layer]
				);
				break;
			default:
				compute_inner_llr(
					LLR[path][cur_layer - 1],
					LayerBlockLen[cur_layer],
					LLR[path][cur_layer],
					NodeType[global_phase][cur_layer],
					BAT[path][cur_layer]
					);
				break;
			}
		
	}
}
*/
void
Polar_decoder::compute_inner_llr(double *_llr_in, uint _length, double *_llr_out, char _node_type, BOOL *_bat_arr)
{
	switch (_node_type)
	{
	case 'f':
		for (int phase = 0; phase < _length; phase++)
			_llr_out[phase] = f_blue(
				_llr_in[phase], 
				_llr_in[phase + _length]);
		break;
	case 'g':
		for (int phase = 0; phase < _length; phase++)
			_llr_out[phase] = g_red (
				_llr_in[phase], 
				_llr_in[phase + _length], 
				_bat_arr[phase]);
		break;
	default:
		std::cout << "Node Type WRONG!\n";
		std::exit(EXIT_FAILURE);
		break;
	}
}
void
Polar_decoder::compute_channel_llr(double *_y_in, uint _length, double *_z_out)
{
	for (int i = 0; i <	N; i++)
		_z_out[i] = 2*_y_in[i]/sigma2;
}
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
int 
comparePM(const void* pm1, const void* pm2)
{
	St_PMF *p1 = (St_PMF *)pm1, 
		   *p2 = (St_PMF *)pm2;
	return p1->metric - p2->metric;	//descending order
}


//--------------------------  data_struct  ----------------------------------
void
Polar_decoder::init_data_struct_of_SCL(uint L)
{
	uint M = log2N;
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
	PathMetricsOfForks = (St_PMF*) malloc ( 2*L*sizeof(St_PMF) ); SPACE_WARNING(PathMetricsOfForks);
	//KeptDescendant -- 2*L
	ForkActiveOrNot = (BOOL **) malloc ( 2 * sizeof(BOOL*) ); SPACE_WARNING(ForkActiveOrNot);
	for (int i = 0; i < 2; i++){
		ForkActiveOrNot[i] = (BOOL *) malloc ( L * sizeof(BOOL) ); SPACE_WARNING(ForkActiveOrNot[i]);
	}
	//ActivePath -- L
	PathActiveOrNot = (BOOL *) malloc ( L*sizeof(BOOL) ); SPACE_WARNING(PathActiveOrNot);

	//set initial value to struct
	set_data_struct_of_SCL(L);
}
void
Polar_decoder::free_data_struct_of_SCL(uint L)
{
	uint M = log2N;
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
	free(PathMetricsOfForks);
	for (int i = 0; i < 2; i++)	{
		free(ForkActiveOrNot[i]);
	} free(ForkActiveOrNot);

	free(PathActiveOrNot);
}
void
Polar_decoder::set_data_struct_of_SCL(uint L)
{
	uint M = log2N;
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
	for (int path = 0; path < L; path++)
		PathMetrics[path] = 0;
	//PathMetricsOfForks -- 2*L
	for (uint u = 0; u < 2; u++){
		for (uint path = 0; path < L; path++)
		{
			uint ind = u*L + path;
			PathMetricsOfForks[ind].u_bit = u;
			PathMetricsOfForks[ind].path = path;
			PathMetricsOfForks[ind].metric = 0;
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
Polar_decoder::show_fixed_struct_of_SCL(uint L)
{
	uint M = log2N;
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
Polar_decoder::show_array_struct_of_SCL(uint L)
{	
	uint M = log2N;
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
Polar_decoder::show_code_struct_of_SCL(uint L, uint _cur_phase)
{
	uint M = log2N;
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
	for (uint u = 0; u < 2; u++){
		for (uint path = 0; path < L; path++){
			uint ind = u*L + path;
			std::cout << PathMetricsOfForks[ind].metric << "\t";
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
Polar_decoder::show_data_struct_of_SCL(uint L, uint _phase)
{
	show_fixed_struct_of_SCL(L);
	show_array_struct_of_SCL(L);
	show_code_struct_of_SCL(L,_phase);
}
