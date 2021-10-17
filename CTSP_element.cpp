#include"CTSP_element.h"

extern int** DIS_MATRIX;
extern int lambda;

void cost_check(vector<int>& seq, int n) {
	float cost = 0;
	for (int i = 0; i < n; i++) {
		if (i < n - 1) {
			cost += DIS_MATRIX[seq[i]][seq[i + 1]];
		}
		else {
			cost += DIS_MATRIX[seq[i]][seq[0]];
		}
	}
	cout << cost << endl;
}

City::City(int id, double x, double y, vector<int>* colors) {
	this->city_id = id;
	this->x = x;
	this->y = y;
	this->colors = colors;
}

City::City(int id, double x, double y){
	this->city_id = id;
	this->x = x;
	this->y = y;
}

City::~City(){
	cout << "delete" << endl;
}

void City::print_color() {
	cout << "color of  city " << this->city_id << " : ";
	for (int i = 0; i < this->colors->size(); i++) {
		cout << this->colors->at(i) << " ";
	}
	cout << endl;
}

TspSolution::TspSolution(int n) {
	this->num = n;
}

void TspSolution::init(){
	for(int i=0; i<this->num; i++){
		this->seq.push_back(i);
		if(i<this->num-1){
			this->cost += DIS_MATRIX[i][i+1];
		}
		else{
			this->cost += DIS_MATRIX[i][0];
		}
	}
}

void TspSolution::print_seq() {
	for (int i = 0; i < this->num - 1; i++) {
		cout << this->seq[i] << " ";
	}
	cout << endl;
}

int get_randint(int low, int up){
	// get random int within [low, up)
	int res;
	res = (rand() % (up-low)) + low;
	return res;
}

void pick_and_insert(int uj, int uj_prim, vector<int> *seq){
	int temp_id = seq->at(uj);
	for(int i=uj; i<uj_prim; i++){
		seq->at(i) = seq->at(i+1);
	}
	seq->at(uj_prim) = temp_id;
}

void TspSolution::GMI(){
	vector<int> seq_pp = this->seq;
	double cost_pp = this->cost;
	int index_uj;
	int index_uj_prim;
	int ns = 5;

	int temp_city_id;
	
	this->omega = false;
	for(int i=0; i<lambda; i++){
		int u0 = get_randint(1, this->num-ns);
		for(int j=0; j<ns-1; j++){
			
			index_uj = u0 + j;
			index_uj_prim = get_randint(index_uj+1, this->num);

			// insert process
			// cout << "---" << index_uj << " " << index_uj_prim << "---" << endl;
			// update cost
			if(index_uj_prim<this->num-1){
				cost_pp = cost_pp - DIS_MATRIX[seq_pp.at(index_uj)][seq_pp.at(index_uj+1)] - DIS_MATRIX[seq_pp.at(index_uj-1)][seq_pp.at(index_uj)] + DIS_MATRIX[seq_pp.at(index_uj-1)][seq_pp.at(index_uj+1)] 
						- DIS_MATRIX[seq_pp.at(index_uj_prim)][seq_pp.at(index_uj_prim+1)] + DIS_MATRIX[seq_pp.at(index_uj_prim)][seq_pp.at(index_uj)] + DIS_MATRIX[seq_pp.at(index_uj)][seq_pp.at(index_uj_prim+1)];
			}
			else{
				cost_pp = cost_pp - DIS_MATRIX[seq_pp.at(index_uj)][seq_pp.at(index_uj+1)] - DIS_MATRIX[seq_pp.at(index_uj-1)][seq_pp.at(index_uj)] + DIS_MATRIX[seq_pp.at(index_uj-1)][seq_pp.at(index_uj+1)] 
						 - DIS_MATRIX[seq_pp.at(index_uj_prim)][seq_pp.at(0)] + DIS_MATRIX[seq_pp.at(index_uj_prim)][seq_pp.at(index_uj)] + DIS_MATRIX[seq_pp.at(index_uj)][seq_pp.at(0)];
			}

			pick_and_insert(index_uj, index_uj_prim, &seq_pp);

			// cout << "cost: " << this->cost << " " << cost_pp << endl;
			if(cost_pp < this->cost){
				cout<< "optimization" << endl;
				for(int k=0; k<this->num; k++){
					this->seq.at(k) = seq_pp.at(k);
				}
				this->cost = cost_pp;
				this->omega = true;
			}
			
			
			// for(int k1=0; k1<this->num; k1++){
			// 	cout<< seq_pp.at(k1) << " ";
			// }
			// cout << endl;
		}
	}
}

void flip_seq(vector<int> &seq, int i, int j){
	int len = j-i;
	vector<int> R_part;
	for(int k1=0; k1 <= len; k1++){
		R_part.push_back(seq[i + len - k1]);
	}

	for(int k2=0; k2<=len; k2++){
		seq[i + k2] = R_part[k2];
	}
}

void TspSolution::EM(){
	double u = rand()/double(RAND_MAX);
	if(!this->omega || u>0.5){
		int u1 = get_randint(1, this->num-1);
		int u2 = get_randint(u1+1, this->num);

		if(u2<this->num-1){
			this->cost = this->cost - DIS_MATRIX[this->seq.at(u1-1)][this->seq.at(u1)] - DIS_MATRIX[this->seq.at(u1)][this->seq.at(u1+1)] + DIS_MATRIX[this->seq.at(u1-1)][this->seq.at(u2)] + DIS_MATRIX[this->seq.at(u2)][this->seq.at(u1+1)]
							- DIS_MATRIX[this->seq.at(u2-1)][this->seq.at(u2)] - DIS_MATRIX[this->seq.at(u2)][this->seq.at(u2+1)] + DIS_MATRIX[this->seq.at(u2-1)][this->seq.at(u1)] + DIS_MATRIX[this->seq.at(u1)][this->seq.at(u2+1)];
		}
		else{
			this->cost = this->cost - DIS_MATRIX[this->seq.at(u1-1)][this->seq.at(u1)] - DIS_MATRIX[this->seq.at(u1)][this->seq.at(u1+1)] + DIS_MATRIX[this->seq.at(u1-1)][this->seq.at(u2)] + DIS_MATRIX[this->seq.at(u2)][this->seq.at(u1+1)]
							- DIS_MATRIX[this->seq.at(u2-1)][this->seq.at(u2)] - DIS_MATRIX[this->seq.at(u2)][this->seq.at(0)] + DIS_MATRIX[this->seq.at(u2-1)][this->seq.at(u1)] + DIS_MATRIX[this->seq.at(u1)][this->seq.at(0)];
		}
		int temp_id = this->seq.at(u1);
		this->seq.at(u1) = this->seq.at(u2);
		this->seq.at(u2) = this->seq.at(u1);
		// cout << "EM: " << u1 << " " << u2 << endl;
		// for(int k1=0; k1<this->num; k1++){
		// 	cout<< this->seq.at(k1) << " ";
		// }
		// cout << endl;
	}
}


bool TspSolution::two_opt(){
	bool flag=false;
	int t = get_randint(0, this->num);
	t = get_randint(0, this->num);
	// cout << t << endl;
	vector<int> R_part;
	double temp_cost = this->cost;
	// temp_cost = temp_cost - DIS_MATRIX[this->seq[t]][this->seq[t+1]] + DIS_MATRIX[this->seq[t+1]][this->seq[0]]
	R_part.insert(R_part.begin(), this->seq.begin(), this->seq.begin()+t);
	R_part.insert(R_part.begin(), this->seq.begin()+t, this->seq.end());
	

	for(int i=1; i<this->num-2; i++){
		for(int j=i+1; j<this->num-1; j++){

			// update cost
			temp_cost = temp_cost - DIS_MATRIX[R_part[i - 1]][R_part[i]] - DIS_MATRIX[R_part[j]][R_part[j + 1]] + DIS_MATRIX[R_part[i - 1]][R_part[j]] + DIS_MATRIX[R_part[i]][R_part[j + 1]];
			

			if (temp_cost < this->cost) {
				flip_seq(R_part, i, j);
				this->seq = R_part;
				this->cost = temp_cost;
				j = i;
				flag = true;
			}
			else {
				temp_cost = this->cost;
			}
			
		}
	}
	return flag;
}

bool TspSolution::three_opt(){
	bool flag=false;
	int t = get_randint(0, this->num);
	int jmax, hmax;
	vector<int> R_part, temp1, temp2;
	double temp_cost = this->cost;
	R_part.insert(R_part.begin(), this->seq.begin(), this->seq.begin()+t);
	R_part.insert(R_part.begin(), this->seq.begin()+t, this->seq.end());
	vector<float> Distiction(4);
	int maxPosition;
	for(int i=0; i<this->num-4; i++){
		temp1 = {this->num-3, i+5};
		jmax = temp1[min_element(temp1.begin(), temp1.end()) - temp1.begin()];
		for(int j=i+1; j<=jmax; j++){
			temp2 = {this->num-2, j+5};
			hmax = temp2[min_element(temp2.begin(), temp2.end()) - temp2.begin()];
			for(int h=j+1; h<=hmax; h++){
				Distiction[0] = DIS_MATRIX[R_part[i]][R_part[i+1]] + DIS_MATRIX[R_part[j]][R_part[j+1]] + DIS_MATRIX[R_part[h]][R_part[h+1]] - DIS_MATRIX[R_part[i]][R_part[j]] - DIS_MATRIX[R_part[j+1]][R_part[h+1]] - DIS_MATRIX[R_part[i+1]][R_part[h]];
				Distiction[1] = DIS_MATRIX[R_part[i]][R_part[i+1]] + DIS_MATRIX[R_part[j]][R_part[j+1]] + DIS_MATRIX[R_part[h]][R_part[h+1]] - DIS_MATRIX[R_part[j]][R_part[h]] - DIS_MATRIX[R_part[i+1]][R_part[h+1]] - DIS_MATRIX[R_part[j+1]][R_part[i]];
				Distiction[2] = DIS_MATRIX[R_part[i]][R_part[i+1]] + DIS_MATRIX[R_part[j]][R_part[j+1]] + DIS_MATRIX[R_part[h]][R_part[h+1]] - DIS_MATRIX[R_part[i]][R_part[h]] - DIS_MATRIX[R_part[i+1]][R_part[j+1]] - DIS_MATRIX[R_part[h+1]][R_part[j]];
				Distiction[3] = DIS_MATRIX[R_part[i]][R_part[i+1]] + DIS_MATRIX[R_part[j]][R_part[j+1]] + DIS_MATRIX[R_part[h]][R_part[h+1]] - DIS_MATRIX[R_part[i]][R_part[j+1]] - DIS_MATRIX[R_part[j]][R_part[h+1]] - DIS_MATRIX[R_part[h]][R_part[i+1]];
				maxPosition = max_element(Distiction.begin(), Distiction.end()) - Distiction.begin();
				if(Distiction[maxPosition]<=0) { continue; }
				else{
					switch (maxPosition)
					{
					case 0:
					{
						flip_seq(R_part, i+1, j);
						flip_seq(R_part, j+1, h);
						break;
					}
					case 1:{
						flip_seq(R_part, j+1, h);
						flip_seq(R_part, i+1, h);
						break;
					}
					case 2:{
						flip_seq(R_part, i+1, j);
						flip_seq(R_part, i+1, h);
						break;
					}
					case 3:{
						flip_seq(R_part, i+1, j);
						flip_seq(R_part, j+1, h);
						flip_seq(R_part, i+1, h);
						break;
					}
					default:
						break;
					}
					h = j;
					temp_cost -= Distiction[maxPosition];
					cout << "three opt" << endl;
					flag = true;
				}
			}
		}
	}
	return flag;
}

void TspSolution::local_search(){
	if(this->omega){
		bool two_improve = this->two_opt();
		if(!two_improve){
			bool three_improve = this->three_opt();
			if(three_improve){
				this->omega = true;
			}
			else{
				this->omega = false;
			}
		}
		else{
			this->omega = true;
		}
	}
}
