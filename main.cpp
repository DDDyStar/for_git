#include<iostream>
#include"read_data.h"
#include"CTSP_element.h"
#include<string>
#include"math.h"
using namespace std;

string filepath = "F:/myVNS/src/eil51.txt";
int CITY_NUM = 51;
int lambda = 5;
int** DIS_MATRIX;

// void delete_city(vector<City*> *c_list){
// 	for(int i=0; i<c_list->size(); i++){
// 		delete c_list->at(i);
// 	}
// }

double calc_distance(City *c1, City *c2){
	double dis = sqrt((c1->x-c2->x)*(c1->x-c2->x) + (c1->y-c2->y)*(c1->y-c2->y));
	return dis;
}

void create_dis_matrix(vector<City*> *c_list){
	DIS_MATRIX = new int *[CITY_NUM];
	for(int i=0; i<CITY_NUM; i++){
		DIS_MATRIX[i] = new int[CITY_NUM];
	}

	for(int i=0; i<CITY_NUM; i++){
		for(int j=0; j<CITY_NUM; j++){
			if(j<i){
				DIS_MATRIX[i][j] = DIS_MATRIX[j][i];
			}
			else if(i==j){
				DIS_MATRIX[i][j] = 0;
			}
			else{
				DIS_MATRIX[i][j] = calc_distance(c_list->at(i), c_list->at(j));
			}
		}
	}
}

int main() {
	vector<City*> city_list;
	read_data(filepath, CITY_NUM, &city_list);
	create_dis_matrix(&city_list);
	srand(2);
	TspSolution ts1(15);

	ts1.init();
	ts1.GMI();
	ts1.local_search();

	cout << "em test" << endl;
	for(int i=0; i<5; i++)
	{
		ts1.EM();
		ts1.print_seq();
		cout << "cost: " << ts1.cost << endl;
	}

	return 0;
}