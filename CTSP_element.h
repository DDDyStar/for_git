#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include"stdlib.h"
using namespace std;


struct City
{
public:
	City(int id, double x, double y, vector<int>* colors);
	City(int id, double x, double y);
	~City();
	int city_id;
	double x;
	double y;
	vector<int>* colors;
	void print_color();
};

struct Salesman
{
public:
	int salesman_id;
	int* cities_id;
	vector<int> seq;
};


class TspSolution {
public:
	TspSolution(int n);
	int num;
	bool omega = false;
	double cost = 0;
	vector<int> seq;

	void init();
	void GMI();
	void EM();
	bool two_opt();
	bool three_opt();
	void print_seq();
	void local_search();

	void VNS();
};
