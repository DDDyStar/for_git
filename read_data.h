#pragma once
#include<iostream>
#include<fstream>
#include<vector>
#include"CTSP_element.h"
#include<string>
using namespace std;

void read_data(string filepath, int city_num, vector<City*> *city_list);