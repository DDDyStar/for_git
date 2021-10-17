#include"read_data.h"

void read_data(string filepath, int city_num, vector<City*> *city_list) {
    ifstream infile;
    infile.open(filepath, ios::in);


    for (int i = 1; i <= city_num; i++)
    {
        int t1;
        double t2, t3;
        infile >> t1 >> t2 >> t3;
        // cout<<t2<<" "<<t3<<endl;
        city_list->push_back(new City(t1, t2, t3));
    }
    infile.close();
}

