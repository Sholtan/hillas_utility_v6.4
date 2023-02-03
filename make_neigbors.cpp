#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
using namespace std;
string line;
char neigh_filename[90], chmod_str[150];
int number_of_pixels;
double dist1;
vector <int> trig37;

// make_neighbours function makes a file "neighbours.IACT%02d" that has structure:
//cluster_number	anode_extended_channel_number	x_coordinate	y_coordinate	number_of_neighbors(from 3 to 6)	first_neighbor's_cluster_number    first_neighbor's_extended_channel_number ...
// returns file's name
                                      //coord_file = 'xy_turn_2019s_EP.txt'   
string make_neighbours(int iact_numb, string coord_file, int *exclud_clust, int *exclud_numb)
{
	ifstream file(coord_file);
	if (!file.is_open()) {
		cout << "Файл не найден" << endl;
		exit(0);
	}
	number_of_pixels = 0;
	while (!file.eof())
	{
		getline(file, line);
		if (file.eof())
		{
			cout << "number of pixels: " << number_of_pixels << endl;
			break;
		}
		else
		{
			number_of_pixels++; // number_of_pixels == number_of_lines
		}
	}
	file.close(); // opened just to count lines
	file.open(coord_file);
	double xy_turn[number_of_pixels][8]; // ! fixed wrong size of xy_turn (was xy_turn[number_of_pixels][7])
	for (int j = 0; j < number_of_pixels; j++)
	{
		getline(file, line);
		stringstream ist(line);
		for (int i = 0; i < 8; i++)
		{
			ist >> xy_turn[j][i];
		}
		//xy_turn[j][0] is cluster_number, xy_turn[j][3] is x_coordinate, xy_turn[j][4] is y_coordinate, xy_turn[j][5] is anode_extended_channel_number,
		// xy_turn[j][6] is dynode_extended_channel_number, xy_turn[j][7] is anode_channel_number,    xy_turn[j][1] and xy_turn[j][2] are not used
		for (int i = 0; i < 28; i++)
		{
			// cout << int(xy_turn[j][0]) << "\t" << exclud_clust[i] << "\t" << int(xy_turn[j][7]) << "\t" << exclud_numb[i] << endl;
			if ((int(xy_turn[j][0]) == exclud_clust[i] && int(xy_turn[j][7]) == exclud_numb[i]))
			{
				xy_turn[j][3] = 1000; // if pixel is excluded by param_file, make their coordinates 1000
				xy_turn[j][4] = 1000;
			}
		}
	}
	file.close();

	sprintf(neigh_filename, "neighbours.IACT%02d", iact_numb + 1);
	ofstream fout(neigh_filename);
	for (int i = 0; i < number_of_pixels; i++)
	{
		if (abs(xy_turn[i][3]) < 100 && abs(xy_turn[i][4]) < 100) // fixed rhs condition (was abs(xy_turn[i][4]) < 100)
		{ 
			vector<int> neigh;
			dist1 = 0;
			for (int j = 0; j < number_of_pixels; j++)
			{
				if (i != j)
				{
					dist1 = pow((pow((xy_turn[i][3] - xy_turn[j][3]), 2) + pow((xy_turn[i][4] - xy_turn[j][4]), 2)), 0.5); // distance between i and j pixels
					if (dist1 < 3.1)
					{ // needsclarification coordinates in xy_turn file are in centimeters???
						neigh.push_back(j);
					}
				}
			}

			fout << xy_turn[i][0] << "\t" << xy_turn[i][5] << "\t" << xy_turn[i][3] << "\t" << xy_turn[i][4] << "\t" << neigh.size() << "\t";
			//cout << "\n NEIGH.size(): " << neigh.size() << "\t";
			for (int j = 0; j < neigh.size(); j++)
			{
				// cout << neigh[j] << "\t";
				fout << xy_turn[neigh[j]][0] << "\t" << xy_turn[neigh[j]][5] << "\t";
			}
			if (neigh.size() < 6)
			{
				for (int j = 6; j > neigh.size(); j--)
				{
					fout << -1 << "\t" << -1 << "\t";
				}
			}
			// cout << endl;
			fout << endl;
			neigh.clear();
		}
	}
	fout.close();
	sprintf(chmod_str, "chmod 777 %s", neigh_filename);
	system(chmod_str);
	return neigh_filename;
}
