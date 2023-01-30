#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

using namespace std;
double delta1, delta2;
//double time_start = 1577907515,	time_end = 1577920226,
//string folder = "jan20";
int UTS = 3, t_start, t_end;
double inf = numeric_limits<double>::infinity();
char BashCommandFolder[150];
vector <double> wobble_time;
vector <string> wobble_path;

vector <string> wobble(string param_wobble_path, double time_start, double time_end) {
	vector <string> vector_path;
	sprintf(BashCommandFolder, "%s%s%s", "readlink -e ", param_wobble_path.c_str(),"pointing_data_* > List_wobble");
	//cout << BashCommandFolder << endl;
	system(BashCommandFolder);

	ifstream fFileList;
	fFileList.open("List_wobble");
	if (!fFileList.is_open()) {
		cout << "wobble.cpp: List_wobble is not found" << endl;
		exit(0);
	}
	vector <string> FileListOuts;
	if (fFileList.is_open()) {
		while(!fFileList.eof()) {
			string path;
			fFileList >> path;
			if(fFileList.eof()) {
				cout << "wobble.cpp: List wobble size: " << wobble_path.size() << endl;
				if(wobble_path.size() == 0) {
					cout << "wobble.cpp: ccd data folder is empty" << endl;
					exit(0);
		 		}
				break;
			}
			wobble_path.push_back(path);
			int year_point = path.find("pointing_data_") + 14;
			//cout << unix_data(UTS, stoi(path.substr(50,4)), stoi(path.substr(55,2)), stoi(path.substr(58,2)), stoi(path.substr(61,2)), stoi(path.substr(64,2)), stoi(path.substr(68,2)), 0,0,0) << endl;
			//cout << unix_data(UTS, stoi(path.substr(48,4)), stoi(path.substr(53,2)), stoi(path.substr(56,2)), stoi(path.substr(59,2)), stoi(path.substr(64,2)), stoi(path.substr(69,2)), 0,0,0)/1e6 << endl;
			//cout << path.substr(year_point,4) << "\t" << path.substr(year_point + 5,2) << "\t" << path.substr(year_point + 8,2) << "\t" << path.substr(year_point + 11,2) << "\t" << path.substr(year_point + 14,2) << "\t" << path.substr(year_point + 17,2) << endl;
			time_cam t(UTS, stoi(path.substr(year_point,4)), stoi(path.substr(year_point + 5,2)), stoi(path.substr(year_point + 8,2)), stoi(path.substr(year_point + 11,2)), stoi(path.substr(year_point + 14,2)), stoi(path.substr(year_point + 17,2)), 0,0,0);
			wobble_time.push_back(t.get_unix_time());
			//cout << wobble_time.size() << endl;
		}
		t_start = -1;
		t_end = -1;
		delta1 = inf;
		delta2 = inf;
		//cout << wobble_time[0] << endl;
		for(int i = 0; i < wobble_time.size(); i++) {
			//cout << time_start << "\t" << wobble_time[i] << endl;
			if(delta1 > abs(time_start - wobble_time[i])) {
				delta1 = abs(time_start - wobble_time[i]);
				t_start = i;
				//cout << t_start << endl;
			}
			if(delta2 > abs(time_end - wobble_time[i]) && wobble_time[i] > time_end) {
				delta2 = abs(time_end - wobble_time[i]);
				t_end = i;
				//cout << t_end << endl;
			}
			if(abs(time_start - wobble_time[i]) < 12*60*60){
				vector_path.push_back(wobble_path[i]);
			}
		}
		cout << "wobble.cpp: delay between first time of ccd and out file name is: " << delta1/60 << " min" << endl;
		if(delta1 > 7200) {cout << "wobble.cpp: WARNING: pointing data file started with delay > 2h from pointing name file" << endl;}
	}
	fFileList.close();
	return vector_path;
}

/*vector <double> unix_t;
   vector <int> vector_h;
   vector <int> vector_m;
   vector <double> vector_s;
   vector <double> vector_err;
   vector <double> vector_sourcexcam;
   vector <double> vector_sourceycam;
   vector <double> vector_starxcam;
   vector <double> vector_starycam;
   vector <double> vector_tracking;
   vector <double> vector_zenit;
   vector <double> vector_good;*/
vector<vector<double> > read_ccd(vector <string> vector_path, double tim_start, double tim_end, bool clean_only){ //читаемые из data_ccd столбцы //unix_time, error_deg, tel_ra, tel_dec, tel_az, tel_el, source_ra, source_dec, source_az, source_el, source_x, source_y, star_x, star_y, tracking, good, weather, alpha_c
	vector<vector<double> > vector_ccd( 18, vector<double> (0));
	int column[18] = {0,4,5,6,7,8,13,14,15,16,17,18,23,24,25,26,27,28}; //читаемые из data_ccd столбцы //unix_time error_deg, tel_ra, tel_dec, tel_az, tel_el, source_ra, source_dec, source_az, 
										//source_el, source_x, source_y, star_x, star_y, tracking, good, weather, alpha_c
	double x;
	string source, line;
	for(int jj = 0; jj < vector_path.size(); jj++) {
		ifstream file_ccd(vector_path[jj]);
		//cout << "wobble.cpp: " << vector_path[jj] << endl;
		if (!file_ccd.is_open()) {
			cout << "wobble.cpp: Файл ccd не найден" << endl;
			exit(0);
		}
		else {
			getline(file_ccd, line);
			while (!file_ccd.eof()) {
				getline(file_ccd, line);
				if(file_ccd.eof()) {
					//cout << "\t\twobble.cpp: ccd_data is empty or ended" << endl;
					break;
				}
				stringstream ist(line);
				double u_time = 0;
				for(int i = 0; i < 29; i++) {
					getline(ist, source, ',');
					if (i == 0){
					u_time = atof(source.c_str());
					}
					if (u_time >= tim_start && u_time <= tim_end){
						//cout << source << endl;
						for(int j = 0; j < 18; j++) {
						//cout << i << "\t" << column[j] << endl;
							if(i == column[j]) {
								//cout << atof(source.c_str()) << endl;
								vector_ccd[j].push_back(atof(source.c_str()));
							}
						}
					}
					else{break;}
				}
			}
		}
		//cout << "size: " << vector_ccd[0].size() << endl;
		file_ccd.close();
		}	
	if (clean_only == 1) {	
		for(int j = 0; j < 18; j++) {	
			vector_ccd[j].push_back(0);	
		}
		vector_ccd[12][0] = 1000;	
		vector_ccd[13][0] = 1000;	
	}
	return vector_ccd;
}
