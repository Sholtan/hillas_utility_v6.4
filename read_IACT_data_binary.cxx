#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <vector>
#include <sys/stat.h>
#include <iterator>
#include <map>
#include "time.cpp"
#include "wobble.cpp"
#include "event.cpp"
#include "CR.cpp"
#include "make_neigbors.cpp"
#include <sys/stat.h>

using namespace std;
char fou[190],save_background_str[190], month_year[70], fou_hillas[190],folder_outs[150], background_folder[190], clean_folder[190], param_path_str[150];
string folder;
	int bsm, cch, ch,ff,f,num,por,trig,n,kkk[6][64][25],pos[6][64][25], Nsos_array[64][25], cluster[25], nn, anti_f[25], jj, jjj,j, x, gg;
int co, number_of_pixels_cam, pix_number[64][25];
double pedp, sigp, b[64][25],bmp[64][25],ped[64][25],sig[64][25], e[25][64], sens[25][64], gain, k_adc, ecode, rel_sens, time_0 = 1, event_unix_time = 0;
string str,srr[25], timet, ped_folder[4] = {"peds", "peds.m3s", "peds.mediana","peds_median_my"}, data_path, out_data_path, hillas_table_name, clean_out_name, param_wobble_path, cleaning_type;
double hour, minute, sec, mksec, mlsec, nsec, time0, x_pos[64][25], y_pos[64][25], tim_start = 0, tim_end = 0, event_delay;
map <int, string> calendar = {{1, "jan"}, {2, "feb"},{3, "mar"},{4, "apr"},{5, "may"},{6, "jun"},{7, "jul"},{8, "aug"},{9, "sep"},{10, "oct"},{11, "nov"},{12, "dec"}};
int exclud_clust[28] = {0}, exclud_numb[28] = {0}, cleaning = -1, ped_param = -1;
char press;
string nsec_time;
bool clean_only, background_marker[64][25], save_background;

double time_start_end_binary(string run_date, string path){
	time_cam t;
	FILE * ptrFile = fopen(path.c_str(), "r" );
	if (ptrFile == NULL)
	{
		cout << "ERROR out file " << run_date << " is not found" << endl << "check out files availability" << endl;
		exit(0);
	}
	else{
		size_t bufferSize = 17;
		unsigned char * buf = (unsigned char*) malloc(sizeof(char) * bufferSize);
		fread(buf, 1, bufferSize, ptrFile);
		int idf = (int)(buf[0] + buf[1]*256);
		if(idf == 3033){
			//cout << (int)buf[8] << "\t" << (int)buf[9] << "\t" << (int)buf[10] << "\t" << (int)(buf[11] + 256*buf[12]) << "\t" << (int)(buf[13] + 256*buf[14]) << "\t" << (int)(buf[15] + 256*buf[16]) << endl;
			t.set_time_binary(run_date, (int)buf[8], (int)buf[9], (int)buf[10],
										(int)(buf[11] + 256*buf[12]), (int)(buf[13] + 256*buf[14]),
										(int)(buf[15] + 256*buf[16]));
			t.get_human_string_data_time();
			fclose(ptrFile);
			return t.get_unix_time();
		}
		else{
			cout << "out file " << run_date << " is broken" << endl << "check input files" << endl;
			fclose(ptrFile);
			exit(0);
		}
	}
	return 0;
}

int fileList(string data_path, string data_folder,string run_numb)
{
	sprintf(BashCommandFolder, "%s%s%s%s%s%s", "readlink -e ", data_path.c_str(), data_folder.c_str(), ".", run_numb.c_str(),"/outsb/*out_* > List_outs");
	system(BashCommandFolder);//taiga_2020/data2020/DATA_IACT02/2020-21.txt.v3/  /k2/DATA_IACT02/2019-20.txt.v3/
	sprintf(BashCommandFolder, "%s%s%s%s%s%s%s%s", "readlink -e ", data_path.c_str(), data_folder.c_str(), ".", run_numb.c_str(),"/", ped_folder[ped_param].c_str(), "/*ped_* > List_peds");
	system(BashCommandFolder);
	return 0;
}

int main(int argc, char **argv)
{
	////////////////////////////////// read parameters
	ifstream pParam;
	pParam.open(argv[1]);
	cout << argv[1] << endl;
	if (!pParam.is_open()) {
		cout << "param file is not found" << endl;
		return 0;
	}
	if(pParam.eof()) {
		cout << "file param is empty" << endl;
		return 0;
	}
	int IACT_numb;
	string datt, factor_file, coord_file;
	double edge1, edge2;
	getline(pParam, line);
	istringstream ist4(line);
	ist4 >> datt >> datt >> IACT_numb;
	getline(pParam, line);
	istringstream ist6(line);
	ist6 >> factor_file;
	getline(pParam, line);
	istringstream ist5(line);
	ist5 >> coord_file;
	getline(pParam, line);
	istringstream ist2(line);
	ist2 >> datt >> ped_param;
	getline(pParam, line);
	istringstream ist1(line);
	ist1 >> datt >> cleaning;
	if(cleaning == 1) {
		cleaning_type = "sig";
	}
	else if(cleaning == 0) {
		cleaning_type = "fix";
	}
	getline(pParam, line);
	istringstream ist(line);
	ist >> datt >> datt >> edge1 >> edge2;
	getline(pParam, line);
	istringstream ist3(line);
	ist3 >> datt >> datt;
	for(int i = 0; i < 28; i++) {
		ist3 >> exclud_clust[i] >> exclud_numb[i];
	}
	getline(pParam, line);
	istringstream ist70(line);
	ist70 >> datt >> clean_only;
	getline(pParam, line);
	istringstream ist700(line);
	ist700 >> datt >> datt >> datt >> datt >> datt >> save_background;
	getline(pParam, line);
	istringstream ist7(line);
	ist7 >> data_path;
	getline(pParam, line);
	istringstream ist8(line);
	ist8 >> out_data_path;
	getline(pParam, line);
	istringstream ist9(line);
	ist9 >> param_wobble_path;
	getline(pParam, line);
	vector <string> FolderList;
	vector <string> RunNumbList;

	string run_numb;
	while(!pParam.eof()) {
		folder = "";
		run_numb = "";
		getline(pParam, line);
		istringstream ist(line);
		ist >> folder >> run_numb;
		if(folder.length() == 0) break;
		FolderList.push_back(folder);
		RunNumbList.push_back(run_numb);
	}
	pParam.close();
	
///////////////////////////////////////////////////// read factors
	ifstream file0(factor_file);
	if (!file0.is_open()) {
		cout << "calibration file is not found" << endl;
		return 0;
	}
	int q = 0;
	if(IACT_numb == 0) {
		for(int i = 0; i < 10; i++) {
			getline(file0, line);
			if(file0.eof()) {
				cout << "calibration file is empty" << endl;
				return 0;
			}
		}
		while(!file0.eof()) {
			getline(file0, line);
			if(!file0.eof()) {
				istringstream ist(line);
				ist >> bsm >> cch >> gain >> gain >> ecode >> rel_sens;
				//cout << "\t" << bsm << "\t" << cch << "\t" << ecode << "\t" << rel_sens << endl;
				if(ecode > 0) {
					e[bsm][cch] = ecode;
					sens[bsm][cch] = rel_sens;
					q++;
				}
				else{
					e[bsm][cch] = 1e9;
					sens[bsm][cch] = -1e9;
				}
			}
		}
	}
	else if(IACT_numb == 1) {
		getline(file0, line);
		while(!file0.eof()) {
			getline(file0, line);
			if(!file0.eof()) {
				istringstream ist(line);
				ist >> bsm >> cch >> ecode >> rel_sens;
				//cout << "\t" << bsm << "\t" << cch << "\t" << ecode << "\t" << rel_sens << endl;
				e[bsm][cch] = ecode;
				sens[bsm][cch] = rel_sens;
				if(ecode > 0) {
					e[bsm][cch] = ecode;
					sens[bsm][cch] = rel_sens;
					q++;
				}
				else{
					e[bsm][cch] = 1e9;
					sens[bsm][cch] = -1e9;
				}
			}
		}
	}
	cout << "file calibration size (28*2*22): " << q << endl;
	file0.close();

	////////////////////////////////// read neighbours table
	ifstream file1(make_neighbours(IACT_numb, coord_file, exclud_clust, exclud_numb));
	if (!file1.is_open()) {
		cout << "file neighbours is not found" << endl;
		return 0;
	}
	q = 0;
	for(int coun = 0; coun < 25; coun++) {
					for(int count = 0; count < 64; count++) {
						pix_number[count][coun]=-1;
						for(int j = 0; j < 6; j++) {
						    kkk[j][count][coun] = -1;
						    pos[j][count][coun] = -1;
						    }
					}
				}
	while(!file1.eof()) {
		int k[6]= {0},sos[6]= {0},kk=0,Nsos=0, ii=0, kx = 0,ix = 0;
		double x=0,y=0;
		getline(file1, line);
		if(file1.eof())
		{
			break;
		}
		istringstream ist(line);
		ist >> kk >> ii;
		pix_number[ii][kk] = q;
		q++;
		ist >> x_pos[ii][kk] >> y_pos[ii][kk] >> Nsos;
		Nsos_array[ii][kk] = Nsos;
		for(int j = 0; j < 6; j++) {
			kx = 0;
			ix = 0;
			ist >> kx >> ix;
			kkk[j][ii][kk] = kx;
			pos[j][ii][kk] = ix;
		}
	}
	number_of_pixels_cam = q;
	cout << "number of using camera pixels: "<< number_of_pixels_cam << endl;
	file1.close();
	////////////////////////////////// read file list
	for ( int jl=0; jl < FolderList.size(); jl++) {
		cout << "open run: " << FolderList[jl] << "." << RunNumbList[jl] << endl;
		ifstream fFileList;
		fileList(data_path, FolderList[jl], RunNumbList[jl]); //create List_outs and List_peds
		////////////////////////////////// create Outs vector
		fFileList.open("List_outs");
		if (!fFileList.is_open()) {
			cout << "List_outs is not found" << endl;
			return 0;
		}
		vector <string> FileListOuts;
		if (fFileList.is_open()) {
			while(!fFileList.eof()) {
				string path;
				fFileList >> path;
				if(fFileList.eof() || path.length() == 0) {
					cout << "List outs size: " << FileListOuts.size() << endl;
					break;
				}
				FileListOuts.push_back(path);
			}
		}
		fFileList.close();

		////////////////////////////////// create Peds vector
		fFileList.open("List_peds");
		if (!fFileList.is_open()) {
			cout << "List_peds is not found" << endl;
			return 0;
		}
		vector <string> FileListPeds;
		if (fFileList.is_open()) {
			while(!fFileList.eof())
			{
				string path;
				fFileList >> path;
				if(fFileList.eof() || path.length() == 0) {
					cout << "List peds size: " << FileListPeds.size() << endl;
					break;
				}
				FileListPeds.push_back(path);
			}
		}
		fFileList.close();

		int dir_err = -1, k = 0;
		while(dir_err == -1){
			if(k == 0){
				sprintf(folder_outs, "%s%s.%s", out_data_path.c_str(), FolderList[jl].c_str(),RunNumbList[jl].c_str());
			}
			else{
				sprintf(folder_outs, "%s%s.%s(%d)", out_data_path.c_str(), FolderList[jl].c_str(),RunNumbList[jl].c_str(),k);
			}
			dir_err = mkdir(folder_outs, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			k++;
		}
		cout << "out directory " << folder_outs << " is created" << endl;
		////////////////////////copy file param to directiry:
		ifstream  src(argv[1], ios::binary);
		sprintf(param_path_str, "%s/%s",  folder_outs, argv[1]);
		ofstream  dst(param_path_str,   ios::binary);
		dst << src.rdbuf();

    dst << src.rdbuf();
		//find start and stop times of outs//////////////////////////////////////
		int List_size = FileListOuts.size();
		cout << "time of run start and end of portions:" << endl;
		tim_start = time_start_end_binary(FolderList[jl], FileListOuts[0]);
		tim_end = 120 + time_start_end_binary(FolderList[jl], FileListOuts[List_size - 1]);
		//sprintf(month_year, "%s%d", calendar[stoi(FolderList[jl].substr(2,2))].c_str(), stoi(FolderList[jl].substr(4,2)));
		//cout << month_year << endl;
		vector<string> vector_wobble = wobble(param_wobble_path.c_str(), tim_start, tim_end);
		cout << "input:" << endl;
		for ( int i=0; i < vector_wobble.size(); i++) {
			cout << i << "\t" << vector_wobble[i] << endl;
		}
		////////////////////////////////////////////////////////////////////////
		vector<vector<double> > vector_ccd = read_ccd(vector_wobble, tim_start, tim_end, clean_only);
		int ccd_id = 0;
		cout << "\t\tnumber of written ccd rows: " << vector_ccd[0].size() << endl;
		////////////////////////////////////////////////////////////////////////
		if(save_background == 1){
			sprintf(background_folder, "%s%s", folder_outs, "/background");
			sprintf(clean_folder, "%s%s", folder_outs, "/clean");
			mkdir(background_folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			mkdir(clean_folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
		sprintf(fou_hillas, "%s/%s.%s_out_hillas_%02.0f_%02.1f%s.csv", folder_outs, FolderList[jl].c_str(),RunNumbList[jl].c_str(), edge1, edge2, cleaning_type.c_str());
		cout << folder_outs << endl;
		ofstream fout_hillas(fou_hillas);
		fout_hillas << "por,event_numb,unix_time,unix_time_long_ns,delta_time,error_deg,tel_az,tel_el,source_az,source_el,CR100phe,CR_portion,numb_pix,size,Xc[0],Yc[0],con2,length[0],width[0],dist[0],dist[1],dist[2],skewness[0],skewness[1],skewness[2],kurtosis,alpha[0],alpha[1],alpha[2],a_axis,b_axis,a_dist[1],b_dist[1],a_dist[2],b_dist[2],tel_ra,tel_dec,source_ra,source_dec,source_x,source_y,tracking,good,star,edge,weather_mark,alpha_c" << endl;
		for ( int i=0; i < List_size; i++) {
			cout << i << "\t" << FileListOuts[i] << endl;
			FILE * ptrFile = fopen(FileListOuts[i].c_str(), "r" );
			if (ptrFile == NULL)
			{
				fputs("Ошибка файла\n", stderr);
				exit(1);
			}
			else
			{
				for(int coun = 0; coun < 25; coun++) {
					for(int count = 0; count < 64; count++) {
						ped[count][coun]=0;
						sig[count][coun]=0;
					}
				}
				int qqq = 0;
				vector <Events> vector_events;
				vector<vector<double> > vector_background( number_of_pixels_cam, vector<double> (0));
				cout << i << "\t" << FileListPeds[i] << endl;
				sprintf(fou, "%s%s%s%s%02.0f%s%02.1f%s%s%03d%s", folder_outs, "/clean/", FolderList[jl].c_str(), ".cleanout_", edge1, "_", edge2, cleaning_type.c_str(), "_", i+1, ".txt");
				ofstream fout(fou);
				sprintf(save_background_str, "%s%s%s%s%02.0f%s%02.1f%s%s%03d%s", folder_outs, "/background/", FolderList[jl].c_str(), ".background_", edge1, "_", edge2, cleaning_type.c_str(), "_", i+1, ".txt");
				ofstream fout_background(save_background_str);
				cout << "output: created files:" << endl << "\t" << fou_hillas << endl;
				if(save_background == 1){
					cout << i << "\t" << fou << endl;
					cout << i << "\t" << save_background_str << endl;
				}
				if(FileListOuts[i].compare(FileListOuts[i].size()-3, 3, FileListPeds[i], FileListPeds[i].size()-3, 3)==0) {
					cout << "por=ped" << endl;
					ifstream DataFilePeds;
					DataFilePeds.open(FileListPeds[i].c_str());
					if (DataFilePeds.is_open()) {
						//cout << "peds read" << endl;
						while (!DataFilePeds.eof())
						{
							getline(DataFilePeds, str);
							//cout << str << endl;
							if(DataFilePeds.eof()) {
								cout << "PED " << i + 1 << " is ended" << endl;
								break;
							}
							istringstream iss(str);
							sigp=-40;
							pedp=-40;
							iss >> ff >> ch >> pedp >> sigp;
							//cout << pedp << "\t" << sigp << endl;
							ped[ch][ff]=pedp;
							if(cleaning == 1) {
								sig[ch][ff]=sigp/(e[ff][ch]*sens[ff][ch]);
							}
							else{
								sig[ch][ff] = 1;
							}
							//cout << ff << "\t" << ch << "\t" << ped[ch][ff] << "\t" << sig[ch][ff] << endl;
						}
						DataFilePeds.close();
						/*if(i == 0) {
							cout << "Check all parameters and press any key to continue" << endl;
							cin >> press;
						}*/
						while(!feof(ptrFile))
    					{
							int *date;
							size_t bufferSize = 3;
							unsigned char * buf = (unsigned char*) malloc(sizeof(char) * bufferSize);
							fread(buf, 1, bufferSize, ptrFile);
							int idf = (int)(buf[0] + buf[1]*256);
							if(feof(ptrFile)) {
								cout << "outs " << i + 1 << " is ended" << endl;
								break;
							}
							else if(idf==3033){ //проверка начала события
								n = (int)buf[2];
								//cout << n << endl;
								size_t bufferSizeCluster = 142*n+3;
								unsigned char * buf_clust = (unsigned char*) malloc(sizeof(char) * bufferSizeCluster);
								fread(buf_clust, 1, bufferSizeCluster, ptrFile);
								if((int)buf_clust[142*n] == (int)buf_clust[142*n + 1] && (int)buf_clust[142*n + 1] == (int)buf_clust[142*n + 2] &&  (int)buf_clust[142*n + 2] == 255){ //проверка конца события 255 255 255
								for (int count = 0; count < 25; count++)
								{
									cluster[count]=0;
									for (int coun = 0; coun < 64; coun++)
									{
										b[coun][count]=0;
										bmp[coun][count]=0;
										background_marker[coun][count]=0;
									}
								}
								time_cam t;
								/*cout << (int)buf_clust[(int)((q-1)*142 + 5)] << "\t" << (int)buf_clust[(int)((q-1)*142 + 6)] << "\t" << (int)buf_clust[(int)((q-1)*142 + 7)] << "\t" << 
								(int)(buf_clust[(int)((q-1)*142 + 8)] + 256*buf_clust[(int)((q-1)*142 + 9)]) << "\t" << (int)(buf_clust[(int)((q-1)*142 + 10)] + 256*buf_clust[(int)((q-1)*142 + 11)]) << "\t" << 
								(int)(buf_clust[(int)((q-1)*142 + 12)] + 256*buf_clust[(int)((q-1)*142 + 13)]) << endl;*/
								t.set_time_binary(FolderList[jl], (int)buf_clust[5], (int)buf_clust[6], (int)buf_clust[7], (int)(buf_clust[8] + 256*buf_clust[9]), 
																  (int)(buf_clust[10] + 256*buf_clust[11]), (int)(buf_clust[12] + 256*buf_clust[13]));
								event_unix_time = t.get_unix_time();
								nsec_time = t.full_unix_nsec();
								//cout << nsec_time << endl;
								timet = t.char_human_string_time();
								if(qqq == 0) {
									time_0 = event_unix_time;
								}
								for( int q=1; q <= n; q++)
								{
									f= (int)buf_clust[(int)((q-1)*142)];
									for(int ib=0;  ib<128; ib+=2) {
										//unsigned char ddd =  buf_clust[(int)((q-1)*142 + 17 + ib)] + buf_clust[(int)((q-1)*142 + 17 + ib + 1)]*256;
										bmp[ib/2][f] = (double)((buf_clust[(int)((q-1)*142 + 14 + ib)] + buf_clust[(int)((q-1)*142 + 15 + ib)]*256)/4);
										//x = (ddd&0x1);   //  1- trigger  0 - not trigger
									}
									for(int ij = 0; ij < 64; ij = ij + 2) {
										if(e[f][ij] > 0 || sens[f][ij] > 0) {
											if((bmp[ij][f]) >= 3000) {
												bmp[ij][f] = (bmp[ij+1][f] - ped[ij+1][f])/(e[f][ij+1]*sens[f][ij+1]);
												//bmp[ij+1][f] = (bmp[ij+1][f] - ped[ij+1][f])/(e[f][ij+1]*sens[f][ij+1]);
												sig[ij][f] = sig[ij+1][f];
											}
											else {
												bmp[ij][f] = (bmp[ij][f] - ped[ij][f])/(e[f][ij]*sens[f][ij]);
												//bmp[ij+1][f] = (bmp[ij+1][f] - ped[ij+1][f])/(e[f][ij+1]*sens[f][ij+1]);
											}
										}
										else{
											bmp[ij][f] = 0;
											bmp[ij+1][f] = 0;
										}
										//sig[ij][f] = sig[ij][f]/(e[f][ij]*sens[f][ij]);
										//sig[ij+1][f] = sig[ij+1][f]/(e[f][ij+1]*sens[f][ij+1]);
									}
								}
								for(f = 1; f <= 25; f++)
								{
									for (int sc = 0; sc < 64; sc = sc + 2)
									{
										if(bmp[sc][f]>edge1*sig[sc][f])
										{
											if((bmp[pos[0][sc][f]][kkk[0][sc][f]]>edge2*sig[pos[0][sc][f]][kkk[0][sc][f]] && bmp[pos[0][sc][f]][kkk[0][sc][f]]>0) ||
												(bmp[pos[1][sc][f]][kkk[1][sc][f]]>edge2*sig[pos[1][sc][f]][kkk[1][sc][f]] && bmp[pos[1][sc][f]][kkk[1][sc][f]]>0) ||
												(bmp[pos[2][sc][f]][kkk[2][sc][f]]>edge2*sig[pos[2][sc][f]][kkk[2][sc][f]] && bmp[pos[2][sc][f]][kkk[2][sc][f]]>0) ||
												(bmp[pos[3][sc][f]][kkk[3][sc][f]]>edge2*sig[pos[3][sc][f]][kkk[3][sc][f]] && bmp[pos[3][sc][f]][kkk[3][sc][f]]>0) ||
												(bmp[pos[4][sc][f]][kkk[4][sc][f]]>edge2*sig[pos[4][sc][f]][kkk[4][sc][f]] && bmp[pos[4][sc][f]][kkk[4][sc][f]]>0) ||
												(bmp[pos[5][sc][f]][kkk[5][sc][f]]>edge2*sig[pos[5][sc][f]][kkk[5][sc][f]] && bmp[pos[5][sc][f]][kkk[5][sc][f]]>0))
											{
												background_marker[sc][f] = 1;
											}

										}
										else if(bmp[sc][f]>edge2*sig[sc][f])
										{
											if((bmp[pos[0][sc][f]][kkk[0][sc][f]]>edge1*sig[pos[0][sc][f]][kkk[0][sc][f]] && bmp[pos[0][sc][f]][kkk[0][sc][f]]>0) ||
												(bmp[pos[1][sc][f]][kkk[1][sc][f]]>edge1*sig[pos[1][sc][f]][kkk[1][sc][f]] && bmp[pos[1][sc][f]][kkk[1][sc][f]]>0) ||
												(bmp[pos[2][sc][f]][kkk[2][sc][f]]>edge1*sig[pos[2][sc][f]][kkk[2][sc][f]] && bmp[pos[2][sc][f]][kkk[2][sc][f]]>0) ||
												(bmp[pos[3][sc][f]][kkk[3][sc][f]]>edge1*sig[pos[3][sc][f]][kkk[3][sc][f]] && bmp[pos[3][sc][f]][kkk[3][sc][f]]>0) ||
												(bmp[pos[4][sc][f]][kkk[4][sc][f]]>edge1*sig[pos[4][sc][f]][kkk[4][sc][f]] && bmp[pos[4][sc][f]][kkk[4][sc][f]]>0) ||
												(bmp[pos[5][sc][f]][kkk[5][sc][f]]>edge1*sig[pos[5][sc][f]][kkk[5][sc][f]] && bmp[pos[5][sc][f]][kkk[5][sc][f]]>0))
											{
												background_marker[sc][f] = 1;
											}
										}
									}
									
									for (int sc = 0; sc < 64; sc = sc + 2)
									{
										if(background_marker[sc][f] == 0 && bmp[sc][f] != 0)
										{
											if(vector_background[pix_number[sc][f]].size() < 10000 && pix_number[sc][f]>=0)
											{
												vector_background[pix_number[sc][f]].push_back(bmp[sc][f]);
											}
											bmp[sc][f] = 0;
										}
									}
								}
								//cout << "ok" << endl;
								Events event;
								vector<vector<double> > vector_pixel( 5, vector<double> (0));
								for (int count = 0; count < 25; count++)
								{
									for (int coun = 0; coun < 64; coun+=2)
									{
										if(bmp[coun][count] != 0) {
											vector_pixel[0].push_back(count);
											vector_pixel[1].push_back(int(coun/2));
											vector_pixel[4].push_back(bmp[coun][count]);
											vector_pixel[2].push_back(x_pos[coun][count]);
											vector_pixel[3].push_back(y_pos[coun][count]);
										}
									}
								}
								event.set_event(i+1, (int)(buf_clust[1] + buf_clust[2]*256 + buf_clust[3]*256*256 + buf_clust[4]*256*256*256), event_unix_time, nsec_time, vector_pixel);
								//cout << (int)(buf_clust[(int)((q-1)*142 + 1)] + buf_clust[(int)((q-1)*142 + 5)]*256 + buf_clust[(int)((q-1)*142 + 3)]*256*256 + buf_clust[(int)((q-1)*142 + 4)]*256*256*256) << endl;
								//event.ccd_id(vector_ccd);
								if(time_0 + 12. <= event.unix_time && event.number_of_pixels > 3) {
									//cout << time_0 << "\t" << event_unix_time << endl;
									if (event.number_of_pixels > 0) {
										//cout << ccd_id << endl;
										ccd_id = event.get_ccd_parameters(ccd_id, vector_ccd);
										event.star_correction(bmp, kkk, pos);
										event.get_hillas();
										event.to_deg();
										event.get_edge(Nsos_array);
										//cout << event.portion << "\t" << event.number << "\t" << event.number_of_pixels << "\t" << FolderList[jl].c_str() << "." << RunNumbList[jl].c_str() << " " << timet << "\t"  << event.size  << "\t" << event.star << " " << event.edge << endl;
										if(save_background == 1){
											fout << event.number << "\t" << event.number_of_pixels << "\t" << timet << "\t" << event.size << endl;
											//cout << event.number << "\t" << event.number_of_pixels << "\t" << timet << "\t" << event.size << endl;
										}
										vector_events.push_back(event);
									}
									if(save_background == 1){
										for (int count = 0; count < vector_pixel[0].size(); count++) {
											for (int coun = 0; coun < 5; coun++) {
												fout << vector_pixel[coun][count] << "\t";
											}
											fout << endl;
										}
									}
								}
								qqq++;
							}
							}
						}
						fclose(ptrFile);
					}
				}
				fout.close();
				double por_cr = write_cr_file(vector_events);
				for(int count = 0; count < vector_events.size(); count++) {
					//"por, event_numb, unix_time, unix time after dot(ns), delta_time, error_deg, tel_az, tel_el, source_az, source_el, CR5sec, CR_portion, numb_pix, size, Xc[0],Yc[0], con2,
					//length[0], width[0], dist[0], dist[1], dist[2], azwidth[1], azwidth[2], miss[1], miss[2], alpha[0], alpha[1], alpha[2], a_axis, b_axis, a_dist[1], b_dist[1], a_dist[2], b_dist[2],
					//tel_ra,tel_dec,source_ra,source_dec,source_x,source_y,tracking,good,star,edge,weather_mark,alpha_c"
					fout_hillas << fixed << vector_events[count].portion << "," << vector_events[count].number << "," << setprecision(6) << vector_events[count].unix_time << "," << setprecision(19) << vector_events[count].nsec_time << "," <<
					        setprecision(2) << vector_events[count].delta << "," << setprecision(2) << vector_events[count].error_deg << "," << setprecision(5) << vector_events[count].tel_az  << "," <<
					        setprecision(3) << vector_events[count].tel_el << "," << setprecision(5) << vector_events[count].source_az << "," << setprecision(3) << vector_events[count].source_el << "," <<
					        setprecision(2) << vector_events[count].cr_sec << "," << setprecision(2) << por_cr << "," << vector_events[count].number_of_pixels << "," <<
					        setprecision(2) << vector_events[count].size << "," << setprecision(6) << vector_events[count].Xc[0] << "," << setprecision(6) << vector_events[count].Yc[0] << "," <<
					        setprecision(2) << vector_events[count].con2 << "," << setprecision(6) << vector_events[count].length[0] << "," <<
					        vector_events[count].width[0] << "," << vector_events[count].dist[0] << "," << vector_events[count].dist[1] << "," <<
					        vector_events[count].dist[2] << "," << setprecision(6) << vector_events[count].skewness[0] << "," <<
					        vector_events[count].skewness[1] << "," << vector_events[count].skewness[2] << "," << vector_events[count].kurtosis << "," <<
					        setprecision(1) << vector_events[count].alpha[0] << "," << vector_events[count].alpha[1] << "," << vector_events[count].alpha[2] << "," <<
					        setprecision(6) << vector_events[count].a_axis[0] << "," << vector_events[count].b_axis[0] << "," << vector_events[count].a_dist[1] << "," <<
					        vector_events[count].b_dist[1] << "," << vector_events[count].a_dist[2] << "," << vector_events[count].b_dist[2] << "," <<
					        setprecision(2) << vector_events[count].tel_ra << "," << vector_events[count].tel_dec << "," << vector_events[count].source_ra << "," <<
					        vector_events[count].source_dec << "," << setprecision(2) << vector_events[count].source_x << "," << vector_events[count].source_y << "," <<
					        vector_events[count].tracking << "," << vector_events[count].good << "," << vector_events[count].star << "," << vector_events[count].edge << "," << 
					        vector_events[count].weather << "," << vector_events[count].alpha_c << endl;
				}
				if(save_background == 1){
					for(f = 1; f <= 25; f++){
						for (int sc = 0; sc < 64; sc = sc + 2){
							if(pix_number[sc][f]!=-1){
								if(vector_background[pix_number[sc][f]].size()>0){
									fout_background << f << "," << int(sc/2) << ",";
									for(jj = 0; jj <= vector_background[pix_number[sc][f]].size(); jj++){
										fout_background << vector_background[pix_number[sc][f]][jj] << ",";
									}
									fout_background << endl;
								}
							}
						}
					}
				}
				fout_background.close();
						
//por, event_numb, unix_time, delta_time, error_deg, altitude, CR5sec, CR_portion, numb_pix, size, Xc[0],Yc[0], con2, length[0], width[0], dist[0], dist[1], dist[2], azwidth[1], azwidth[2], miss[1], miss[2], alpha[0], alpha[1], alpha[2], source_x, source_y, source_ra, source_dec, tracking, good
				vector_events.clear();
			}
		}
		fout_hillas.close();

	}
	return 0;
}
