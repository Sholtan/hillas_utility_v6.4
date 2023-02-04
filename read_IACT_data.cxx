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
// y:
#include <set>

using namespace std;
char fou[190], save_background_str[190], month_year[70], fou_hillas[190], folder_outs[150],
	background_folder[190], clean_folder[190], param_path_str[150];

string folder;

int bsm, cch, ch, ff, clstr_number, event_number, por, trig, mtrx_neigh_cls_n[6][64][25], mtrx_neigh_chnl_n[6][64][25], matrix_neighbors_n[64][25],
	cluster[25], nn, anti_f[25], jj, jjj, j, x, gg;

int IACT_numb, co, number_of_pixels_cam, pix_number[64][25];

double pedp, sigp, amps_mtrx_s_ev[64][25] /*amplitudes matrix for single event*/, ped[64][25],
	sig[64][25], e[25][64] /*single photoelectron amplitude in ADC codes*/,
	sens[25][64] /*relative sensitivity of pixel*/, gain, k_adc, ecode,
	rel_sens, very_first_event_in_file_unix_time = 1, event_unix_time = 0;
double edge1, edge2;
string str, srr[25], event_time, ped_folder[4] = {"peds", "peds.m3s", "peds.mediana", "peds_median_my"},
							data_path, out_data_path, hillas_table_name, clean_out_name, param_wobble_path, cleaning_type;

double hour, minute, sec, mksec, mlsec, nsec, time0, x_pos[64][25], y_pos[64][25], first_file_first_event_unix_timestamp = 0,
																				   last_file_first_event_unix_timestamp_plus120 = 0, event_delay;

map<int, string> calendar = {{1, "jan"}, {2, "feb"}, {3, "mar"}, {4, "apr"}, {5, "may"}, {6, "jun"}, {7, "jul"}, {8, "aug"}, {9, "sep"}, {10, "oct"}, {11, "nov"}, {12, "dec"}};

int exclud_clust[28] = {0}, exclud_numb[28] = {0}, cleaning = -1, ped_param = -1;

char press;

string nsec_time, factor_file, coord_file;

bool clean_only, background_marker[64][25], save_background;

vector<string> dates_to_process; // date list
vector<string> runs_to_process;	 // run list

// ylines:
// set<int> y_set_times;
// vector<int> y_vector_times;

double time_start_end(string run_date, string path); // returns timestamp of given date and time of first event of file at given path
int parse_files_abs_path(string data_path, string data_folder, string run_numb);
void read_config_file(string config_file_name); // initializes global strings, date list and run list
void read_factors_file(string factor_file_name); // initializes e[25][64], sens[25][64]
void write_neighbors_matrices(string xy_coords_file_name);   // initializes mtrx_neigh_cls_n, mtrx_neigh_chnl_n
void write_amp_ped_abs_path(vector<string> *amp_files_abs_paths, 
	vector<string> *ped_files_abs_paths, int date_run_index); // writes amp_files_abs_paths vector and ped_files_abs_paths vector
void make_out_dir(char **config_file_name, int date_run_index); // creates out directory, copies config file there
void read_pointing_data(vector<vector<double>> *vector_ccd, int date_run_index, 
	vector<string> *amp_files_abs_paths, int number_of_portions_to_process);// fills vector_ccd with "pointing_data_*" files data





int main(int argc, char **argv)
{
	cout << "START" << endl;
	read_config_file(argv[1]);
	read_factors_file(factor_file);
	write_neighbors_matrices(coord_file);

	vector<string> amp_files_abs_paths;
	vector<string> ped_files_abs_paths;

	for (int date_run_index = 0; date_run_index < dates_to_process.size(); date_run_index++) // loop through date_run
	{
		amp_files_abs_paths.clear();
		ped_files_abs_paths.clear();
		write_amp_ped_abs_path(&amp_files_abs_paths, &ped_files_abs_paths, date_run_index);
		make_out_dir(argv, date_run_index);
		int number_of_portions_to_process = amp_files_abs_paths.size();
		int ccd_id = 0;
		vector<vector<double>> vector_ccd;
		read_pointing_data(&vector_ccd,date_run_index, &amp_files_abs_paths,
			number_of_portions_to_process);
		if (save_background == 1)
		{
			sprintf(background_folder, "%s%s", folder_outs, "/background");
			sprintf(clean_folder, "%s%s", folder_outs, "/clean");
			mkdir(background_folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			mkdir(clean_folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}

		sprintf(fou_hillas, "%s/%s.%s_out_hillas_%02.0f_%02.1f%s.csv", folder_outs, dates_to_process[date_run_index].c_str(), runs_to_process[date_run_index].c_str(), edge1, edge2, cleaning_type.c_str());
		ofstream fout_hillas(fou_hillas);
		// write header:
		fout_hillas << "por,event_numb,unix_time,unix_time_long_ns,delta_time,error_deg,tel_az,tel_el,source_az,source_el,CR100phe,CR_portion,numb_pix,size,Xc[0],Yc[0],con2,length[0],width[0],dist[0],dist[1],dist[2],skewness[0],skewness[1],skewness[2],kurtosis,alpha[0],alpha[1],alpha[2],a_axis,b_axis,a_dist[1],b_dist[1],a_dist[2],b_dist[2],tel_ra,tel_dec,source_ra,source_dec,source_x,source_y,tracking,good,star,edge,weather_mark,alpha_c" << endl;
		for (int i = 0; i < number_of_portions_to_process; i++)
		{ // loop over files to process
			cout << "file to process:\t" << i << "\t" << amp_files_abs_paths[i] << endl;
			ifstream DataFileOuts;
			DataFileOuts.open(amp_files_abs_paths[i].c_str());
			if (DataFileOuts.is_open())
			{
				for (int coun = 0; coun < 25; coun++)
				{
					for (int count = 0; count < 64; count++)
					{
						ped[count][coun] = 0;
						sig[count][coun] = 0;
					}
				}
				
				vector<Events> vector_events;
				vector<vector<double>> vector_background(number_of_pixels_cam, vector<double>(0));

				cout << i << "\t" << ped_files_abs_paths[i] << endl;

				sprintf(fou, "%s%s%s%s%02.0f%s%02.1f%s%s%03d%s", folder_outs, "/clean/", dates_to_process[date_run_index].c_str(),
						".cleanout_", edge1, "_", edge2, cleaning_type.c_str(), "_", i + 1, ".txt");
				ofstream fout(fou);

				sprintf(save_background_str, "%s%s%s%s%02.0f%s%02.1f%s%s%03d%s", folder_outs, "/background/",
						dates_to_process[date_run_index].c_str(), ".background_", edge1, "_", edge2,
						cleaning_type.c_str(), "_", i + 1, ".txt");
				ofstream fout_background(save_background_str);
				cout << "output: created files:" << endl
					 << "\t" << fou_hillas << endl;
				if (save_background == 1)
				{
					cout << i << "\t" << fou << endl;
					cout << i << "\t" << save_background_str << endl;
				}
				if (amp_files_abs_paths[i].compare(amp_files_abs_paths[i].size() - 3, 3, ped_files_abs_paths[i], ped_files_abs_paths[i].size() - 3, 3) == 0)
				{ // checking if corresponding ped file number matches out file number
					cout << "por=ped" << endl;
					ifstream DataFilePeds;
					DataFilePeds.open(ped_files_abs_paths[i].c_str());
					if (DataFilePeds.is_open())
					{
						// cout << "peds read" << endl;
						while (!DataFilePeds.eof())
						{
							getline(DataFilePeds, str);
							// cout << str << endl;
							if (DataFilePeds.eof())
							{
								cout << "PED " << i + 1 << " is ended" << endl;
								break;
							}
							istringstream iss(str);
							sigp = -40;
							pedp = -40;
							iss >> ff >> ch >> pedp >> sigp; // ff is cluster number, ch is channel number, pedp is pedestal
							// cout << pedp << "\t" << sigp << endl;
							ped[ch][ff] = pedp;
							if (cleaning == 1)
							{
								// nextlineneedsclarification
								sig[ch][ff] = sigp / (e[ff][ch] * sens[ff][ch]); // e[ff][ch] is single photoelectron amplitude in ADC codes
							}
							else
							{
								sig[ch][ff] = 1;
							}
							// cout << ff << "\t" << ch << "\t" << ped[ch][ff] << "\t" << sig[ch][ff] << endl;
						}
						DataFilePeds.close();
						/*if(i == 0) {
							cout << "Check all parameters and press any key to continue" << endl;
							cin >> press;
						}*/

						int very_first_event_in_file_marker = 0;
						while (!DataFileOuts.eof())
						{  // loop through events in amp file
							int *date;
							getline(DataFileOuts, str);
							if (DataFileOuts.eof())
							{
								cout << "outs " << i + 1 << " is ended" << endl;
								break;
							}
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								
								// here amp_mtrx filling starts
								istringstream iss(str);
								int clasters_in_event_n;
								iss >> clasters_in_event_n;
								// cout << clasters_in_event_n << endl;
								for (int count = 0; count < 25; count++)
								{
									cluster[count] = 0; // needsclarification pre-initialization
									for (int coun = 0; coun < 64; coun++)
									{
										amps_mtrx_s_ev[coun][count] = 0;				// amps from out file, pre-initialization
										background_marker[coun][count] = 0; // pre-initialization, background_marker[i][j] will be set to 1 for image and boundary pixels
									}
								}
								for (int clstr_iter_in_event = 1; clstr_iter_in_event <= clasters_in_event_n; clstr_iter_in_event++)
								{
									getline(DataFileOuts, srr[clstr_iter_in_event]);
									istringstream ist(srr[clstr_iter_in_event]);
									ist >> clstr_number >> event_number >> event_time;
									if (clstr_iter_in_event == 1)
									{
										time_cam t;
										t.set_time(dates_to_process[date_run_index], event_time);
										// cout <<  date[0] << "\t" << date[1] << "\t" << date[2] << "\t" << date[3] << ":" <<  date[4] << ":" << date[5] << "," << date[6] << "." << date[7] << "."<< date[8] << endl;
										event_unix_time = t.get_unix_time(); // unix_time in seconds,  needsclarification

										// cout << setprecision(10) << "\n\n LOOK: " << event_unix_time << "\n\n";   // y
										//  exit(1); // y

										nsec_time = t.full_unix_nsec(); // unix_time in nanoseconds

										// cout <<setprecision(6) << fixed << event_unix_time << "\t" << a->GetSec() << "." << a->GetNanoSec() << endl;
										if (very_first_event_in_file_marker == 0)
										{ // only very first line in the file
											// cout << very_first_event_in_file_marker << "\t" << event_unix_time << endl;
											very_first_event_in_file_unix_time = event_unix_time;
										}
									}
									cluster[clstr_iter_in_event] = clstr_number;
									for (int ch_n = 0; ch_n < 64; ch_n = ch_n + 8)
									{
										getline(DataFileOuts, str);
										istringstream ist(str);
										ist >> amps_mtrx_s_ev[ch_n][clstr_number] >> x >> // every second input is trigger_status, not used
											amps_mtrx_s_ev[ch_n + 1][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 2][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 3][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 4][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 5][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 6][clstr_number] >> x >> 
											amps_mtrx_s_ev[ch_n + 7][clstr_number] >> x;
									}
									for (int ch_n = 0; ch_n < 64; ch_n = ch_n + 2)
									{
										if (e[clstr_number][ch_n] > 0 || sens[clstr_number][ch_n] > 0)
										{
											if ((amps_mtrx_s_ev[ch_n][clstr_number]) >= 3000)
											{
												amps_mtrx_s_ev[ch_n][clstr_number] = (amps_mtrx_s_ev[ch_n + 1][clstr_number] - ped[ch_n + 1][clstr_number]) / (e[clstr_number][ch_n + 1] * sens[clstr_number][ch_n + 1]); // every second is anode channel
												// amps_mtrx_s_ev[ch_n+1][clstr_number] = (amps_mtrx_s_ev[ch_n+1][clstr_number] - ped[ch_n+1][clstr_number])/(e[clstr_number][ch_n+1]*sens[clstr_number][ch_n+1]);
												sig[ch_n][clstr_number] = sig[ch_n + 1][clstr_number];
											}
											else
											{
												amps_mtrx_s_ev[ch_n][clstr_number] = (amps_mtrx_s_ev[ch_n][clstr_number] - ped[ch_n][clstr_number]) / (e[clstr_number][ch_n] * sens[clstr_number][ch_n]);
												// amps_mtrx_s_ev[ch_n+1][clstr_number] = (amps_mtrx_s_ev[ch_n+1][clstr_number] - ped[ch_n+1][clstr_number])/(e[clstr_number][ch_n+1]*sens[clstr_number][ch_n+1]);
											}
										}
										else
										{
											amps_mtrx_s_ev[ch_n][clstr_number] = 0; // if ((e[clstr_number][ch_n] > 0 || sens[clstr_number][ch_n] > 0) == false), channel is set to zero
											amps_mtrx_s_ev[ch_n + 1][clstr_number] = 0;
										}
										// sig[ch_n][clstr_number] = sig[ch_n][clstr_number]/(e[clstr_number][ch_n]*sens[clstr_number][ch_n]);
										// sig[ch_n+1][clstr_number] = sig[ch_n+1][clstr_number]/(e[clstr_number][ch_n+1]*sens[clstr_number][ch_n+1]);
									}
								}
								// here amp_mtrx filling ends 
















								for (int f = 1; f <= 25; f++)
								{
									for (int sc = 0; sc < 64; sc = sc + 2)
									{
										if (amps_mtrx_s_ev[sc][f] > edge1 * sig[sc][f]) // edge1 == 14, edge2 == 70,    sig is sigma_of_channel_ped/(e * sens), sig == 1 when (cleaning == 0)
										{									 // above condition checks for picture threshold

											// amps_mtrx_s_ev is a matrix of amps, amps_mtrx_s_ev[channel_n][cluster_n] gives you amp
											// mtrx_neigh_chnl_n is a matrix of neighbor's channel number
											// mtrx_neigh_chnl_n[0][channel_n][cluster_n] gives channel_n of first neighbor
											// mtrx_neigh_cls_n[0][channel_n][cluster_n] gives cluster_n of first neighbor
											if ((amps_mtrx_s_ev[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] > 
													edge2 * sig[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] > 0))
											{								  // above condition checks if one of the neighbors is above picture boundary threshold
												background_marker[sc][f] = 1; // needsclarification, what does it mean background_marker==1???
											}
										}
										else if (amps_mtrx_s_ev[sc][f] > edge2 * sig[sc][f])
										{ // checks for picture boundary threshold
											if ((amps_mtrx_s_ev[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[0][sc][f]][mtrx_neigh_cls_n[0][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[1][sc][f]][mtrx_neigh_cls_n[1][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[2][sc][f]][mtrx_neigh_cls_n[2][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[3][sc][f]][mtrx_neigh_cls_n[3][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[4][sc][f]][mtrx_neigh_cls_n[4][sc][f]] > 0) ||
												(amps_mtrx_s_ev[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] > 
													edge1 * sig[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] && 
													amps_mtrx_s_ev[mtrx_neigh_chnl_n[5][sc][f]][mtrx_neigh_cls_n[5][sc][f]] > 0))
											{								  // checks neighbors for picture threshold
												background_marker[sc][f] = 1; //  needsclarification, maybe when background_marker remain 0, that channel is background
											}
										}
									}

									for (int sc = 0; sc < 64; sc = sc + 2)
									{
										if (background_marker[sc][f] == 0 && amps_mtrx_s_ev[sc][f] != 0)
										{																					   // vector<vector<double> > vector_background( number_of_pixels_cam, vector<double> (0));
											if (vector_background[pix_number[sc][f]].size() < 10000 && pix_number[sc][f] >= 0) // marktofind needsclarification
											{
												vector_background[pix_number[sc][f]].push_back(amps_mtrx_s_ev[sc][f]); // vector_background will have all background signal, the order of pixels in vector_background is from order of rows from neighbors file "neighbours.IACT01"
											}
											amps_mtrx_s_ev[sc][f] = 0;
										}
									}
								}
								// cout << "ok" << endl;
								Events event;
								vector<vector<double>> vector_pixel(5, vector<double>(0));
								for (int count = 0; count < 25; count++)
								{
									for (int coun = 0; coun < 64; coun += 2)
									{
										if (amps_mtrx_s_ev[coun][count] != 0)
										{
											vector_pixel[0].push_back(count);			   // first vector is cluster_number_
											vector_pixel[1].push_back(int(coun / 2));	   // second vector is channel_number/2
											vector_pixel[4].push_back(amps_mtrx_s_ev[coun][count]);   // third vector is signal amp (already divided by (e*sens))
											vector_pixel[2].push_back(x_pos[coun][count]); // fourth vector is x position of pixel
											vector_pixel[3].push_back(y_pos[coun][count]); // fifth vector is y position of pixel
										}
									}
								}
								// i+1 is file_number, event_number is event_number
								event.set_event(i + 1, event_number, event_unix_time, nsec_time, vector_pixel);
								// event.ccd_id(vector_ccd);

								// yel lines:
								/*if (!(y_set_times.find(very_first_event_in_file_unix_time) != y_set_times.end()))
								{
									y_set_times.insert(very_first_event_in_file_unix_time);
									y_vector_times.push_back(very_first_event_in_file_unix_time);
								}
								cout << "\n\n\n";
								cout << setprecision(10) << "lets look at very_first_event_in_file_unix_time   : " << very_first_event_in_file_unix_time << endl;
								cout << setprecision(10) << "and event.unix_time   : " << event.unix_time << endl;
								cout << setprecision(10) << "very_first_event_in_file_unix_time + 12.          : " << very_first_event_in_file_unix_time + 12 << endl;
								cout << setprecision(10) << "event.unix_time-very_first_event_in_file_unix_time: " << event.unix_time-very_first_event_in_file_unix_time << endl;
								cout << "\n\n\n";*/

								if (very_first_event_in_file_unix_time + 12. <= event.unix_time && event.number_of_pixels > 3)
								{
									// cout << very_first_event_in_file_unix_time << "\t" << event_unix_time << endl;
									if (event.number_of_pixels > 0)
									{
										// cout << ccd_id << endl;
										ccd_id = event.get_ccd_parameters(ccd_id, vector_ccd);
										event.star_correction(amps_mtrx_s_ev, mtrx_neigh_cls_n, mtrx_neigh_chnl_n);
										event.get_hillas();
										event.to_deg();
										event.get_edge(matrix_neighbors_n);
										// cout << event.portion << "\t" << event.number << "\t" << event.number_of_pixels << "\t" << dates_to_process[date_run_index].c_str() << "." << runs_to_process[date_run_index].c_str() << " " << event_time << "\t"  << event.size  << "\t" << event.star << " " << event.edge << endl;
										if (save_background == 1)
										{
											fout << event.number << "\t" << event.number_of_pixels << "\t" << event_time << "\t" << event.size << endl;
										}
										vector_events.push_back(event);
									}
									if (save_background == 1)
									{
										for (int count = 0; count < vector_pixel[0].size(); count++)
										{
											for (int coun = 0; coun < 5; coun++)
											{
												fout << vector_pixel[coun][count] << "\t";
											}
											fout << endl;
										}
									}
								}
								very_first_event_in_file_marker++;
						}
					}
				}
				fout.close();
				double por_cr = write_cr_file(vector_events);
				for (int count = 0; count < vector_events.size(); count++)
				{
					//"por, event_numb, unix_time, unix time after dot(ns), delta_time, error_deg, tel_az, tel_el, source_az, source_el, CR5sec, CR_portion, numb_pix, size, Xc[0],Yc[0], con2,
					// length[0], width[0], dist[0], dist[1], dist[2], azwidth[1], azwidth[2], miss[1], miss[2], alpha[0], alpha[1], alpha[2], a_axis, b_axis, a_dist[1], b_dist[1], a_dist[2], b_dist[2],
					// tel_ra,tel_dec,source_ra,source_dec,source_x,source_y,tracking,good,star,edge,weather_mark,alpha_c"
					fout_hillas << fixed << vector_events[count].portion << "," << vector_events[count].number << "," << setprecision(6) << vector_events[count].unix_time << "," << vector_events[count].nsec_time << "," << setprecision(2) << vector_events[count].delta << "," << setprecision(2) << vector_events[count].error_deg << "," << setprecision(5) << vector_events[count].tel_az << "," << setprecision(3) << vector_events[count].tel_el << "," << setprecision(5) << vector_events[count].source_az << "," << setprecision(3) << vector_events[count].source_el << "," << setprecision(2) << vector_events[count].cr_sec << "," << setprecision(2) << por_cr << "," << vector_events[count].number_of_pixels << "," << setprecision(2) << vector_events[count].size << "," << setprecision(6) << vector_events[count].Xc[0] << "," << setprecision(6) << vector_events[count].Yc[0] << "," << setprecision(2) << vector_events[count].con2 << "," << setprecision(6) << vector_events[count].length[0] << "," << vector_events[count].width[0] << "," << vector_events[count].dist[0] << "," << vector_events[count].dist[1] << "," << vector_events[count].dist[2] << "," << setprecision(6) << vector_events[count].skewness[0] << "," << vector_events[count].skewness[1] << "," << vector_events[count].skewness[2] << "," << vector_events[count].kurtosis << "," << setprecision(1) << vector_events[count].alpha[0] << "," << vector_events[count].alpha[1] << "," << vector_events[count].alpha[2] << "," << setprecision(6) << vector_events[count].a_axis[0] << "," << vector_events[count].b_axis[0] << "," << vector_events[count].a_dist[1] << "," << vector_events[count].b_dist[1] << "," << vector_events[count].a_dist[2] << "," << vector_events[count].b_dist[2] << "," << setprecision(2) << vector_events[count].tel_ra << "," << vector_events[count].tel_dec << "," << vector_events[count].source_ra << "," << vector_events[count].source_dec << "," << setprecision(2) << vector_events[count].source_x << "," << vector_events[count].source_y << "," << vector_events[count].tracking << "," << vector_events[count].good << "," << vector_events[count].star << "," << vector_events[count].edge << "," << vector_events[count].weather << "," << vector_events[count].alpha_c << endl;
				}
				if (save_background == 1)
				{
					for (int f = 1; f <= 25; f++)
					{
						for (int sc = 0; sc < 64; sc = sc + 2)
						{
							if (pix_number[sc][f] != -1)
							{
								if (vector_background[pix_number[sc][f]].size() > 0)
								{
									fout_background << f << "," << int(sc / 2) << ",";
									for (jj = 0; jj <= vector_background[pix_number[sc][f]].size(); jj++)
									{
										fout_background << vector_background[pix_number[sc][f]][jj] << ",";
									}
									fout_background << endl;
								}
							}
						}
					}
				}
				fout_background.close();

				// por, event_numb, unix_time, delta_time, error_deg, altitude, CR5sec, CR_portion, numb_pix, size, Xc[0],Yc[0], con2, length[0], width[0], dist[0], dist[1], dist[2], azwidth[1], azwidth[2], miss[1], miss[2], alpha[0], alpha[1], alpha[2], source_x, source_y, source_ra, source_dec, tracking, good
				vector_events.clear();
			}
		}
		fout_hillas.close();
	}
	// ylines:
	/*cout << "y_vector_times.size():    " << y_vector_times.size() << endl;
	for (auto it = y_vector_times.begin(); it != y_vector_times.end(); ++it)
		cout << ' ' << *it << endl;*/
	return 0;
}

double time_start_end(string run_date, string path)
{
	int x, *date;
	string str, first_event_time;
	ifstream DataFileOuts;
	DataFileOuts.open(path.c_str());
	if (!DataFileOuts.is_open())
	{
		cout << "ERROR out file " << run_date << " is not found" << endl
			 << "check out files availability" << endl;
		exit(0);
	}
	else
	{
		getline(DataFileOuts, str);
		if (DataFileOuts.eof())
		{
			cout << "first out file " << run_date << " is empty" << endl
				 << "check size of out files" << endl;
			exit(0);
		}
		else
		{
			getline(DataFileOuts, str);
			istringstream ist(str);
			// cout << str << endl;
			ist >> x >> x >> first_event_time;
			time_cam t;
			t.set_time(run_date, first_event_time);
			DataFileOuts.close();
			t.print_human_string_data_time();
			// cout << (2000 + date[2]) << "." << date[1] << "." << date[0] << "\t" << date[3] << ":" << date[4] << ":" << date[5] << "," << date[6] << "." << date[7] << "." << date[8] << endl;
			return t.get_unix_time(); // returns unix_time accurate to seconds
		}
	}
	return 0;
}

int parse_files_abs_path(string data_path, string data_folder, string run_numb)
{
	sprintf(BashCommandFolder, "%s%s%s%s%s%s", "readlink -e ", data_path.c_str(), data_folder.c_str(), ".", run_numb.c_str(), "/outs/*out_* > List_outs");
	// parses all *out_* files_abs_paths into List_outs file
	system(BashCommandFolder); // taiga_2020/data2020/DATA_IACT02/2020-21.txt.v3/  /k2/DATA_IACT02/2019-20.txt.v3/
	sprintf(BashCommandFolder, "%s%s%s%s%s%s%s%s", "readlink -e ", data_path.c_str(), data_folder.c_str(), ".", run_numb.c_str(), "/", ped_folder[ped_param].c_str(), "/*ped_* > List_peds");
	// parses all *ped_* files_abs_path from mediana folder into List_peds,readlink -e .../data_folder/030122.01/peds.mediana/*ped_* > List_peds
	system(BashCommandFolder);
	return 0;
}

void read_config_file(string config_file_name)
{
	ifstream pParam;
	pParam.open(config_file_name);
	cout << config_file_name << endl;
	if (!pParam.is_open())
	{
		cout << "param file is not found" << endl; // проверка config файла
		exit(1);
	}
	if (pParam.eof())
	{
		cout << "file param is empty" << endl;
		exit(1);
	}

	string junk_str; // for reading strings that won't be used

	getline(pParam, line);
	istringstream ist4(line);
	ist4 >> junk_str >> junk_str >> IACT_numb; // IACT number 0

	getline(pParam, line);
	istringstream ist6(line);
	ist6 >> factor_file; // factors file name

	getline(pParam, line);
	istringstream ist5(line);
	ist5 >> coord_file; // XY coordinate file

	getline(pParam, line);
	istringstream ist2(line);
	ist2 >> junk_str >> ped_param; // to determine what kind of peds to use

	getline(pParam, line);
	istringstream ist1(line);
	ist1 >> junk_str >> cleaning; // type of cleaning
	if (cleaning == 1)
	{
		cleaning_type = "sig";
	}
	else if (cleaning == 0)
	{
		cleaning_type = "fix";
	}

	getline(pParam, line);
	istringstream ist(line);
	ist >> junk_str >> junk_str >> edge1 >> edge2; // cleaning thresholds	14	7

	getline(pParam, line);
	istringstream ist3(line);
	ist3 >> junk_str >> junk_str;
	for (int i = 0; i < 28; i++)
	{
		ist3 >> exclud_clust[i] >> exclud_numb[i]; // excluded pixels
	}

	getline(pParam, line);
	istringstream ist70(line);
	ist70 >> junk_str >> clean_only; // do or dont use pointing data

	getline(pParam, line);
	istringstream ist700(line);
	ist700 >> junk_str >> junk_str >> junk_str >> junk_str >> junk_str >> save_background; // make clean and background files 1

	getline(pParam, line);
	istringstream ist7(line);
	ist7 >> data_path; // path to data folder

	getline(pParam, line);
	istringstream ist8(line);
	ist8 >> out_data_path; // path to output folder

	getline(pParam, line);
	istringstream ist9(line);
	ist9 >> param_wobble_path; // path to wobble folder

	getline(pParam, line); // skips "date run" line

	string run_numb;
	while (!pParam.eof())
	{
		folder = "";
		run_numb = "";
		getline(pParam, line);
		istringstream ist(line);
		ist >> folder >> run_numb; //  date run
		if (folder.length() == 0)
			break;
		dates_to_process.push_back(folder);	 // list of dates
		runs_to_process.push_back(run_numb); // list of runs
	}
	pParam.close();
}

void read_factors_file(string factor_file_name)
{
	ifstream file0(factor_file); //  opens file with k-adc and 1pe
	if (!file0.is_open())
	{
		cout << "calibration file is not found" << endl;
		exit(1);
	}

	int count_single_e_amp_that_more_than_zero = 0; // how many lines from 'factors_*.txt' was with (ecode > 0)
	if (IACT_numb == 0)
	{
		for (int i = 0; i < 10; i++)
		{ // reads header, first 10 lines
			getline(file0, line);
			if (file0.eof())
			{
				cout << "calibration file is empty" << endl;
				exit(1);
			}
		}
		while (!file0.eof())
		{
			getline(file0, line);
			if (!file0.eof())
			{
				istringstream ist(line);
				ist >> bsm >> cch >> gain >> gain >> ecode >> rel_sens; //  overwrites gain?? => skips gain column in file

				if (ecode > 0)
				{						 // single photoelectron amplitude in ADC codes
					e[bsm][cch] = ecode; // bsm - cluster, cch - channel
					sens[bsm][cch] = rel_sens;
					count_single_e_amp_that_more_than_zero++;
				}
				else
				{
					e[bsm][cch] = 1e9;
					sens[bsm][cch] = -1e9;
				}
			}
		}
	}
	else if (IACT_numb == 1)
	{
		getline(file0, line);
		while (!file0.eof())
		{
			getline(file0, line);
			if (!file0.eof())
			{
				istringstream ist(line);
				ist >> bsm >> cch >> ecode >> rel_sens;
				// cout << "\t" << bsm << "\t" << cch << "\t" << ecode << "\t" << rel_sens << endl;
				e[bsm][cch] = ecode;
				sens[bsm][cch] = rel_sens;
				if (ecode > 0)
				{
					e[bsm][cch] = ecode;
					sens[bsm][cch] = rel_sens;
					count_single_e_amp_that_more_than_zero++;
				}
				else
				{
					e[bsm][cch] = 1e9;
					sens[bsm][cch] = -1e9;
				}
			}
		}
	}
	cout << "file calibration size (28*2*22): " << count_single_e_amp_that_more_than_zero << endl;
	cout << "how many lines from 'factors_*.txt' was with (ecode > 0)" << count_single_e_amp_that_more_than_zero << endl;
	file0.close();
}

void write_neighbors_matrices(string xy_coords_file_name)
{
	int counter_of_pixels;
	ifstream neigh_file(make_neighbours(IACT_numb, xy_coords_file_name, exclud_clust, exclud_numb)); // makes and opens file with neighbours for every pixel
	if (!neigh_file.is_open())
	{
		cout << "file neighbours is not found" << endl;
		exit(1);
	}
	counter_of_pixels = 0; // counter of ... pixels
	for (int cluster_iter = 0; cluster_iter < 25; cluster_iter++)
	{
		for (int channel_iter = 0; channel_iter < 64; channel_iter++)
		{
			pix_number[channel_iter][cluster_iter] = -1;
			for (int neighbor = 0; neighbor < 6; neighbor++)
			{
				mtrx_neigh_cls_n[neighbor][channel_iter][cluster_iter] = -1;
				mtrx_neigh_chnl_n[neighbor][channel_iter][cluster_iter] = -1;
			}
		}
	}
	while (!neigh_file.eof())
	{
		int cluster_number_t = 0, number_of_neighbors = 0, anode_extended_channel_number_t = 0,
			neighbor_cluster_number = 0, neighbor_ext_channel_number = 0;
		double x = 0, y = 0;
		getline(neigh_file, line);
		if (neigh_file.eof())
		{
			break;
		}
		istringstream ist(line);
		ist >> cluster_number_t >> anode_extended_channel_number_t;
		pix_number[anode_extended_channel_number_t][cluster_number_t] = counter_of_pixels; // anode_extended_channel_number_t here is only even numbers, => odd elements remain == -1
		counter_of_pixels++;
		ist >> x_pos[anode_extended_channel_number_t][cluster_number_t] >> y_pos[anode_extended_channel_number_t][cluster_number_t] >> number_of_neighbors;
		matrix_neighbors_n[anode_extended_channel_number_t][cluster_number_t] = number_of_neighbors;
		for (int neighbor = 0; neighbor < 6; neighbor++)
		{
			neighbor_cluster_number = 0;
			neighbor_ext_channel_number = 0;
			ist >> neighbor_cluster_number >> neighbor_ext_channel_number;
			mtrx_neigh_cls_n[neighbor][anode_extended_channel_number_t][cluster_number_t] = neighbor_cluster_number;
			mtrx_neigh_chnl_n[neighbor][anode_extended_channel_number_t][cluster_number_t] = neighbor_ext_channel_number;
		}
	}
	number_of_pixels_cam = counter_of_pixels; // counted all lines in "neighbours.IACT01", i.e. pixels with neighbors
	cout << "number of using camera pixels: " << number_of_pixels_cam << endl;
	neigh_file.close();
}

void write_amp_ped_abs_path(vector<string> *amp_files_abs_paths, vector<string> *ped_files_abs_paths, int date_run_index)
{
	cout << "open run: " << dates_to_process[date_run_index] << "." << runs_to_process[date_run_index] << endl;
	ifstream fFileList;
	parse_files_abs_path(data_path, dates_to_process[date_run_index], runs_to_process[date_run_index]); // make "List_outs" file and "List_peds" file (files with parsed abs_path of amp files)
	fFileList.open("List_outs");
	if (!fFileList.is_open())
	{
		cout << "List_outs is not found" << endl;
		exit(1);
	}

	if (fFileList.is_open())
	{
		while (!fFileList.eof())
		{
			string path;
			fFileList >> path;
			if (fFileList.eof() || path.length() == 0)
			{
				cout << "List outs size: " << amp_files_abs_paths->size() << endl;
				break;
			}
			amp_files_abs_paths->push_back(path);
		}
	}
	fFileList.close();

	////////////////////////////////// create Peds vector
	fFileList.open("List_peds");
	if (!fFileList.is_open())
	{
		cout << "List_peds is not found" << endl;
		exit(1);
	}

	if (fFileList.is_open())
	{
		while (!fFileList.eof())
		{
			string path;
			fFileList >> path;
			if (fFileList.eof() || path.length() == 0)
			{
				cout << "List peds size: " << ped_files_abs_paths->size() << endl;
				break;
			}
			ped_files_abs_paths->push_back(path);
		}
	}
	fFileList.close();
}

void make_out_dir(char **config_file_name, int date_run_index)
{
	int dir_err = -1, k = 0; // k - out_directory_version, dir_err - flag for mkdir() function
	while (dir_err == -1)
	{ // program doesn't overwrites out_files, opens available 030122.01(k) directory
		if (k == 0)
		{
			sprintf(folder_outs, "%s%s.%s", out_data_path.c_str(), dates_to_process[date_run_index].c_str(), runs_to_process[date_run_index].c_str());
		}
		else
		{
			sprintf(folder_outs, "%s%s.%s(%d)", out_data_path.c_str(), dates_to_process[date_run_index].c_str(), runs_to_process[date_run_index].c_str(), k);
		}
		dir_err = mkdir(folder_outs, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // 2d argument is permission mode
		k++;
	}
	cout << "out directory " << folder_outs << " is created" << endl;
	////////////////////////copy file param to directiry:
	ifstream src(config_file_name[1], ios::binary);
	sprintf(param_path_str, "%s/%s", folder_outs, config_file_name[1]);
	ofstream dst(param_path_str, ios::binary);
	dst << src.rdbuf(); // copies src(parram_file) contents to dst(param_file at out directory)

	dst << src.rdbuf(); // ???
}

void read_pointing_data(vector<vector<double>> *vector_ccd, int date_run_index, vector<string> *amp_files_abs_paths,int number_of_portions_to_process)
{
	// wobble start
	first_file_first_event_unix_timestamp = time_start_end(dates_to_process[date_run_index], (*amp_files_abs_paths)[0]); // ,
	last_file_first_event_unix_timestamp_plus120 = 120 + time_start_end(dates_to_process[date_run_index], (*amp_files_abs_paths)[number_of_portions_to_process - 1]);
	// needsclarification: why +120 seconds plus first event time, if portion's duration is not 120 seconds???
	// sprintf(month_year, "%s%d", calendar[stoi(dates_to_process[date_run_index].substr(2,2))].c_str(), stoi(dates_to_process[date_run_index].substr(4,2)));
	// cout << month_year << endl;
	vector<string> vector_wobble = wobble(param_wobble_path.c_str(),					 // parses "pointing_data_*" files
										  first_file_first_event_unix_timestamp,		 // whose time is not 12 hours less or more
										  last_file_first_event_unix_timestamp_plus120); // from first_file_first_event_unix_timestamp
	cout << "input:" << endl;
	for (int i = 0; i < vector_wobble.size(); i++)
	{
		cout << i << "\t" << vector_wobble[i] << endl; // pointing_data files that will be used...
	}
	////////////////////////////////////////////////////////////////////////
	*vector_ccd = read_ccd(vector_wobble,
						  first_file_first_event_unix_timestamp, last_file_first_event_unix_timestamp_plus120,
						  clean_only);
	// concatenates pointing_data files whose time is between
	// first_file_first_event_unix_timestamp and last_file_first_event_unix_timestamp_plus120
	// vector_ccd[column][row]

	
	cout << "\t\tnumber of written ccd rows: " << vector_ccd[0].size() << endl;
	// wobble finish
}