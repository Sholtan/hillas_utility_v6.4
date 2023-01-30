/*event.cpp */
#include <string>
#include <vector>
#define Pi 3.141592653589793238462643

class Events {

public:
vector <int> pixel_cluster;
vector <int> pixel_number;
vector <double> pixel_amp;
vector <double> pixel_x;
vector <double> pixel_y;
int portion;
int number;
double cr_sec;
string nsec_time;
int number_of_pixels;
double unix_time;
double error_deg;
double tel_az;
double tel_el;
double tel_ra;
double tel_dec;
double source_ra;
double source_dec;
double source_az;
double source_el;
double source_x;
double source_y;
double star_x;
double star_y;
int tracking;
int good;
double weather;
double alpha_c;
double delta;
int id;
void set_event(int por, int event_number, double tim, string nst, vector<vector<double> > pixel_parameter){
	portion = por;
	number = event_number;
	unix_time = tim;
	nsec_time = nst;
	for (int i = 0; i < pixel_parameter[0].size(); i++) {
		pixel_cluster.push_back(pixel_parameter[0][i]);
		pixel_number.push_back(pixel_parameter[1][i]);
		pixel_amp.push_back(pixel_parameter[4][i]);
		pixel_x.push_back(pixel_parameter[2][i]);
		pixel_y.push_back(pixel_parameter[3][i]);
	}
	number_of_pixels = pixel_amp.size();
}

int get_ccd_parameters(int ccd_id, vector<vector<double> > ccd_data){
	double ccd_unix_time;
	delta = inf;
	delta1 = inf;
	id = 0;
	for (int i = ccd_id; i < ccd_data[0].size(); i++) {
		delta1 = abs(unix_time - ccd_data[0][i]);
		//cout << delta << "\t" << delta1 << endl;
		//cout << i << endl;
		if(delta < delta1){
		    id = i-1;
		    break;
		}
		delta = delta1;
	}
		//if(ccd_data[0][i] > unix_time && abs(unix_time - ccd_data[0][i]) > delta) {
		//	break;
		//}
	ccd_unix_time = ccd_data[0][id];
	error_deg = ccd_data[1][id];
	tel_ra =    ccd_data[2][id];
	tel_dec =   ccd_data[3][id];
	tel_az =    ccd_data[4][id];
	tel_el =    ccd_data[5][id];
	source_ra = ccd_data[6][id];
	source_dec =ccd_data[7][id];
	source_az = ccd_data[8][id];
	source_el = ccd_data[9][id];
	source_x =  ccd_data[10][id];
	source_y =  ccd_data[11][id];
	star_x =    ccd_data[12][id];
	star_y =    ccd_data[13][id];
	tracking =  ccd_data[14][id];
	good =      ccd_data[15][id];
	weather =   ccd_data[16][id];
	alpha_c =   ccd_data[17][id];
	return id;
}

int star;
void star_correction(double b[64][25], int kkk[6][64][25], int pos[6][64][25]){
	star = 0;
	for (int k = 0; k < pixel_amp.size(); k++) {
		if(sqrt(pow(pixel_x[k] - star_x,2)+pow(pixel_y[k] - star_y,2)) <= 1.5) {
			int sos = 0;
			for(int ss = 0; ss < 6; ss++) {
				if(pos[ss][pixel_number[k]][pixel_cluster[k]] > 0 && kkk[ss][pixel_number[k]][pixel_cluster[k]] > 0 &&
				   b[pos[ss][pixel_number[k]][pixel_cluster[k]]][kkk[ss][pixel_number[k]][pixel_cluster[k]]] > 0) {
					sos++;
				}
			}
			if(sos > 1) {
				pixel_amp[k] = pixel_amp[k]/sos;
				b[pixel_number[k]][pixel_cluster[k]] = pixel_amp[k]/sos;
			}
			else{
				pixel_amp[k] = 0;
				b[pixel_number[k]][pixel_cluster[k]] = 0;
			}
			//cout << pixel_x[k] << "\t" << star_x << "\t" << pixel_y[k] << "\t"  << star_y << endl;
			star = 1;
		}
	}
}
bool edge = 0;
void get_edge(int Ns[64][25]){
	for (int i = 0; i < pixel_number.size(); i++) {
		if(Ns[pixel_number[i]][pixel_cluster[i]] != 6) {
			edge = 1;
			break;
		}
	}
}

double Xc[3];
double Yc[3];
double length[3];
double width[3];
double azwidth[3];
double miss[3];
double dist[3];
double alpha[3];
double con2;
double con1;
double amp_max;
double size;
double skewness[3];
double kurtosis;
double a_axis[3];
double b_axis[3];
double a_dist[3];
double b_dist[3];


void to_deg(){
	Xc[0] = 0.1206*Xc[0];
	Yc[0] = 0.1206*Yc[0];
	length[0] = 0.1206*length[0];
	width[0] = 0.1206*width[0];
	dist[0] = 0.1206*dist[0];
	dist[1] = 0.1206*dist[1];
	dist[2] = 0.1206*dist[2];
	azwidth[0] = 0.1206*azwidth[0];
	azwidth[1] = 0.1206*azwidth[1];
	azwidth[2] =  0.1206*azwidth[2];
	miss[0] =  0.1206*miss[0];
	miss[1] =  0.1206*miss[1];
	miss[2] =  0.1206*miss[2];
	source_x = 0.1206*source_x;
	source_y = 0.1206*source_y;
	star_x = 0.1206*star_x;
	star_y = 0.1206*star_y;
	b_axis[0] = 0.1206*b_axis[0];
	b_axis[1] = 0.1206*b_axis[1];
	b_axis[2] = 0.1206*b_axis[2];
	b_dist[0] = 0.1206*b_dist[0];
	b_dist[1] = 0.1206*b_dist[1];
	b_dist[2] = 0.1206*b_dist[2];
	kurtosis = 0.1206*kurtosis;
};

void set_cr(double countrate){
	cr_sec = countrate;
};
void get_hillas(){
	double Xc2[3];
	double Yc2[3];
	double XYc[3];
	double sigx2[3];
	double sigy2[3];
	double sigxy[3];
	double d[3];
	double z[3];
	double U[3];
	double V[3];
	double bm[3];
	double skew;
	double amp2 = 0;
	double longitudinal1 = 0;
	double longitudinal2 = 0;
	double longitudinal3 = 0;
	double longitudinal4 = 0;
	double length_sig = 0;
	//static double hillas[18]; //size, Xc[0],Yc[0], con2, length[0], width[0], dist[0], dist[1], dist[2], azwidth[0],
	//azwidth[1], azwidth[2], miss[0], miss[1], miss[2], alpha[0], alpha[1], alpha[2]
	skew = -1000;
	size = 0;
	for(int i = 0; i < pixel_amp.size(); i++) {
		size+=pixel_amp[i];
	}

	amp2 = 0;
	amp_max = 0;
	con2 = 0;
	con1 = 0;
	size = 0;
	for(int i = 0; i < 3; i++) {
		Xc[i] = 0;
		Yc[i] = 0;
		a_axis[i] = 0;
		b_axis[i] = 0;
		a_dist[i] = 0;
		b_dist[i] = 0;
		Xc2[i]=0;
		Yc2[i]=0;
		XYc[i]=0;
		sigx2[i] = 0;
		sigy2[i] = 0;
		sigxy[i] = 0;
		d[i]=0;
		z[i]=0;
		length[i]=0;
		width[i] = 0;
		azwidth[i]=0;
		U[i]=0;
		V[i]=0;
		bm[i]=0;
		miss[i]=0;
		dist[i]=0;
		alpha[i]=0;
	}
	for(int k = 0; k < pixel_amp.size(); k++) {
		Xc[0] += (pixel_amp[k]*(pixel_x[k]));
		Yc[0] += (pixel_amp[k]*(pixel_y[k]));
		Xc2[0] += (pixel_amp[k]*pow(pixel_x[k],2));
		Yc2[0] += (pixel_amp[k]*pow(pixel_y[k],2));
		XYc[0] += (pixel_amp[k]*(pixel_x[k])*(pixel_y[k]));

		Xc[1] += (pixel_amp[k]*(pixel_x[k]-source_x));
		Yc[1] += (pixel_amp[k]*(pixel_y[k]-source_y));
		Xc2[1] += (pixel_amp[k]*pow(pixel_x[k]-source_x,2));
		Yc2[1] += (pixel_amp[k]*pow(pixel_y[k]-source_y,2));
		XYc[1] += (pixel_amp[k]*(pixel_x[k]-source_x)*(pixel_y[k]-source_y));

		Xc[2] += (pixel_amp[k]*(pixel_x[k]+source_x));
		Yc[2] += (pixel_amp[k]*(pixel_y[k]+source_y));
		Xc2[2] += (pixel_amp[k]*pow(pixel_x[k]+source_x,2));
		Yc2[2] += (pixel_amp[k]*pow(pixel_y[k]+source_y,2));
		XYc[2] += (pixel_amp[k]*(pixel_x[k]+source_x)*(pixel_y[k]+source_y));

		if(amp_max <= pixel_amp[k]) {
			amp2 = amp_max;
			amp_max = pixel_amp[k];
		}
		else if(amp_max > pixel_amp[k] && pixel_amp[k] > amp2) {
			amp2 = pixel_amp[k];
		}

		size+=pixel_amp[k];
	}

	con2 = (amp_max + amp2)/size;
	con1 = amp_max/size;
	for(int i = 0; i < 3; i++) {

		Xc[i] = Xc[i]/size;
		Yc[i] = Yc[i]/size;
		Xc2[i]=Xc2[i]/size;
		Yc2[i]=Yc2[i]/size;
		XYc[i]=XYc[i]/size;
		sigx2[i] = Xc2[i] - pow(Xc[i],2);
		sigy2[i] = Yc2[i] - pow(Yc[i],2);
		sigxy[i] = XYc[i] - (Xc[i])*(Yc[i]);
		d[i]=sigy2[i]-sigx2[i];
		a_axis[i] = (d[i] + sqrt(pow(d[i],2) + 4*pow(sigxy[i], 2)))/(2*sigxy[i]);
		b_axis[i] = Yc[i] - (a_axis[i]*Xc[i]);
		z[i]=sqrt(pow(d[i],2)+4.*pow(sigxy[i],2));
		length[i]=sqrt((sigx2[i] + sigy2[i] + z[i])/2.);
		if(sigx2[i] + sigy2[i] - z[i] >= 0) {
			width[i]=sqrt((sigx2[i] + sigy2[i] - z[i])/2.);
		}
		else{
			width[i] = 0;
		}
		dist[i] = sqrt(pow(Xc[i],2)+pow(Yc[i],2));
		azwidth[i] = sqrt((pow(Xc[i],2)*Yc2[i]-2.*Xc[i]*Yc[i]*XYc[i]+Xc2[i]*pow(Yc[i],2))/pow(dist[i],2));
		if(z[i] != 0) {
			U[i]=(1. + (d[i]/z[i]));
			V[i]=2. - U[i];
			bm[i] = (0.5*(U[i]*pow(Xc[i],2)+V[i]*pow(Yc[i],2)))-(2.*sigxy[i]*Xc[i]*Yc[i]/z[i]);
		}
		else{
			U[i] = nan("");
			bm[i] = nan("");
			V[i] = nan("");
			miss[i] = nan("");
			alpha[i] = nan("");
		}
		if(bm[i] < 0) {
			bm[i]=nan("");
			miss[i]=nan("");
			alpha[i]=nan("");
		}
		else if(bm[i] >= 0) {
			miss[i]=sqrt(bm[i]);
			alpha[i]=asin(miss[i]/dist[i])*(180./Pi);
		}
	}

	if(length[0] != 0){
		for(int k = 0; k < pixel_amp.size(); k++) {
				longitudinal1 += pixel_amp[k]*pow((pixel_x[k] - Xc[0]) * cos(atan(a_axis[0])) + (pixel_y[k] - Yc[0]) * sin(atan(a_axis[0])),1);
				longitudinal2 += pixel_amp[k]*pow((pixel_x[k] - Xc[0]) * cos(atan(a_axis[0])) + (pixel_y[k] - Yc[0]) * sin(atan(a_axis[0])),2);
				longitudinal3 += pixel_amp[k]*pow((pixel_x[k] - Xc[0]) * cos(atan(a_axis[0])) + (pixel_y[k] - Yc[0]) * sin(atan(a_axis[0])),3);
				longitudinal4 += pixel_amp[k]*pow((pixel_x[k] - Xc[0]) * cos(atan(a_axis[0])) + (pixel_y[k] - Yc[0]) * sin(atan(a_axis[0])),4);
			}
			longitudinal1 = longitudinal1/size;
			longitudinal2 = longitudinal2/size;
			longitudinal3 = longitudinal3/size;
			longitudinal4 = longitudinal4/size;
			length_sig = pow(longitudinal2 - pow(longitudinal1, 2), 0.5);
			skew = (longitudinal3 - 3*longitudinal2*longitudinal1 + 2*pow(longitudinal1, 2))/pow(length_sig,3);
			kurtosis = (longitudinal4 - 4*longitudinal3*longitudinal1 + 6*longitudinal2*pow(longitudinal1, 2) - 3*pow(longitudinal1, 4)/pow(length_sig,4));
			
			if(((0 - Xc[0] < 0 && skew > 0) || (0 - Xc[0] > 0 && skew < 0))){
			    skewness[0] = abs(skew);
			}
			else{
			    skewness[0] = -abs(skew);
			}
			if(((source_x - Xc[0] < 0 && skew > 0) || (source_x - Xc[0] > 0 && skew < 0))){
			    skewness[1] = abs(skew);
			}
			else{
			    skewness[1] = -abs(skew);
			}
			if(((-source_x - Xc[0] < 0 && skew > 0) || (-source_x - Xc[0] > 0 && skew < 0))){
			    skewness[2] = abs(skew);
			}
			else{
			    skewness[2] = -abs(skew);
			}
	}
	else{
		skewness[0] = NAN;
		skewness[1] = NAN;
		skewness[2] = NAN;
		kurtosis = NAN;
	}

	a_dist[0] = (Yc[0])/(Xc[0]);
	b_dist[0] = 0;
	a_dist[1] = (Yc[0]-source_y)/(Xc[0]-source_x);
	b_dist[1] = source_y - a_dist[1]*source_x;
	a_dist[2] = (Yc[0]+source_y)/(Xc[0]+source_x);
	b_dist[2] = -source_y + a_dist[2]*source_x;
}

private:

};
