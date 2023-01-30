#include <vector>
#include <iterator>
#include <algorithm>


double write_cr_file(vector <Events> &vector_event){
	int step = 30; //sec
	int number_events;
	double min_time = vector_event.front().unix_time;
	double max_time = vector_event.back().unix_time;
	double por_time = max_time - min_time;
	vector <double> vector_number_of_events;
	vector <double> vector_events;
	for(double i_grid = min_time; i_grid <= max_time; i_grid+=step) {
		number_events = 0;
		int id = 0;
		for(int e_time = 0; e_time < vector_event.size(); e_time++) {
			if(vector_event[e_time].unix_time >= i_grid && vector_event[e_time].unix_time < i_grid+step) {
			    if (vector_event[e_time].size > 100){
				number_events++;
			    }
				vector_events.push_back(e_time);
			}
			if(vector_event[e_time].unix_time >= i_grid+step) {
				id = e_time;
				break;
			}
		}
		for(int e_time = 0; e_time < vector_events.size(); e_time++) {
			//cout << number_events << "\t" << step << "\t" << (double)number_events/step << endl;
			vector_event[vector_events[e_time]].set_cr((double)number_events/step);
			//cout << vector_event[vector_events[e_time]].number << "\t" << vector_events[e_time] << "\t" << vector_event[vector_events[e_time]].cr_sec << endl;

		}
		vector_events.clear();
	}
	cout << "CR portion is " << vector_event.size()/por_time << " Hz" << endl;
	return vector_event.size()/por_time;
}
