#include "time_measure.h"

#define MEASURE

static std::chrono::system_clock::time_point init_time;

void initializeTime(){
	init_time = std::chrono::system_clock::now();
}

double elapsedTime(){
	std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
	return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(cur_time - init_time).count() / 1000.0);
}

void writeTime(std::string file_name, std::string func, double t_start, double t_end){
	std::ofstream ofs(file_name, std::ios::app);
	ofs<<std::left<<std::setw(30)<<func<<" : ";
	ofs<<(t_end - t_start)/1000<<" [s]"<<std::endl;
	ofs.close();
}