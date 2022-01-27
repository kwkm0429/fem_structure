#include <chrono>

#define MEASURE

static std::chrono::system_clock::time_point init_time;

void initializeTime(){
	init_time = std::chrono::system_clock::now();
}

double elapsedTime(){
	std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
	return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(cur_time - init_time).count() / 1000.0);
}