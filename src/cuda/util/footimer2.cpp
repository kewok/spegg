#include <sys/time.h>
#include <iostream>
#include <util/footimer2.h>

using namespace std;

void memory_diagnostic()
	{
	size_t memFree = 0;
	size_t memTotal = 0;

	CUresult cu_res = cuMemGetInfo(&memFree, &memTotal);
	float fraction_available = ((float) memFree)/((float) memTotal);

	std::cout <<  "Memory Stuff: " << fraction_available << " MemFree: " << memFree << " MemTotal: " << memTotal << std::endl;
	}

void footimer2::start()
	{		
	gettimeofday(&start_time, NULL);
	}

void footimer2::stop()
	{
	gettimeofday(&stop_time, NULL);
	}

void footimer2::printTime()
	{
	time_t start_mseconds = start_time.tv_sec*1000 + start_time.tv_usec/1000;
	time_t stop_mseconds = stop_time.tv_sec*1000 + stop_time.tv_usec/1000;
	time_t diff_mseconds = stop_mseconds - start_mseconds;
	
	time_t seconds = diff_mseconds/1000;
	time_t mseconds = diff_mseconds%1000;
	
	cout << "Time elapsed: " << seconds << "s, " << mseconds << "ms" << endl;
	}

void footimer2::uprintTime()
	{
	time_t start_useconds = start_time.tv_sec*1000000 + start_time.tv_usec;
	time_t stop_useconds = stop_time.tv_sec*1000000 + stop_time.tv_usec;
	time_t diff_useconds = stop_useconds - start_useconds;
	
	time_t seconds = diff_useconds/1000000;
	time_t mseconds = (diff_useconds%1000000)/1000;
	time_t useconds = diff_useconds%1000;
	
	cout << "Time elapsed: " << seconds << "s, " << mseconds << "ms, " << useconds << "us" << endl;
	}

int footimer2::getElapsed()
	{
	time_t start_mseconds = start_time.tv_sec*1000 + start_time.tv_usec/1000;
	time_t stop_mseconds = stop_time.tv_sec*1000 + stop_time.tv_usec/1000;
	time_t diff_mseconds = stop_mseconds - start_mseconds;
	
	return (int)diff_mseconds;
	}

