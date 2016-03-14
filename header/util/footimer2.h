#ifndef FOOTIMER2_H
#define FOOTIMER2_H

#include <sys/time.h>
#include <cuda.h>
#include <stdio.h>

void memory_diagnostic();

struct footimer2
	{
	public:
		void start();
		void stop();
		void printTime();
		void uprintTime();
	
		int getElapsed();
	private:
		struct timeval start_time;
		struct timeval stop_time;
	};

#endif
