#include "HillClimbers_Simulator.h"
#include <util/footimer2.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{	
	footimer2 timer;
	timer.start();

	HillClimbers_Simulator *HillClimbers_model;
	HillClimbers_model = new HillClimbers_Simulator();
	HillClimbers_model->run();
	delete HillClimbers_model;

	timer.stop();
	timer.uprintTime();
	}
