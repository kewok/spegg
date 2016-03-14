#include "HillClimbers_Simulator.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{
	HillClimbers_Simulator *HillClimbers_model;
	HillClimbers_model = new HillClimbers_Simulator();
	HillClimbers_model->run();
	delete HillClimbers_model;
	}
