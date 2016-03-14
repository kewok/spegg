#include "coevolutionSimulator.h"
#include <util/footimer2.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{
	footimer2 timer;
	timer.start();

	Coevolution_Simulator *Coevolution_model;
	Coevolution_model = new Coevolution_Simulator();
	Coevolution_model->run();
	delete Coevolution_model;

	timer.stop();
	timer.uprintTime();
	}
