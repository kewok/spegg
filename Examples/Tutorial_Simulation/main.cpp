#include "Penguin_Drift_Simulator.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{
	Penguin_Drift_Simulator *Penguin_model;
	Penguin_model = new Penguin_Drift_Simulator();
	Penguin_model->run();
	delete Penguin_model;
	}
