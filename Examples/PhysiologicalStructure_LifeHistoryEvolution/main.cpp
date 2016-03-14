#include "Fish_Simulator.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{
	Fish_Simulator *Fish_model;
	Fish_model = new Fish_Simulator();
	Fish_model->run();
	delete Fish_model;
	}
