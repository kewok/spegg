#include <species/update/updatebehavior.h>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

class update_HillClimbers : public UpdateBehavior
	{
	class inds_stochastic *species;
	// Constructor
	public:
		update_HillClimbers(inds_stochastic *species) 
		 	{
			this->species = species;
		
			// Copy the constants 
			this->size = species->size;
			this->Number_of_Demes = species->Num_Demes;
			}

	void update();

	protected:

	// variable names
	int MORTALITY_PHENOTYPE_INDEX;
	int size;
	int Number_of_Demes;
	void survive();
	};
