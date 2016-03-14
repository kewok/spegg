#include <species/update/updatebehavior.h>
#include <util/amplify.h>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

class update_Penguins : public UpdateBehavior
	{
	class inds_stochastic *species;
	// Constructor
	public:
	    	update_Penguins(inds_stochastic *species) 
	    	 	{
	    		this->species = species;
	    
	    		// Copy the constants 
	    		this->size = species->size;
	    		this->Number_of_Demes = species->Num_Demes;
	    		}
	    
	    	void update();

	protected:
		void update_color();

	    	// variable names
	    	int MORTALITY_PHENOTYPE_INDEX;
	    	int CROWN_COLOR_INDEX;
	    	int size;
	    	int Number_of_Demes;
	    	void survive();
    };


struct crown_color_updater
	{
	float *crown_phenotypes;
	float *target_crown_color;

	crown_color_updater(float* crown_phens, float* targ_crown_color) : crown_phenotypes(crown_phens), target_crown_color(targ_crown_color)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's index
	1: individual's deme
	2: the rate at which crown color decays for the individual's deme.
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_index = thrust::get<0>(t);
		int ind_deme = thrust::get<1>(t);
		float delta_crown_color = thrust::get<2>(t);

		float old_crown_color = crown_phenotypes[ind_index];
		crown_phenotypes[ind_index] = old_crown_color - delta_crown_color*(old_crown_color - target_crown_color[ind_deme]);
		}   
	};
