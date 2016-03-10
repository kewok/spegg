# In addition to these header files, you will need to add the header file for the species you are simulating, e.g., 
# -I./header/species/my_species
HEADERS = -I./header -I./header/environ -I./header/math -I./header/species -I./header/util -I./header/species/add_kids -I./header/species/movement -I./header/species/update

# Potentially add optimization instructions: e.g., to optimize on SM_3.0 GPU architecture
CFLAGS = -O3 
CFLAGS += -arch=sm_30

# Aliases for directories; you can potentially add more here if you have further directories, e.g.,
# MY_SPECIES = ${SPECIES}/my_species
SRC = ./src
CUDA = ${SRC}/cuda
SPECIES = ${CUDA}/species
MATH = ${CUDA}/math
ENVIRON = ${CUDA}/environ
ADD_KIDS = ${SPECIES}/add_kids
UPDATE = ${SPECIES}/update
MOVEMENT = ${SPECIES}/movement
UTIL = ${CUDA}/util
OBJDIR = objdir

all: $(OBJDIR) a.out

# The location of your cuda libraries. Customize based on your system specific configuration.
export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

# Instructions to create the executable. Note not all the code base is used here.
a.out : $(OBJDIR)/main.o $(OBJDIR)/inds.o $(OBJDIR)/inds_stochastic.o $(OBJDIR)/inds_stochastic_migratory.o $(OBJDIR)/thrust_prob_table.o $(OBJDIR)/ConfigFile.o $(OBJDIR)/thrust_prob_table_demes.o $(OBJDIR)/mating_thrust_prob_table_demes.o $(OBJDIR)/environment.o $(OBJDIR)/amplify.o $(OBJDIR)/random_variables_functions.o $(OBJDIR)/reduce_by_key_with_zeroes.o $(OBJDIR)/statistics_class.o $(OBJDIR)/deme_specific_data_class.o $(OBJDIR)/parents_class.o $(OBJDIR)/neonates_class.o $(OBJDIR)/genetic_deme_specific_data.o $(OBJDIR)/one_dim_two_dim.o $(OBJDIR)/remove_duplicate_pairs.o $(OBJDIR)/which_function.o $(OBJDIR)/gather_values_by_deme.o $(OBJDIR)/histogram.o $(OBJDIR)/Sampling_Input.o $(OBJDIR)/species_specific_mate_sampling_rules.o $(OBJDIR)/Sampling_Event.o $(OBJDIR)/Sample_without_Replacement_1Pass.o $(OBJDIR)/Sample_With_Replacement.o $(OBJDIR)/genotype_phenotype_map_parameters.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/demographic_statistics_class.o $(OBJDIR)/Simulation_Class.o $(OBJDIR)/MigrationBehavior.o 
	nvcc -O3 -lcurand -lrt -lcuda -lconfig++ $(OBJDIR)/main.o $(OBJDIR)/inds.o $(OBJDIR)/inds_stochastic.o $(OBJDIR)/inds_stochastic_migratory.o $(OBJDIR)/thrust_prob_table.o $(OBJDIR)/ConfigFile.o $(OBJDIR)/thrust_prob_table_demes.o $(OBJDIR)/mating_thrust_prob_table_demes.o  $(OBJDIR)/environment.o $(OBJDIR)/amplify.o $(OBJDIR)/random_variables_functions.o $(OBJDIR)/reduce_by_key_with_zeroes.o $(OBJDIR)/statistics_class.o $(OBJDIR)/deme_specific_data_class.o $(OBJDIR)/parents_class.o $(OBJDIR)/neonates_class.o $(OBJDIR)/genetic_deme_specific_data.o $(OBJDIR)/one_dim_two_dim.o $(OBJDIR)/remove_duplicate_pairs.o $(OBJDIR)/which_function.o $(OBJDIR)/gather_values_by_deme.o $(OBJDIR)/histogram.o $(OBJDIR)/Sampling_Input.o $(OBJDIR)/species_specific_mate_sampling_rules.o $(OBJDIR)/Sampling_Event.o $(OBJDIR)/Sample_without_Replacement_1Pass.o $(OBJDIR)/Sample_With_Replacement.o $(OBJDIR)/genotype_phenotype_map_parameters.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/demographic_statistics_class.o $(OBJDIR)/Simulation_Class.o $(OBJDIR)/MigrationBehavior.o 

$(OBJDIR):
	mkdir $(OBJDIR)

# Compilation instructions:

$(OBJDIR)/main.o : main.cpp 
	nvcc -c $(CFLAGS) ${HEADERS} main.cpp -o $(OBJDIR)/main.o 

$(OBJDIR)/ConfigFile.o : ${SRC}/ConfigFile.cpp
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/ConfigFile.cpp -o $(OBJDIR)/ConfigFile.o

$(OBJDIR)/Simulation_Class.o : ${SRC}/Simulation_Class.cpp
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Simulation_Class.cpp -o $(OBJDIR)/Simulation_Class.o
	
$(OBJDIR)/file_util.o : ${UTIL}/file_util.cpp
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/file_util.cpp -o $(OBJDIR)/file_util.o

$(OBJDIR)/UpdateBehavior.o : ${UPDATE}/UpdateBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UPDATE}/UpdateBehavior.cu -o $(OBJDIR)/UpdateBehavior.o

$(OBJDIR)/MigrationBehavior.o : ${MOVEMENT}/MigrationBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MOVEMENT}/MigrationBehavior.cu -o $(OBJDIR)/MigrationBehavior.o

$(OBJDIR)/inds.o : ${SPECIES}/inds.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${SPECIES}/inds.cu -o $(OBJDIR)/inds.o

$(OBJDIR)/inds_stochastic.o : ${SPECIES}/inds_stochastic.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${SPECIES}/inds_stochastic.cu -o $(OBJDIR)/inds_stochastic.o

$(OBJDIR)/inds_stochastic_migratory.o : ${SPECIES}/inds_stochastic_migratory.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${SPECIES}/inds_stochastic_migratory.cu -o $(OBJDIR)/inds_stochastic_migratory.o

$(OBJDIR)/environment.o : ${ENVIRON}/environment.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ENVIRON}/environment.cu -o $(OBJDIR)/environment.o

$(OBJDIR)/mating_thrust_prob_table_demes.o : ${MATH}/mating_thrust_prob_table_demes.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/mating_thrust_prob_table_demes.cu -o $(OBJDIR)/mating_thrust_prob_table_demes.o

$(OBJDIR)/thrust_prob_table_demes.o : ${MATH}/thrust_prob_table_demes.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/thrust_prob_table_demes.cu -o $(OBJDIR)/thrust_prob_table_demes.o
	
$(OBJDIR)/thrust_prob_table.o : ${MATH}/thrust_prob_table.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/thrust_prob_table.cu -o $(OBJDIR)/thrust_prob_table.o
	
$(OBJDIR)/histogram.o :  ${MATH}/histogram.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/histogram.cu -o $(OBJDIR)/histogram.o

$(OBJDIR)/footimer2.o : ${UTIL}/footimer2.cpp
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/footimer2.cpp -o $(OBJDIR)/footimer2.o

$(OBJDIR)/amplify.o : ${UTIL}/amplify.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/amplify.cu -o $(OBJDIR)/amplify.o

$(OBJDIR)/reduce_by_key_with_zeroes.o : ${UTIL}/reduce_by_key_with_zeroes.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/reduce_by_key_with_zeroes.cu -o $(OBJDIR)/reduce_by_key_with_zeroes.o

$(OBJDIR)/random_variables_functions.o : ${MATH}/random_variables_functions.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/random_variables_functions.cu -o $(OBJDIR)/random_variables_functions.o

$(OBJDIR)/statistics_class.o : ${MATH}/statistics_class.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MATH}/statistics_class.cu -o $(OBJDIR)/statistics_class.o

$(OBJDIR)/demographic_statistics_class.o : ${MATH}/demographic_statistics_class.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${MATH}/demographic_statistics_class.cu -o $(OBJDIR)/demographic_statistics_class.o

$(OBJDIR)/deme_specific_data_class.o : ${SPECIES}/deme_specific_data_class.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${SPECIES}/deme_specific_data_class.cu -o $(OBJDIR)/deme_specific_data_class.o

$(OBJDIR)/parents_class.o : ${ADD_KIDS}/parents_class.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/parents_class.cu -o $(OBJDIR)/parents_class.o

$(OBJDIR)/neonates_class.o : ${ADD_KIDS}/neonates_class.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/neonates_class.cu -o $(OBJDIR)/neonates_class.o

$(OBJDIR)/genotype_phenotype_map.o : ${ADD_KIDS}/genotype_phenotype_map.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${ADD_KIDS}/genotype_phenotype_map.cu -o $(OBJDIR)/genotype_phenotype_map.o

$(OBJDIR)/genetic_deme_specific_data.o : ${ADD_KIDS}/genetic_deme_specific_data.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/genetic_deme_specific_data.cu -o $(OBJDIR)/genetic_deme_specific_data.o

$(OBJDIR)/genotype_phenotype_map_parameters.o : ${ADD_KIDS}/genotype_phenotype_map_parameters.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/genotype_phenotype_map_parameters.cu -o $(OBJDIR)/genotype_phenotype_map_parameters.o

$(OBJDIR)/assortative_mating_parents_class.o : ${ADD_KIDS}/assortative_mating_parents_class.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${ADD_KIDS}/assortative_mating_parents_class.cu -o $(OBJDIR)/assortative_mating_parents_class.o

$(OBJDIR)/assortative_mating_neonates_class.o: ${ADD_KIDS}/assortative_mating_neonates_class.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/assortative_mating_neonates_class.cu -o $(OBJDIR)/assortative_mating_neonates_class.o

$(OBJDIR)/one_dim_two_dim.o : ${UTIL}/one_dim_two_dim.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/one_dim_two_dim.cu -o $(OBJDIR)/one_dim_two_dim.o

$(OBJDIR)/remove_duplicate_pairs.o : ${UTIL}/remove_duplicate_pairs.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/remove_duplicate_pairs.cu -o $(OBJDIR)/remove_duplicate_pairs.o

$(OBJDIR)/which_function.o : ${UTIL}/which_function.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/which_function.cu -o $(OBJDIR)/which_function.o

$(OBJDIR)/gather_values_by_deme.o :  ${UTIL}/gather_values_by_deme.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UTIL}/gather_values_by_deme.cu -o $(OBJDIR)/gather_values_by_deme.o

$(OBJDIR)/Sampling_Input.o : ${UTIL}/Sampling_Input.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${UTIL}/Sampling_Input.cu -o $(OBJDIR)/Sampling_Input.o

$(OBJDIR)/Sampling_Event.o : ${UTIL}/Sampling_Event.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${UTIL}/Sampling_Event.cu -o $(OBJDIR)/Sampling_Event.o

$(OBJDIR)/Sample_without_Replacement_1Pass.o : ${UTIL}/Sample_without_Replacement_1Pass.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${UTIL}/Sample_without_Replacement_1Pass.cu -o $(OBJDIR)/Sample_without_Replacement_1Pass.o

$(OBJDIR)/Sample_With_Replacement.o : ${UTIL}/Sample_With_Replacement.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${UTIL}/Sample_With_Replacement.cu -o $(OBJDIR)/Sample_With_Replacement.o

$(OBJDIR)/species_specific_mate_sampling_rules.o : ${ADD_KIDS}/species_specific_mate_sampling_rules.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${ADD_KIDS}/species_specific_mate_sampling_rules.cu -o $(OBJDIR)/species_specific_mate_sampling_rules.o

## Add more compilation instructions for species-specific code, e.g.,

# $(OBJDIR)/mySpecies.o : ${MY_SPECIES}/mySpecies.cu
#	nvcc  -c $(CFLAGS) ${HEADERS}  ${MY_SPECIES}/mySpecies.cu -o $(OBJDIR)/mySpecies.o


# Code related to the clean command.
.PHONY : all clean
clean :
	@rm -rf $(OBJDIR)
