HEADERS = -I./header

CFLAGS = -Xcompiler
CFLAGS += -O3
CFLAGS += -arch=sm_30

SRC = ./src
CUDA = ./src/cuda
SPECIES = ${CUDA}/species
MATH = ${CUDA}/math
ENVIRON = ${CUDA}/environ
ADD_KIDS = ${SPECIES}/add_kids
UPDATE = ${SPECIES}/update
MOVEMENT = ${SPECIES}/movement
UTIL = ${CUDA}/util

OBJDIR = objdir

all: $(OBJDIR) a.out

export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

a.out : $(OBJDIR)/inds.o $(OBJDIR)/inds_stochastic.o $(OBJDIR)/inds_stochastic_migratory.o $(OBJDIR)/thrust_prob_table.o $(OBJDIR)/ConfigFile.o $(OBJDIR)/thrust_prob_table_demes.o $(OBJDIR)/mating_thrust_prob_table_demes.o $(OBJDIR)/environment.o $(OBJDIR)/amplify.o $(OBJDIR)/random_variables_functions.o $(OBJDIR)/reduce_by_key_with_zeroes.o $(OBJDIR)/statistics_class.o $(OBJDIR)/deme_specific_data_class.o $(OBJDIR)/parents_class.o $(OBJDIR)/neonates_class.o $(OBJDIR)/genetic_deme_specific_data.o $(OBJDIR)/determine_mortality.o $(OBJDIR)/MigrationFunctions.o $(OBJDIR)/one_dim_two_dim.o $(OBJDIR)/remove_duplicate_pairs.o $(OBJDIR)/which_function.o $(OBJDIR)/gather_values_by_deme.o $(OBJDIR)/histogram.o $(OBJDIR)/footimer2.o $(OBJDIR)/Sampling_Input.o $(OBJDIR)/species_specific_mate_sampling_rules.o $(OBJDIR)/Sampling_Event.o $(OBJDIR)/Sample_without_Replacement_1Pass.o $(OBJDIR)/Sample_With_Replacement.o $(OBJDIR)/genotype_phenotype_map_parameters.o  $(OBJDIR)/demographic_statistics_class.o $(OBJDIR)/Simulation_Class.o 

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/ConfigFile.o : ${SRC}/ConfigFile.cpp
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/ConfigFile.cpp -o $(OBJDIR)/ConfigFile.o

$(OBJDIR)/Simulation_Class.o : ${SRC}/Simulation_Class.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Simulation_Class.cu -o $(OBJDIR)/Simulation_Class.o
		
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

$(OBJDIR)/genotype_phenotype_map_parameters.o : ${ADD_KIDS}/genotype_phenotype_map_parameters.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/genotype_phenotype_map_parameters.cu -o $(OBJDIR)/genotype_phenotype_map_parameters.o

$(OBJDIR)/genetic_deme_specific_data.o : ${ADD_KIDS}/genetic_deme_specific_data.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/genetic_deme_specific_data.cu -o $(OBJDIR)/genetic_deme_specific_data.o

$(OBJDIR)/prey_fluctuations.o : ${PREY_VARIABLES}/prey_fluctuations.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${PREY_VARIABLES}/prey_fluctuations.cu -o $(OBJDIR)/prey_fluctuations.o

$(OBJDIR)/assortative_mating_parents_class.o : ${ADD_KIDS}/assortative_mating_parents_class.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${ADD_KIDS}/assortative_mating_parents_class.cu -o $(OBJDIR)/assortative_mating_parents_class.o

$(OBJDIR)/assortative_mating_neonates_class.o: ${ADD_KIDS}/assortative_mating_neonates_class.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${ADD_KIDS}/assortative_mating_neonates_class.cu -o $(OBJDIR)/assortative_mating_neonates_class.o

$(OBJDIR)/determine_mortality.o: ${UPDATE}/determine_mortality.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${UPDATE}/determine_mortality.cu -o $(OBJDIR)/determine_mortality.o

$(OBJDIR)/MigrationFunctions.o: ${MOVEMENT}/MigrationFunctions.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${MOVEMENT}/MigrationFunctions.cu -o $(OBJDIR)/MigrationFunctions.o

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

.PHONY : all clean
clean :
	@rm -rf $(OBJDIR)
