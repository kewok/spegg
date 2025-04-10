# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(spegg LANGUAGES CXX CUDA)

# Find CUDA package
find_package(CUDA REQUIRED)

# Set Required C++ Standard Ensure the appropriate C++ standard is set explicitly
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set C and C++ compilers explicitly
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

# Set CUDA compiler explicitly
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

# Set cuda architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES 75-real)

# Set where all the output would go
set(FULL_OUTPUT "${CMAKE_SOURCE_DIR}/objdir2")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${FULL_OUTPUT}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${FULL_OUTPUT}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${FULL_OUTPUT}/..")

# Add include directories for headers
include_directories(${CMAKE_SOURCE_DIR}/header)
include_directories(${CMAKE_SOURCE_DIR}/header/environ)
include_directories(${CMAKE_SOURCE_DIR}/header/math)
include_directories(${CMAKE_SOURCE_DIR}/header/species)
include_directories(${CMAKE_SOURCE_DIR}/header/species/add_kids)
include_directories(${CMAKE_SOURCE_DIR}/header/species/movement)
include_directories(${CMAKE_SOURCE_DIR}/header/species/update)
include_directories(${CMAKE_SOURCE_DIR}/header/util)

# Add CUDA include directories explicitly
include_directories(/usr/local/cuda/include)
include_directories(${CUDA_INCLUDE_DIRS})

# Define source directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(CUDA ${SRC}/cuda)
set(ENVIRON ${CUDA}/environ)
set(MATH ${CUDA}/math)
set(SPECIES ${CUDA}/species)
set(UTIL ${CUDA}/util)
set(ADD_KIDS ${SPECIES}/add_kids)
set(MOVEMENT ${SPECIES}/movement)
set(UPDATE ${SPECIES}/update)

# Add library - shared objects
add_library(Simulation_Class STATIC ${SRC}/Simulation_Class.cu)
add_library(inds STATIC ${SPECIES}/inds.cu)
add_library(inds_stochastic STATIC ${SPECIES}/inds_stochastic.cu)
add_library(inds_stochastic_migratory STATIC ${SPECIES}/inds_stochastic_migratory.cu)
add_library(environment STATIC ${ENVIRON}/environment.cu)
add_library(mating_thrust_prob_table_demes STATIC ${MATH}/mating_thrust_prob_table_demes.cu)
add_library(thrust_prob_table_demes STATIC ${MATH}/thrust_prob_table_demes.cu)
add_library(thrust_prob_table STATIC ${MATH}/thrust_prob_table.cu)
add_library(histogram STATIC ${MATH}/histogram.cu)
add_library(footimer2 STATIC ${UTIL}/footimer2.cpp)
add_library(amplify STATIC ${UTIL}/amplify.cu)
add_library(reduce_by_key_with_zeroes STATIC ${UTIL}/reduce_by_key_with_zeroes.cu)
add_library(Shuffle_by_Key STATIC ${UTIL}/Shuffle_by_Key.cu)
add_library(random_variables_functions STATIC ${MATH}/random_variables_functions.cu)
add_library(statistics_class STATIC ${MATH}/statistics_class.cu)
add_library(demographic_statistics_class STATIC ${MATH}/demographic_statistics_class.cu)
add_library(deme_specific_data_class STATIC ${SPECIES}/deme_specific_data_class.cu)
add_library(parents_class STATIC ${ADD_KIDS}/parents_class.cu)
add_library(neonates_class STATIC ${ADD_KIDS}/neonates_class.cu)
# add_library(genotype_phenotype_map SHARED ${ADD_KIDS}/genotype_phenotype_map.cu)
add_library(genotype_phenotype_map_parameters STATIC ${ADD_KIDS}/genotype_phenotype_map_parameters.cu)
add_library(genetic_deme_specific_data STATIC ${ADD_KIDS}/genetic_deme_specific_data.cu)
# add_library(prey_fluctuations SHARED ${PREY_VARIABLES}/prey_fluctuations.cu)
# add_library(assortative_mating_parents_class SHARED ${ADD_KIDS}/assortative_mating_parents_class.cu)
# add_library(assortative_mating_neonates_class SHARED ${ADD_KIDS}/assortative_mating_neonates_class.cu)
add_library(determine_mortality STATIC ${UPDATE}/determine_mortality.cu)
add_library(MigrationFunctions STATIC ${MOVEMENT}/MigrationFunctions.cu)
add_library(one_dim_two_dim STATIC ${UTIL}/one_dim_two_dim.cu)
add_library(remove_duplicate_pairs STATIC ${UTIL}/remove_duplicate_pairs.cu)
add_library(which_function STATIC ${UTIL}/which_function.cu)
add_library(gather_values_by_deme STATIC ${UTIL}/gather_values_by_deme.cu)
add_library(Sampling_Input STATIC ${UTIL}/Sampling_Input.cu)
add_library(Sampling_Event STATIC ${UTIL}/Sampling_Event.cu)
add_library(Sample_without_Replacement_1Pass STATIC ${UTIL}/Sample_without_Replacement_1Pass.cu)
add_library(Sample_without_Replacement_byDeme STATIC ${UTIL}/Sample_without_Replacement_byDeme.cu)
add_library(Sample_With_Replacement STATIC ${UTIL}/Sample_With_Replacement.cu)
add_library(Sample_without_Replacement_Test STATIC ${UTIL}/Sample_without_Replacement_Test.cu)
add_library(species_specific_mate_sampling_rules STATIC ${ADD_KIDS}/species_specific_mate_sampling_rules.cu)

# Enable separable compilation for all CUDA targets
set_target_properties(Simulation_Class inds inds_stochastic inds_stochastic_migratory environment mating_thrust_prob_table_demes 
    thrust_prob_table_demes thrust_prob_table histogram amplify reduce_by_key_with_zeroes Shuffle_by_Key random_variables_functions 
    statistics_class demographic_statistics_class deme_specific_data_class parents_class neonates_class genotype_phenotype_map_parameters 
    genetic_deme_specific_data determine_mortality MigrationFunctions one_dim_two_dim remove_duplicate_pairs which_function gather_values_by_deme 
    Sampling_Input Sampling_Event Sample_without_Replacement_1Pass Sample_without_Replacement_byDeme Sample_With_Replacement 
    Sample_without_Replacement_Test species_specific_mate_sampling_rules PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON)

# Resolve external libraries
# Find the libraries
find_library(LIBCONFIG_LIB config++ REQUIRED)
find_library(CURAND_LIB curand REQUIRED HINTS /usr/local/cuda/lib64)
find_library(CUDA_LIB cuda REQUIRED HINTS /usr/local/cuda/lib64)
find_library(RT_LIB rt REQUIRED)

# Store the list of external libraries in a variable
set(EXTERNAL_LIBS
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDA_LIB}
    ${RT_LIB}
)

# Link target libraries
# Link the external libraries to each shared object file
target_link_libraries(Simulation_Class ${EXTERNAL_LIBS})
target_link_libraries(inds ${EXTERNAL_LIBS} deme_specific_data_class environment)
target_link_libraries(inds_stochastic ${EXTERNAL_LIBS})
target_link_libraries(inds_stochastic_migratory ${EXTERNAL_LIBS})
target_link_libraries(environment ${EXTERNAL_LIBS})
target_link_libraries(mating_thrust_prob_table_demes ${EXTERNAL_LIBS})
target_link_libraries(thrust_prob_table_demes ${EXTERNAL_LIBS})
target_link_libraries(thrust_prob_table ${EXTERNAL_LIBS})
target_link_libraries(histogram ${EXTERNAL_LIBS})
target_link_libraries(footimer2 ${EXTERNAL_LIBS})
target_link_libraries(amplify ${EXTERNAL_LIBS})
target_link_libraries(reduce_by_key_with_zeroes ${EXTERNAL_LIBS})
target_link_libraries(Shuffle_by_Key ${EXTERNAL_LIBS})
target_link_libraries(random_variables_functions ${EXTERNAL_LIBS})
target_link_libraries(statistics_class ${EXTERNAL_LIBS} histogram reduce_by_key_with_zeroes gather_values_by_deme)
target_link_libraries(demographic_statistics_class ${EXTERNAL_LIBS})
target_link_libraries(deme_specific_data_class ${EXTERNAL_LIBS} genetic_deme_specific_data)
target_link_libraries(parents_class ${EXTERNAL_LIBS})
target_link_libraries(neonates_class ${EXTERNAL_LIBS} deme_specific_data_class amplify mating_thrust_prob_table_demes)
target_link_libraries(genotype_phenotype_map_parameters ${EXTERNAL_LIBS})
target_link_libraries(genetic_deme_specific_data ${EXTERNAL_LIBS})
target_link_libraries(determine_mortality ${EXTERNAL_LIBS})
target_link_libraries(MigrationFunctions ${EXTERNAL_LIBS})
target_link_libraries(one_dim_two_dim ${EXTERNAL_LIBS})
target_link_libraries(remove_duplicate_pairs ${EXTERNAL_LIBS})
target_link_libraries(which_function ${EXTERNAL_LIBS})
target_link_libraries(gather_values_by_deme ${EXTERNAL_LIBS})
target_link_libraries(Sampling_Input ${EXTERNAL_LIBS})
target_link_libraries(Sampling_Event ${EXTERNAL_LIBS})
target_link_libraries(Sample_without_Replacement_1Pass ${EXTERNAL_LIBS})
target_link_libraries(Sample_without_Replacement_byDeme ${EXTERNAL_LIBS})
target_link_libraries(Sample_With_Replacement ${EXTERNAL_LIBS})
target_link_libraries(Sample_without_Replacement_Test ${EXTERNAL_LIBS})
target_link_libraries(species_specific_mate_sampling_rules ${EXTERNAL_LIBS})

# Set the appropriate compile options for CUDA files
# set_target_properties(a.out PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Add a custom clean target to remove the object directory and executable
# By renaming the target to `clean_custom`, you avoid conflicting with the reserved `clean` target.
add_custom_target(clean_obj
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${FULL_OUTPUT}
    COMMENT "Cleaning object directory"
)
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}