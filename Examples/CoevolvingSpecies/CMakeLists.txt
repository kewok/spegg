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
include_directories(${CMAKE_SOURCE_DIR}/../../header)
include_directories(${CMAKE_SOURCE_DIR}/../../header)
include_directories(${CMAKE_SOURCE_DIR}/../../header/environ)
include_directories(${CMAKE_SOURCE_DIR}/../../header/math)
include_directories(${CMAKE_SOURCE_DIR}/../../header/species)
include_directories(${CMAKE_SOURCE_DIR}/../../header/species/add_kids)
include_directories(${CMAKE_SOURCE_DIR}/../../header/species/movement)
include_directories(${CMAKE_SOURCE_DIR}/../../header/species/update)
include_directories(${CMAKE_SOURCE_DIR}/../../header/util)

# Add CUDA include directories explicitly
include_directories(/usr/local/cuda/include)
include_directories(${CUDA_INCLUDE_DIRS})
# Resolve external libraries
# Find the libraries
find_library(LIBCONFIG_LIB config++ REQUIRED)
find_library(CURAND_LIB curand REQUIRED HINTS /usr/local/cuda/lib64)
find_library(CUDA_LIB cuda REQUIRED HINTS /usr/local/cuda/lib64)
find_library(RT_LIB rt REQUIRED)

# Define source directories
set(SRC ${CMAKE_SOURCE_DIR}/src)

# Add library - shared objects
add_library(genotype_phenotype_map STATIC ${SRC}/genotype_phenotype_map.cu)
add_library(UpdateBehavior STATIC ${SRC}/UpdateBehavior.cu)
add_library(coevolutionSimulator STATIC ${SRC}/coevolutionSimulator.cu)
add_library(MigrationBehavior STATIC ${SRC}/MigrationBehavior.cu)
add_library(coevolvingSpecie STATIC ${SRC}/coevolvingSpecie.cu)
add_library(coevolvingSpeciesParents STATIC ${SRC}/coevolvingSpeciesParents.cu)
add_library(update_coevolvingSpecies STATIC ${SRC}/update_coevolvingSpecies.cu)
add_library(migrate_coevolvingSpecies STATIC ${SRC}/migrate_coevolvingSpecies.cu)
add_library(myResource STATIC ${SRC}/myResource.cu)
add_library(myConsumer STATIC ${SRC}/myConsumer.cu)
add_library(myConsumer_genotype_phenotype_maps STATIC ${SRC}/myConsumer_genotype_phenotype_maps.cu)
add_library(myResource_genotype_phenotype_maps STATIC ${SRC}/myResource_genotype_phenotype_maps.cu)
add_library(deme_specific_data_class_coevolvingSpecies STATIC ${SRC}/deme_specific_data_class_coevolvingSpecies.cu)

# Add executables
add_executable(${PROJECT_NAME} main.cpp)

# Enable separable compilation for all CUDA targets
set_target_properties(${PROJECT_NAME} 
        genotype_phenotype_map 
        UpdateBehavior 
        coevolutionSimulator 
        MigrationBehavior 
        coevolvingSpecie 
        coevolvingSpeciesParents 
        update_coevolvingSpecies 
        migrate_coevolvingSpecies 
        myResource 
        myConsumer 
        myConsumer_genotype_phenotype_maps 
        myResource_genotype_phenotype_maps 
        deme_specific_data_class_coevolvingSpecies 
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Find shared object files from the codebase
set(CODEBASE ${CMAKE_SOURCE_DIR}/../../objdir2)
file(GLOB CODEBASE_LIBS "${CODEBASE}/*.a")

# Link target libraries - NOTE: YOU NEED TO LINK LIBRARIES THAT DEPEND ON EACH OTHER!
target_link_libraries(genotype_phenotype_map ${CODEBASE_LIBS})
target_link_libraries(UpdateBehavior ${CODEBASE_LIBS})
target_link_libraries(coevolutionSimulator ${CODEBASE_LIBS})
target_link_libraries(MigrationBehavior ${CODEBASE_LIBS})
target_link_libraries(coevolvingSpecie ${CODEBASE_LIBS})
target_link_libraries(coevolvingSpeciesParents ${CODEBASE_LIBS})
target_link_libraries(update_coevolvingSpecies ${CODEBASE_LIBS})
target_link_libraries(migrate_coevolvingSpecies ${CODEBASE_LIBS})
target_link_libraries(myResource ${CODEBASE_LIBS})
target_link_libraries(myConsumer ${CODEBASE_LIBS})
target_link_libraries(myConsumer_genotype_phenotype_maps ${CODEBASE_LIBS})
target_link_libraries(myResource_genotype_phenotype_maps ${CODEBASE_LIBS})
target_link_libraries(deme_specific_data_class_coevolvingSpecies ${CODEBASE_LIBS})

target_link_libraries(${PROJECT_NAME}
    genotype_phenotype_map 
    UpdateBehavior 
    coevolutionSimulator 
    MigrationBehavior 
    coevolvingSpecie 
    coevolvingSpeciesParents 
    update_coevolvingSpecies 
    migrate_coevolvingSpecies 
    myResource 
    myConsumer 
    myConsumer_genotype_phenotype_maps 
    myResource_genotype_phenotype_maps 
    deme_specific_data_class_coevolvingSpecies 
    ${CODEBASE_LIBS}
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDA_LIB}
    ${RT_LIB}
)

# Add library search path
link_directories(${CODEBASE})

add_custom_target(clean_obj
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${FULL_OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E remove_file ${CMAKE_SOURCE_DIR}/spegg
    COMMENT "Cleaning object directory"
)
