# In addition to these header files, you will need to add the header file for the species you are simulating, e.g., 
# -I/.header/species/my_species
HEADERS = -I./../../header/ -I./header 

CFLAGS = -O3 
CFLAGS += -arch=sm_30

SRC = ./src

OBJDIR_SPEGG = ./../../objdir
OBJDIR = objdir

all: $(OBJDIR) a.out

export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

a.out : $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/Fish_Habitat.o $(OBJDIR)/Fish_genotype_phenotype_maps.o $(OBJDIR)/Fish_Simulator.o $(OBJDIR)/Fish.o $(OBJDIR)/prey_class.o $(OBJDIR)/prey_deme_specific_data.o $(OBJDIR)/update_Fish.o $(OBJDIR)/Fish_Parents.o
	nvcc -O3 -lcurand -lrt -lcuda -lconfig++ $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/Fish_Habitat.o $(OBJDIR)/Fish_genotype_phenotype_maps.o $(OBJDIR)/Fish_Simulator.o $(OBJDIR)/Fish.o $(OBJDIR)/prey_class.o $(OBJDIR)/prey_deme_specific_data.o $(OBJDIR)/update_Fish.o $(OBJDIR)/Fish_Parents.o

$(OBJDIR):
	mkdir $(OBJDIR)

# Compilation instructions:

$(OBJDIR)/main.o : main.cpp 
	nvcc -c $(CFLAGS) ${HEADERS} main.cpp -o $(OBJDIR)/main.o

$(OBJDIR)/UpdateBehavior.o : ${SRC}/UpdateBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${SRC}/UpdateBehavior.cu -o $(OBJDIR)/UpdateBehavior.o

$(OBJDIR)/genotype_phenotype_map.o : ${SRC}/genotype_phenotype_map.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/genotype_phenotype_map.cu -o $(OBJDIR)/genotype_phenotype_map.o

$(OBJDIR)/Fish_Habitat.o : ${SRC}/Fish_Habitat.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Fish_Habitat.cu -o $(OBJDIR)/Fish_Habitat.o

$(OBJDIR)/prey_deme_specific_data.o : ${SRC}/prey_deme_specific_data.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/prey_deme_specific_data.cu -o $(OBJDIR)/prey_deme_specific_data.o

$(OBJDIR)/prey_class.o : ${SRC}/prey_class.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/prey_class.cu -o $(OBJDIR)/prey_class.o

$(OBJDIR)/Fish_genotype_phenotype_maps.o : ${SRC}/Fish_genotype_phenotype_maps.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Fish_genotype_phenotype_maps.cu -o $(OBJDIR)/Fish_genotype_phenotype_maps.o

$(OBJDIR)/Fish.o : ${SRC}/Fish.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Fish.cu -o $(OBJDIR)/Fish.o

$(OBJDIR)/update_Fish.o : ${SRC}/update_Fish.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/update_Fish.cu -o $(OBJDIR)/update_Fish.o

$(OBJDIR)/Fish_Parents.o : ${SRC}/Fish_Parents.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Fish_Parents.cu -o $(OBJDIR)/Fish_Parents.o

$(OBJDIR)/Fish_Simulator.o : ${SRC}/Fish_Simulator.cpp
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Fish_Simulator.cpp -o $(OBJDIR)/Fish_Simulator.o

# Code related to the clean command.
.PHONY : all clean
clean :
	@rm -rf $(OBJDIR) a.out
