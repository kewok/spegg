HEADERS = -I./../../header/ -I./header

CFLAGS = -Xcompiler
CFLAGS += -O3 
CFLAGS += -arch=sm_30
CFLAGS += -lgsl
CFLAGS += -lgslcblas

SRC = ./src

OBJDIR_SPEGG = ./../../objdir
OBJDIR = objdir

all: $(OBJDIR_SPEGG) $(OBJDIR) a.out

export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

a.out : $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/coevolvingSpecie.o $(OBJDIR)/coevolvingSpeciesParents.o $(OBJDIR)/update_coevolvingSpecies.o $(OBJDIR)/MigrationBehavior.o $(OBJDIR)/migrate_coevolvingSpecies.o $(OBJDIR)/coevolutionSimulator.o $(OBJDIR)/myResource.o $(OBJDIR)/myConsumer.o  $(OBJDIR)/myConsumer_genotype_phenotype_maps.o $(OBJDIR)/myResource_genotype_phenotype_maps.o $(OBJDIR)/deme_specific_data_class_coevolvingSpecies.o 
	nvcc -O3 -lcurand -lrt -lcuda -lconfig++ -lgsl -lgslcblas $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/coevolvingSpecie.o $(OBJDIR)/coevolvingSpeciesParents.o $(OBJDIR)/update_coevolvingSpecies.o $(OBJDIR)/MigrationBehavior.o $(OBJDIR)/migrate_coevolvingSpecies.o $(OBJDIR)/coevolutionSimulator.o $(OBJDIR)/myResource.o $(OBJDIR)/myConsumer.o $(OBJDIR)/myConsumer_genotype_phenotype_maps.o $(OBJDIR)/myResource_genotype_phenotype_maps.o $(OBJDIR)/deme_specific_data_class_coevolvingSpecies.o
$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/main.o : main.cpp 
	nvcc -c $(CFLAGS) ${HEADERS} main.cpp -o $(OBJDIR)/main.o 

$(OBJDIR)/UpdateBehavior.o : ${SRC}/UpdateBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/UpdateBehavior.cu -o $(OBJDIR)/UpdateBehavior.o

$(OBJDIR)/MigrationBehavior.o : ${SRC}/MigrationBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/MigrationBehavior.cu -o $(OBJDIR)/MigrationBehavior.o

$(OBJDIR)/genotype_phenotype_map.o : ${SRC}/genotype_phenotype_map.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/genotype_phenotype_map.cu -o $(OBJDIR)/genotype_phenotype_map.o

$(OBJDIR)/coevolvingSpecie.o : ${SRC}/coevolvingSpecie.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/coevolvingSpecie.cu -o $(OBJDIR)/coevolvingSpecie.o

$(OBJDIR)/coevolvingSpeciesParents.o : ${SRC}/coevolvingSpeciesParents.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/coevolvingSpeciesParents.cu -o $(OBJDIR)/coevolvingSpeciesParents.o

$(OBJDIR)/update_coevolvingSpecies.o : ${SRC}/update_coevolvingSpecies.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/update_coevolvingSpecies.cu -o $(OBJDIR)/update_coevolvingSpecies.o

$(OBJDIR)/migrate_coevolvingSpecies.o : ${SRC}/migrate_coevolvingSpecies.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/migrate_coevolvingSpecies.cu -o $(OBJDIR)/migrate_coevolvingSpecies.o

$(OBJDIR)/coevolutionSimulator.o : ${SRC}/coevolutionSimulator.cpp
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/coevolutionSimulator.cpp -o $(OBJDIR)/coevolutionSimulator.o

$(OBJDIR)/myResource.o : ${SRC}/myResource.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/myResource.cu -o $(OBJDIR)/myResource.o

$(OBJDIR)/myConsumer.o : ${SRC}/myConsumer.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/myConsumer.cu -o $(OBJDIR)/myConsumer.o

$(OBJDIR)/myConsumer_genotype_phenotype_maps.o : ${SRC}/myConsumer_genotype_phenotype_maps.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/myConsumer_genotype_phenotype_maps.cu -o $(OBJDIR)/myConsumer_genotype_phenotype_maps.o

$(OBJDIR)/myResource_genotype_phenotype_maps.o : ${SRC}/myResource_genotype_phenotype_maps.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/myResource_genotype_phenotype_maps.cu -o $(OBJDIR)/myResource_genotype_phenotype_maps.o

$(OBJDIR)/deme_specific_data_class_coevolvingSpecies.o : ${SRC}/deme_specific_data_class_coevolvingSpecies.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/deme_specific_data_class_coevolvingSpecies.cu -o $(OBJDIR)/deme_specific_data_class_coevolvingSpecies.o

.PHONY : all clean
clean :
	@rm -rf $(OBJDIR) a.out
