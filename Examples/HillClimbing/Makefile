HEADERS = -I./../../header -I./header

CFLAGS = -Xcompiler
CFLAGS += -O3
CFLAGS += -arch=sm_35

SRC = ./src

OBJDIR_SPEGG = ./../../objdir
OBJDIR = objdir

all: $(OBJDIR_SPEGG) $(OBJDIR) a.out

export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

a.out : $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/HillClimbers.o $(OBJDIR)/HillClimberParents.o $(OBJDIR)/update_HillClimbers.o $(OBJDIR)/HillClimbers_Simulator.o $(OBJDIR)/HillClimbers_genotype_phenotype_maps.o 
	nvcc -Xcompiler -O3 -lcurand -lrt -lcuda -lconfig++ $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/HillClimbers.o $(OBJDIR)/HillClimberParents.o $(OBJDIR)/update_HillClimbers.o $(OBJDIR)/HillClimbers_Simulator.o $(OBJDIR)/HillClimbers_genotype_phenotype_maps.o 
$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/main.o : main.cpp 
	nvcc -c $(CFLAGS) ${HEADERS} main.cpp -o $(OBJDIR)/main.o 

$(OBJDIR)/UpdateBehavior.o : ${SRC}/UpdateBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/UpdateBehavior.cu -o $(OBJDIR)/UpdateBehavior.o

$(OBJDIR)/genotype_phenotype_map.o : ${SRC}/genotype_phenotype_map.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/genotype_phenotype_map.cu -o $(OBJDIR)/genotype_phenotype_map.o

$(OBJDIR)/HillClimbers.o : ${SRC}/HillClimbers.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/HillClimbers.cu -o $(OBJDIR)/HillClimbers.o

$(OBJDIR)/HillClimberParents.o : ${SRC}/HillClimberParents.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/HillClimberParents.cu -o $(OBJDIR)/HillClimberParents.o

$(OBJDIR)/update_HillClimbers.o : ${SRC}/update_HillClimbers.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/update_HillClimbers.cu -o $(OBJDIR)/update_HillClimbers.o

$(OBJDIR)/HillClimbers_Simulator.o : ${SRC}/HillClimbers_Simulator.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/HillClimbers_Simulator.cu -o $(OBJDIR)/HillClimbers_Simulator.o

$(OBJDIR)/HillClimbers_genotype_phenotype_maps.o : ${SRC}/HillClimbers_genotype_phenotype_maps.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/HillClimbers_genotype_phenotype_maps.cu -o $(OBJDIR)/HillClimbers_genotype_phenotype_maps.o

.PHONY : all clean
clean :
	@rm -rf $(OBJDIR) a.out
