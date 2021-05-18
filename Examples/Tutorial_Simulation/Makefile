HEADERS = -I./../../header -I./header

CFLAGS = -Xcompiler
CFLAGS += -O3

SRC = ./src

# If necessary, update the location associated with the directory where the objdir folder of the main sPEGG code base resides:
OBJDIR_SPEGG = ./../../objdir

OBJDIR = objdir

all: $(OBJDIR_SPEGG) $(OBJDIR) a.out

export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

a.out : $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/Penguins.o $(OBJDIR)/update_Penguins.o $(OBJDIR)/Penguins_genotype_phenotype_maps.o $(OBJDIR)/Penguin_Parents.o $(OBJDIR)/Penguin_Drift_Simulator.o
	nvcc -Xcompiler -O3 -lcurand -lrt -lcuda -lconfig++ $(OBJDIR_SPEGG)/*.o $(OBJDIR)/main.o $(OBJDIR)/genotype_phenotype_map.o $(OBJDIR)/UpdateBehavior.o $(OBJDIR)/Penguins.o $(OBJDIR)/update_Penguins.o $(OBJDIR)/Penguins_genotype_phenotype_maps.o $(OBJDIR)/Penguin_Parents.o $(OBJDIR)/Penguin_Drift_Simulator.o

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/main.o : main.cpp 
	nvcc -c $(CFLAGS) ${HEADERS} main.cpp -o $(OBJDIR)/main.o 

$(OBJDIR)/UpdateBehavior.o : ${SRC}/UpdateBehavior.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/UpdateBehavior.cu -o $(OBJDIR)/UpdateBehavior.o

$(OBJDIR)/genotype_phenotype_map.o : ${SRC}/genotype_phenotype_map.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/genotype_phenotype_map.cu -o $(OBJDIR)/genotype_phenotype_map.o

$(OBJDIR)/Penguins.o : ${SRC}/Penguins.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/Penguins.cu -o $(OBJDIR)/Penguins.o

$(OBJDIR)/Penguin_Parents.o : ${SRC}/Penguin_Parents.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/Penguin_Parents.cu -o $(OBJDIR)/Penguin_Parents.o

$(OBJDIR)/update_Penguins.o : ${SRC}/update_Penguins.cu
	nvcc  -c $(CFLAGS) ${HEADERS}  ${SRC}/update_Penguins.cu -o $(OBJDIR)/update_Penguins.o

$(OBJDIR)/Penguins_genotype_phenotype_maps.o : ${SRC}/Penguins_genotype_phenotype_maps.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Penguins_genotype_phenotype_maps.cu -o $(OBJDIR)/Penguins_genotype_phenotype_maps.o

$(OBJDIR)/Penguin_Drift_Simulator.o : ${SRC}/Penguin_Drift_Simulator.cu
	nvcc -c $(CFLAGS) ${HEADERS} ${SRC}/Penguin_Drift_Simulator.cu -o $(OBJDIR)/Penguin_Drift_Simulator.o

.PHONY : all clean
clean :
	@rm -rf $(OBJDIR) a.out
