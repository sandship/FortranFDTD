FC=gfortran
OMP=fopenmp
OPT=O3
DEP=utility.f90 config.f90 load_data.f90 calc.f90 output.f90 main.f90
TARGET=fdtd.out


$(TARGET): $(DEP)
	$(FC) -$(OMP) -$(OPT) $(DEP) -o $(TARGET)