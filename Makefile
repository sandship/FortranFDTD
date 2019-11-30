FC=gfortran
OMP=fopenmp
OPT=O3
DEP=utility.f90 config.f90 load_data.f90 calc.f90 output.f90 main.f90
TARGET=fdtd.exe
DEBUG=fbacktrace
VALIDATION_MEM=fbounds-check

$(TARGET): $(DEP)
	$(FC) -$(OMP) -$(OPT) $(DEP) -o $(TARGET) -$(DEBUG) -$(VALIDATION_MEM)
