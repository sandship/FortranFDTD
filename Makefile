FC=gfortran
OMP=-fopenmp
OPT=-O3

DEBUG=-fbacktrace
VALIDATION_MEM=-fbounds-check
LAPACK=-I/usr/local/include -llapack95 -llapack -lblas

TARGET_FDTD=fdtd.out
DEP_FDTD=utility.f90 config.f90 load_data.f90 calc.f90 output.f90 main.f90

$(TARGET_FDTD): $(DEP_FDTD)
	$(FC) $(OMP) $(OPT) $(DEP_FDTD) -o $(TARGET_FDTD) $(DEBUG) $(VALIDATION_MEM)
