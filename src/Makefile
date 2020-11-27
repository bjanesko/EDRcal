# If using Intel fortran (ifort)
#OPT = -O2 -qopenmp -qopenmp-link=static -threads -qopt-matmul
#FC = ifort

# If using gfortran (ifort) with Openmpi
 OPT = -O3 -ffree-line-length-0 -std=legacy -fopenmp
 FC = mpifort

LIB =
EXE = EDRcal

objects = functions.o define.o IO.o funct3D.o main.o

default : $(objects) 
	$(FC) $(OPT) $(objects) $(LIB) -o $(EXE)
	@echo " ------------------------------------------------------ " 
	@echo "       EDRcal has been built successfully. Enjoy."      
	@echo " ------------------------------------------------------ " 

clean:
	rm -f $(EXE) *.o *.mod


#Define modules that used by other files

define.o : define.f90
	$(FC) $(OPT) -c define.f90

functions.o : functions.f90 define.o
	$(FC) $(OPT) -c functions.f90

modules = define.o functions.o
	
#Others

main.o : main.f90 $(modules)
	$(FC) $(OPT) -c main.f90

IO.o : IO.f90 $(modules)
	$(FC) $(OPT) -c IO.f90

funct3D.o : funct3D.f90 $(modules)
	$(FC) $(OPT) -c funct3D.f90

