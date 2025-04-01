.SUFFIXES: .f .f90 .o

F90 = gfortran # compile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS =   mod_types.o SUBROUTINES_F77.o math_functions.o coordinateTransf.o helper_functions.o  main.o
OBJS_Test =  mod_types.o SUBROUTINES_F77.o math_functions.o coordinateTransf.o helper_functions.o testingFunc.o  testing.o
OBJS_Par =   mod_types.o SUBROUTINES_F77.o math_functions.o coordinateTransf.o helper_functions.o  parallel_example.o
all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	./ejec.x

par : $(OBJS_Par)
	$(F90) $(OBJS_Par) -o paralel.x -fopenmp
	./paralel.x

test: $(OBJS_Test)
	$(F90) $(OBJS_Test) -o ejec_test.x
	.\ejec_test.x

develop: $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	./ejec.x

build : $(OBJS)
	$(F90) $(OBJS) -o ejec.x

clean:
	del *.exe
	del *.x
	del *.o
	del *.mod

build-run-test: $(OBJS_Test)
	$(F90) $(OBJS_Test) -o ejec_test.x
	./ejec_test.x


run: 
	./ejec.x




$(OBJS) :
.f90.o:
	$(F90) $(CMPLFLG) $<
.f.o:
	$(F90) $(CMPLFLG) $<
