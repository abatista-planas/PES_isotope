.SUFFIXES: .f .f90 .o

F90 = gfortran # compile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS =   SUBROUTINES_F77.o mathFunc.o coordinateTransf.o helperFunc.o testingFunc.o  main.o

OBJS_Test =   SUBROUTINES_F77.o mathFunc.o coordinateTransf.o helperFunc.o testingFunc.o  testing.o

all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	rm *.o 
	rm *.mod 

develop-test: $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	./ejec_test.x
	./ejec.x

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
