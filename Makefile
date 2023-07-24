.SUFFIXES: .f .f90 .o

F90 = gfortran # compile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS = main.o Int_to_Cart_User.o Int_to_Cart.o  SUBROUTINES_F77.o mathFunc.o coordinateTransf.o helperFunc.o

OBJS_Test = testing.o  Int_to_Cart_User.o Int_to_Cart.o  SUBROUTINES_F77.o mathFunc.o coordinateTransf.o helperFunc.o

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
