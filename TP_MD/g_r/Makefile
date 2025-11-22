calculo_dist_radiosMAKEFILE = Makefile
exe   = main
fcomp = gfortran
flags = -O3
OBJS  = main.o

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(fcomp) -c $(flags) $<

$(exe): $(OBJS) $(calculo_dist_radiosMAKEFILE)
	$(fcomp) $(flags) -o $(exe) $(OBJS)

clean:
	rm ./*.o ./*.mod

main.o: main.f90

