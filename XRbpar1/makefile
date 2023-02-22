
OBJS = main.o input.o star1.o star2.o diskgeom.o diskflux.o lightcurve.o output.o utility.o fitdata.o diagnose.o

XRbinary :  $(OBJS)
	gcc -g -o a.out $(OBJS) -lm

main.o :  main.c header.h
	gcc -c main.c

input.o :  input.c header.h
	gcc -c input.c

star1.o :  star1.c header.h
	gcc -c star1.c

star2.o :  star2.c header.h
	gcc -c star2.c

diskgeom.o :  diskgeom.c header.h
	gcc -c diskgeom.c

diskflux.o :  diskflux.c header.h
	gcc -c diskflux.c

lightcurve.o :  lightcurve.c header.h
	gcc -c lightcurve.c

output.o :  output.c header.h
	gcc -c output.c

utility.o :  utility.c header.h
	gcc -c utility.c

fitdata.o :  fitdata.c header.h
	gcc -c fitdata.c

diagnose.o : diagnose.c header.h
	gcc -c diagnose.c
