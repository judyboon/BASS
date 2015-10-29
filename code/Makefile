CC=g++
LINKER=-L../Libs/GSL/lib/
INCLUDE=-I../Libs/GSL/include/ -I./
LIBS=-lgsl -lgslcblas -lm
objects = main.o funWrapper.o gig_par.o RngStream.o BASS.o

main : $(objects)
	$(CC) $(objects) -o main $(LINKER) $(LIBS) $(INCLUDE)

main.o : main.cpp BASS.h
	gcc -c main.cpp $(INCLUDE)

funWrapper.o : funWrapper.cpp funWrapper.h
	gcc -c funWrapper.cpp $(INCLUDE)

gig_par.o : gig_par.cpp gig_par.h RngStream.h
	gcc -c gig_par.cpp $(INCLUDE)

RngStream.o : RngStream.cpp RngStream.h
	gcc -c RngStream.cpp $(INCLUDE)

BASS.o : BASS.cpp BASS.h gig_par.h
	gcc -c BASS.cpp $(INCLUDE)

#main.o : BASS.h
#funWrapper.o : funWrapper.h
#gig_par.o : gig_par.h RngStream.h
#RngStream.o : RngStream.h
#BASS.o : BASS.h gig_par.h


clean : 
	rm main *.o
