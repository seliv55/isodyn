CPP = g++

CFLAGS = -c -fPIC -std=c++11  -static -O3

LFLAGS = -lm

DEPENDPATH = -I. -I../include

objects = IntegratorT.o \
	StiffIntegratorT.o \
	NonStiffIntegratorT.o \
	decsol.o

integT.a: $(objects)
	rm -f integT.a
	ar r integT.a *.o

IntegratorT.o: IntegratorT.cpp IntegratorT.h StiffIntegratorT.o \
				NonStiffIntegratorT.o
	$(CPP) $(CFLAGS) $(DEPENDPATH) IntegratorT.cpp
	
StiffIntegratorT.o: StiffIntegratorT.cpp StiffIntegratorT.h decsol.o
	$(CPP) $(CFLAGS) $(DEPENDPATH) StiffIntegratorT.cpp

NonStiffIntegratorT.o: NonStiffIntegratorT.cpp NonStiffIntegratorT.h
	$(CPP) $(CFLAGS) $(DEPENDPATH) NonStiffIntegratorT.cpp

decsol.o: decsol.cpp decsol.h
	$(CPP) $(CFLAGS) $(DEPENDPATH) decsol.cpp

clean:
	rm -f *.o *.a *~
