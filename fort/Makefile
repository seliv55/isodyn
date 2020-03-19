FFLG = -c -O3

SRC = daux.f  ddassl.f  dlinpk.f

OBJECTS = $(SRC:.f=.o)

CEL = fort.a
all: $(CEL)
$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.f.o:
	gfortran $(FFLG) $< -o $@

clean:
	rm -f *.o *~ $(CEL)

