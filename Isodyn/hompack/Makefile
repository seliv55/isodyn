FFLG = -static -O3 -c

SRC = HOMPACK90.f LAPACK.f

OBJECTS = $(SRC:.f=.o)

CEL = homol.a
all: $(CEL)
$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.f.o:
	gfortran $(FFLG) $< -o $@

clean:
	rm -f *.o *~ $(CEL)

