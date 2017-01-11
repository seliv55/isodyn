FFLG = -O3 -c

SRC = daux.f  ddassl.f  dlinpk.f dogleg.f  dpmpar.f  enorm.f  fdjac1.f  hybrd.f hybrd1.f  qform.f  qrfac.f  r1mpyq.f  r1updt.f

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

