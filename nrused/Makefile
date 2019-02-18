CFLAGS = -c -fPIC -std=c++11  -static -O3 -I../include

SRC = odeint.cpp pythag.cpp pzextr.cpp ludcmp.cpp lubksb.cpp rkqs.cpp rkck.cpp stifbs.cpp stiff.cpp simprsparse.cpp bsstep.cpp mmid.cpp sort2a.cpp sprsin.cpp sprsax.cpp sprstx.cpp asolve.cpp atimes.cpp snrm.cpp linbcg.cpp gaussj.cpp svdcmp.cpp svbksb.cpp \
hqr.cpp elmhes.cpp balanc.cpp \
newt.cpp mnewt.cpp  fdjac.cpp  fmin.cpp  lnsrch.cpp
 
CEL = nr.a

OBJECTS = $(SRC:.cpp=.o)

all: $(CEL)
$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o 

.cpp.o: 
	g++ $(CFLAGS) $< -o $@

clean:	
	rm -f *.o *~ $(CEL)
