CFLAGS = -c -O3 -I../include -I../integrT

SRC = coord.cpp  solvers.cpp  svd.cpp

OBJECTS = $(SRC:.cpp=.o)

CEL = svd.a

all: $(CEL)

$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.cpp.o:
	g++ $(CFLAGS) $< -o $@


clean:
	rm -f *.o *~ $(CEL)

