CFLAGS = -static -O3 -c -I../include

SRC = lab.cpp integrator.cpp distr.cpp nv.cpp stat.cpp tk.cpp readexp.cpp enzfn.cpp ssc.cpp sst.cpp

OBJECTS = $(SRC:.cpp=.o)

CEL = mod.a
all: $(CEL)
$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.cpp.o:
	g++ $(CFLAGS) $< -o $@


clean:
	rm -f *.o *~ $(CEL)

