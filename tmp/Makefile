CFLG = -c -O3 -I../include

SRC = readexp.cpp 

COBJ = $(SRC:.cpp=.o)

CEL = setmod.out

all: $(CEL)

$(CEL): $(COBJ)
	g++ $(COBJ) -o $(CEL)

.cpp.o:
	g++ -std=c++11 $(CFLG) $< -o $@

clean:
	rm -f setmod.out *.o *~
	
