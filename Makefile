CFLG = -c -O3 -Iinclude -IintegrT

LFLG = -O3 -lgfortran

SRC = main.cpp

LIBS = con512tpl/mod.a svd/svd.a nrused/nr.a integrT/integT.a fort/fort.a hompack/homol.a 

COBJ = $(SRC:.cpp=.o)

CEL = isodyn.out

all: $(CEL)

$(CEL): $(COBJ) $(LIBS)
	g++ $(COBJ) $(LIBS) $(LFLG) -o $(CEL)
#	g++ -O3 -lgfortran main.o svd/svd.a /opt/NAG/fll6a22dfl/lib/libnag_nag.a con512tpl/mod.a nrused/nr.a integrT/integT.a
#	g++ -O3 -lgfortran -m64 main.o svd/svd.a con512tpl/mod.a nrused/nr.a integrT/integT.a -o isodyn.out
#	g++ -O3 -lf52 -lf2c -u MAIN__ main.o svd/svd.a con512tpl/mod.a nrused/nr.a integrT/integT.a -o isodyn.out

.cpp.o:
	g++ $(CFLG) $< -o $@

con512tpl/mod.a:
	make -C con512tpl

svd/svd.a:
	make -C svd

integrT/integT.a:
	make -C integrT

nrused/nr.a:
	make -C nrused

clean:
	rm -f isodyn.out *.o *~
	make clean -C con512tpl
	make clean -C svd
#	make clean -C nrused
#	make clean -C integrT

