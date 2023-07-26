include ~/example.mk

CC=mpic++

LDIR =

OBJ1 = anaEul.o
OBJ2 = anaHyb.o
OBJ3 = annulus.o

%.o: %.cpp
	$(CC) -O0 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

aEul: $(OBJ1)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
aHyb: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
Annu: $(OBJ3)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: aEul aHyb Annu

run: all
	mpirun -np 4 /Active3d 

.PHONY: clean all run

clean:
	rm -f *.o *~ core aEul aHyb Annu
