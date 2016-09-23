PROGRAM	= PSO
CXX	= g++
CFLAGS	= -std=c++11 -Wall

$(PROGRAM): main.o ps.o solution.o particle.o
	$(CXX) -o $(PROGRAM) $(CFLAGS) main.o ps.o solution.o particle.o

.PHONY:	clean

clean:
	rm -f *.o $(PROGRAM)
