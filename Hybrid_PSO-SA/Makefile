PROGRAM	= HYB
CXX	= g++
CFLAGS	= -std=c++11 -Wall

$(PROGRAM): main.o solution.o particle.o
	$(CXX) -o $(PROGRAM) $(CFLAGS) main.o solution.o particle.o

.PHONY:	clean

clean:
	rm -f *.o $(PROGRAM)
