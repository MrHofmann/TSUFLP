PROGRAM	= SA
CXX	= g++
CFLAGS	= -std=c++11 -Wall

$(PROGRAM): main.o sa.o solution.o
	$(CXX) -o $(PROGRAM) $(CFLAGS) main.o sa.o solution.o

.PHONY:	clean

clean:
	rm -f *.o $(PROGRAM)
