CPP = g++
OFLAG = -o 
CPPFLAGS= -std=c++14 -I/usr/include/eigen3 -Wall -Wextra -O2
.SUFFIXES : .o .c++ 
.c++.o : 
	$(CPP) $(CPPFLAGS) -c $<


main:  main.cpp signal.hpp sismograph.hpp
	$(CPP) $(CPPFLAGS) main.cpp $(OFLAG) a.out
run: main
	./a.out

clean:
	rm *.o
