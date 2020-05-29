CC = g++
OPT = -O3 -DNDEBUG
CPPFLAGS = -Wall -std=c++17 $(OPT)
PTHEIGEN = /usr/local/include/eigen3
PTHKDTR = /usr/local/include
INC = -I../include -I$(PTHEIGEN) -I$(PTHKDTR)

RunAndTumble : RunAndTumble.o
	$(CC) $(CPPFLAGS) -o RunAndTumble RunAndTumble.o
	rm RunAndTumble.o

RunAndTumble.o : RunAndTumble.cpp
	$(CC) $(CPPFLAGS) $(INC) -c $<

clean :
	rm -f RunAndTumble.o RunAndTumble
