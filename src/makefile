# Makefile for LEM3D
# Author : Kam-wing John Wong
# Email  : kwjw2@cam.ac.uk
# Date   : 12 Nov 2013

BDIR = /home/john/Build
TARGET = $(BDIR)/lem

#Compiler
CC = g++

# Include and library paths
EINC= -I/home/kwjw2 -I/home/john/voro/include/voro++ -I/home/john/Software/voro++-0.4.6/src
ELIB= -L/home/john/voro/lib -L/home/john/Software/voro++-0.4.6/src

# Flags for the compiler
CFLAGS = -std=c++0x -Wall -g -O3 -fopenmp -DNDEBUG -msse3 $(EINC) $(ELIB)

# Libraries used
LIBS = -fopenmp -lvoro++

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, $(BDIR)/%.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h, *.hpp)

$(BDIR)/%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) -Wall $(ELIB) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
