# Makefile for LEM3D
# Author : Kam-wing John Wong
# Email  : kwjw2@cam.ac.uk
# Date   : 21 Sept 2016

# Tell make that these are phony targets
.PHONY: all clean

# Build all of the executable files
all:
	$(MAKE) -C src

# Clean up the executable files
clean:
	$(MAKE) -C src clean

