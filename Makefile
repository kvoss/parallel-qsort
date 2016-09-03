# C-Makefile for Sarek
#
# Assumes that you have run:
# $ module add mpich/pgi
# 
# Note: I'm using $ module add openmpi/gcc
# - K. Voss, July 2009

CC     = mpicc
CFLAGS = -O2 -DTEST
LIBS = 

# Define the application object files and target name
#   APPOBJ = list of required object files
#   APP    = name of target executable

APPOBJ = p-qsort.o
APP = p-qsort

$(APP): $(APPOBJ)
	$(CC) $(CFLAGS) -o $(APP) $(APPOBJ) $(LIBS)

.PHONY: .
clean:
	/bin/rm -f $(APP) $(APPOBJ) 
