# Makefile for the K-D tree.
#
# make compile	 ----> Compile the source code.
# make debug	 ----> Show both the output and debug messages during the
# 				       construction of the tree.
# make output    ----> Show the binary tree (and nothing else).
# make time      ----> Produces an executable which prints only the time taken
#				 	   to build the tree.
# make leaks	 ----> Produces an executable which can be used to find find
#					   memory leaks.
# make mpidebug  ----> Enables debugging via gdb (the rank of the debugged
#					   process must be set via the environment variable
#					   MPI_DEBUG_RANK).
# make file      ----> Save the k-d tree as a CSV file in the current directory,
#					   in order of increasing level (left to right).
# make clean	 ----> Clean the folder (removes executable/precompiled headers
#		   	           files)

.DEFAULT_GOAL := all
FLAGS = -std=c++20 -O3 -Wall -Wextra -fopenmp
UTILITIES = main.cpp
COMPILER = g++
OUT = KDTreeBuilder.out

all:
	$(COMPILER) $(FLAGS) $(ADDITIONAL_FLAGS) $(UTILITIES) $(SRC) -o $(OUT)

clean:
	rm *.out
