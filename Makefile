# Makefile for Periodic Orbit Computation in Hamiltonian Systems

CC = gcc
CFLAGS = -O2 -Wall
TARGET = main

# Source files
SRC = main.c FUNCTIONS.c Pr1Ex2.c Pr1Ex3.c Pr1Ex4.c Pr3Ex2.c Pr4Ex2.c Pr4Ex3.c Pr5Ex2.c Pr5Ex4.c Pr5ExFinal.c E2.c
HEADERS = FuncionsQR.h

# Build target
all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Clean executable and object files
clean:
	rm -f $(TARGET) *.o
