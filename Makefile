# Makefile 

CC	= gcc
CFLAGS	= -std=c99 -Wall
LIBS	= -lm
DEPS	= info.h
OBJ	= main.o nsf.o rwa.o
TARGET	= sim_rwa

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

sim_rwa: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean: 
	rm -f $(OBJ) $(TARGET)
