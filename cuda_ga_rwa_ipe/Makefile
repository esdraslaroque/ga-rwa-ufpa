# Makefile 

CC	= nvcc
#CFLAGS	= -O2 -arch=sm_21 -I/home/esdras/NVIDIA_CUDA-8.0_Samples/common/inc/
CFLAGS	= -arch=sm_21 -Wno-deprecated-gpu-targets -O2 -I/home/esdras/NVIDIA_CUDA-8.0_Samples/common/inc/
LIBS	= -lm -lcurand
DEPS	= info.h
OBJ	= main.o nsf.o rwa.o curand.o
#OBJ	= main.o nsf.o curand.o
TARGET	= sim_cudarwa

%.o: %.cu $(DEPS)
	$(CC) -c -dc -o $@ $< $(CFLAGS)

sim_cudarwa: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean: 
	rm -f $(OBJ) $(TARGET)
