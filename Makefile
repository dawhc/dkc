SRC_DIR = src
FLAGS = -g -O3 -Wall -fopenmp
SOURCES = $(SRC_DIR)/graph.h $(SRC_DIR)/graph.c $(SRC_DIR)/dkc.h $(SRC_DIR)/dkc.c $(SRC_DIR)/main.c
TARGET = dkc

$(TARGET) : $(SOURCES)
	mpicc -o $@ $(FLAGS) $(SOURCES)
clean :
	rm $(TARGET)

