CC := g++

SDIR := src
ODIR := obj
IDIR := inc
TARGET := main

SOURCES := $(wildcard ${SDIR}/*.cpp)
OBJECTS := $(patsubst $(SDIR)/%, $(ODIR)/%, $(SOURCES:.cpp=.o))

np=4
CFLAGS := -std=c++17 -g -Wall -Ofast -fopenmp -DN_THREADS=${np}

INC := -I $(IDIR)

$(TARGET): $(OBJECTS)
	$(CC) $^ -o $(TARGET) -fopenmp

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(ODIR)/*.o $(TARGET)

.PHONY: clean
