CC := g++

SDIR := src
ODIR := obj
IDIR := inc
TARGET := main

SOURCES := $(wildcard ${SDIR}/*.cpp)
OBJECTS := $(patsubst $(SDIR)/%, $(ODIR)/%, $(SOURCES:.cpp=.o))

CFLAGS := -std=c++14 -Wall

INC := -I $(IDIR)

$(TARGET): $(OBJECTS)
	$(CC) $^ -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	rm -f $(ODIR)/*.o $(TARGET)

.PHONY: clean
