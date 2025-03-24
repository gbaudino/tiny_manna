# Compiler
CXX ?= clang++

RUNS ?= 5

# Default Flags
OPTFLAGS ?= -O3 -march=native -ffast-math -funroll-loops -flto
CXXFLAGS = $(OPTFLAGS) -Wall -Wextra -std=c++23
CPPFLAGS =
LDFLAGS =

# Binary file
TARGET = tiny_manna

# Files
SOURCES = tiny_manna.cpp
OBJS = $(patsubst %.cpp, %.o, $(SOURCES))

# Rules
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o *.dat

.PHONY: run perf

run: $(TARGET)
	./$(TARGET)

perf: $(TARGET)
	perf stat -r $(RUNS) ./$(TARGET)
