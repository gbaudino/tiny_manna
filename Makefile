# Compiler
CXX ?= g++

# Default Flags
OPTFLAGS ?= -g -O2
CXXFLAGS = $(OPTFLAGS) -Wall -Wextra -std=c++17
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
	rm -f $(TARGET) *.o *.dat $(TARGET).s $(TARGET).ii $(TARGET).gcda $(TARGET).bc $(TARGET).l* $(TARGET).r* $(TARGET).w*

.PHONY: run perf

run: $(TARGET)
	./$(TARGET)

perf: $(TARGET)
	perf record -F 1000 -g -- ./$(TARGET)
	perf report -i perf.data

force: clean $(TARGET)

clang: 
	make $(TARGET) CXX=clang++

perf_file: force
	perf stat -r 4 -o performance.txt --append -e cache-references,cache-misses,instructions,cycles,task-clock ./$(TARGET)