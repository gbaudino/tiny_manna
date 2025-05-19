# Compiler
CXX ?= clang++

# Default Flags
OPTFLAGS ?= -O3 -march=native -mavx2 -g -fno-omit-frame-pointer -fopenmp
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
	rm -f $(TARGET) *.o $(TARGET).s $(TARGET).ii $(TARGET).gcda $(TARGET).bc $(TARGET).l* $(TARGET).r* $(TARGET).w*

clean_old:
	rm -f old_$(TARGET) *.o *.dat old_$(TARGET).s old_$(TARGET).ii old_$(TARGET).gcda old_$(TARGET).bc old_$(TARGET).l* old_$(TARGET).r* old_$(TARGET).w*

.PHONY: run perf force clang

run: $(TARGET)
	./$(TARGET)

run_old: old_$(TARGET)
	./old_$(TARGET)

force: clean $(TARGET)

perf: force
	perf record -F 1000 -g -- ./$(TARGET)
	perf report -i perf.data

clang: clean
	make $(TARGET) CXX=clang++
	./$(TARGET)

perf_file: force
	perf stat -r 3 -o performance.txt --append -e cache-references,cache-misses,instructions,cycles,task-clock ./$(TARGET)

max:
	grep -o 'max=[0-9]\+' sand.dat | cut -d= -f2 | sort -n | tail -1
