NVCC = nvcc
CFLAGS = -std=c++17 -O3 -arch=sm_61 --use_fast_math --extra-device-vectorization -lineinfo
LIBS =

# Targets
TARGET = manna_cuda
SOURCE = manna_cuda.cu

all: $(TARGET)

$(TARGET): $(SOURCE) params.h
	$(NVCC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LIBS)

# Profiling targets
profile: $(TARGET)
	nvprof ./$(TARGET)

profile-detailed: $(TARGET)
	nvprof --print-gpu-trace --print-api-trace ./$(TARGET)

profile-metrics: $(TARGET)
	nvprof --metrics achieved_occupancy,gld_throughput,gst_throughput,branch_efficiency ./$(TARGET)

# Debugging
debug: CFLAGS += -g -G
debug: $(TARGET)

# Clean
clean:
	rm -f $(TARGET)

# Test different block sizes
test-blocks: $(TARGET)
	@echo "Testing different block sizes..."
	@for bs in 128 256 512 1024; do \
		echo "Block size: $$bs"; \
		sed -i "s/block_size(.*)/block_size($$bs)/" $(SOURCE); \
		$(MAKE) $(TARGET) > /dev/null 2>&1; \
		./$(TARGET) | grep "Granos/us"; \
		echo ""; \
	done

.PHONY: all profile profile-detailed profile-metrics debug clean test-blocks