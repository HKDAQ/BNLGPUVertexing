
INC	:= -I$(CUDA_ROOT)/include -I.
LIB	:= -L$(CUDA_ROOT)/lib64 -lcudart

# NVCCFLAGS	:= -lineinfo -arch=sm_20 --ptxas-options=-v --use_fast_math
# NVCCFLAGS	:= -lineinfo -arch=sm_50 --ptxas-options=-v --use_fast_math
NVCCFLAGS     := -lineinfo -arch=$(GPU_ARCH_TYPE) --ptxas-options=-v --use_fast_math -Xptxas="-dlcm=ca"

all:	daq_code daq_nhits inspect_gpu

debug: daq_code_debug daq_nhits_debug inspect_gpu_debug

daq_code:	daq_code.cu library_daq.h Makefile
	nvcc daq_code.cu -o daq_code $(INC) $(NVCCFLAGS) $(LIB)

daq_nhits:	daq_nhits.cu library_daq.h Makefile
	nvcc daq_nhits.cu -o daq_nhits $(INC) $(NVCCFLAGS) $(LIB)

inspect_gpu:	inspect_gpu.cu Makefile
	nvcc inspect_gpu.cu -o inspect_gpu $(INC) $(NVCCFLAGS) $(LIB)

daq_code_debug:	daq_code.cu library_daq.h Makefile
	nvcc -g -G daq_code.cu -o daq_code_debug $(INC) $(NVCCFLAGS) $(LIB)

daq_nhits_debug:	daq_nhits.cu library_daq.h Makefile
	nvcc -g -G daq_nhits.cu -o daq_nhits_debug $(INC) $(NVCCFLAGS) $(LIB)

inspect_gpu_debug:	inspect_gpu.cu Makefile
	nvcc -g -G inspect_gpu.cu -o inspect_gpu_debug $(INC) $(NVCCFLAGS) $(LIB)

clean:
	rm -f daq_code daq_nhits inspect_gpu daq_code_debug daq_nhits_debug inspect_gpu_debug
