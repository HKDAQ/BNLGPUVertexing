//
// include files
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <helper_cuda.h>
#include <sys/time.h>
#include <library_daq.h>

// CUDA = Computer Device Unified Architecture

__global__ void kernel_correct_times(unsigned int *ct);




//
// main code
//

int main(int argc, const char **argv)
{


  /////////////////////
  // initialise card //
  /////////////////////
  findCudaDevice(argc, argv);


  // initialise CUDA timing
  bool use_timing = true;
  if( use_timing ){
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  }
  cudaEventCreate(&total_start);
  cudaEventCreate(&total_stop);
  float elapsed_parameters, elapsed_pmts, elapsed_detector, elapsed_vertices,
    elapsed_threads, elapsed_tof, elapsed_memory_tofs_dev, elapsed_tofs_copy_dev,
    elapsed_input, elapsed_memory_dev, elapsed_copy_dev, elapsed_kernel, 
    elapsed_threads_candidates, elapsed_candidates_memory_dev, elapsed_candidates_kernel,
    elapsed_candidates_copy_host, elapsed_coalesce, elapsed_gates, elapsed_free, elapsed_total,
    elapsed_tofs_free, elapsed_reset;
  bool use_verbose = true;


  ////////////////////
  // inspect device //
  ////////////////////
  print_gpu_properties();



  ///////////////////////
  // define parameters //
  ///////////////////////
  if( use_timing )
    start_c_clock();
  distance_between_vertices = 500.; // cm
  time_step_size  = 10; // ns
  threshold_number_of_pmts = 45;
  coalesce_time = 500.; // ns
  trigger_gate_up = 950.0; // ns
  trigger_gate_down = -400.0 -200; // ns
  output_txt = false;
  if( use_verbose ){
    printf(" --- user parameters \n");
    printf(" distance between test vertices = %f cm \n", distance_between_vertices);
    printf(" time step size = %d ns \n", time_step_size);
    printf(" threshold_number_of_pmts = %d \n", threshold_number_of_pmts);
    printf(" coalesce_time = %f ns \n", coalesce_time);
    printf(" trigger_gate_up = %f ns \n", trigger_gate_up);
    printf(" trigger_gate_down = %f ns \n", trigger_gate_down);
  }
  if( use_timing )
    elapsed_parameters = stop_c_clock();




  ////////////////
  // read PMTs  //
  ////////////////
  if( use_timing )
    start_c_clock();
  if( !read_the_pmts() ) return 0;
  if( use_timing )
    elapsed_pmts = stop_c_clock();


  /////////////////////
  // read detector ////
  /////////////////////
  if( use_timing )
    start_c_clock();
  if( !read_the_detector() ) return 0;
  if( use_timing )
    elapsed_detector = stop_c_clock();




  ////////////////////////
  // make test vertices //
  ////////////////////////
  if( use_timing )
    start_c_clock();
  make_test_vertices();
  if( use_timing )
    elapsed_vertices = stop_c_clock();



  //////////////////////////////
  // table of times_of_flight //
  //////////////////////////////
  if( use_timing )
    start_c_clock();
  make_table_of_tofs();
  if( use_timing )
    elapsed_tof = stop_c_clock();

  if( use_timing )
    start_cuda_clock();
  allocate_tofs_memory_on_device();
  if( use_timing )
    elapsed_memory_tofs_dev = stop_cuda_clock();


  if( use_timing )
    start_cuda_clock();
  fill_tofs_memory_on_device();
  if( use_timing )
    elapsed_tofs_copy_dev = stop_cuda_clock();


  ////////////////
  // read input //
  ////////////////
  if( use_timing )
    start_c_clock();
  if( !read_the_input() ) return 0;
  if( use_timing )
    elapsed_input = stop_c_clock();
  

  allocate_candidates_memory_on_host();


  ////////////////////////////////////////////////
  // set number of blocks and threads per block //
  ////////////////////////////////////////////////
  if( use_timing )
    start_c_clock();
  if( !setup_threads_for_tof() ) return 0;
  if( use_timing )
    elapsed_threads = stop_c_clock();


  start_total_cuda_clock();
  ///////////////////////////////
  // allocate memory on device //
  ///////////////////////////////
  if( use_timing )
    start_cuda_clock();
  allocate_correct_memory_on_device();
  if( use_timing )
    elapsed_memory_dev = stop_cuda_clock();


  //////////////////////////////////////
  // copy input into device variables //
  //////////////////////////////////////
  if( use_timing )
    start_cuda_clock();
  fill_correct_memory_on_device();
  if( use_timing )
    elapsed_copy_dev = stop_cuda_clock();



  ////////////////////
  // execute kernel //
  ////////////////////
  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- execute kernel \n");
  kernel_correct_times<<<number_of_kernel_blocks,number_of_threads_per_block>>>(device_n_pmts_per_time_bin);
  getLastCudaError("correct_kernel execution failed\n");
  if( use_timing )
    elapsed_kernel = stop_cuda_clock();



  /////////////////////////////////////
  // find candidates above threshold //
  /////////////////////////////////////
  if( use_timing )
    start_c_clock();
  if( !setup_threads_to_find_candidates() ) return 0;
  if( use_timing )
    elapsed_threads_candidates = stop_c_clock();


  if( use_timing )
    start_cuda_clock();
  allocate_candidates_memory_on_device();
  if( use_timing )
    elapsed_candidates_memory_dev = stop_cuda_clock();

  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- execute candidates kernel \n");
  kernel_find_vertex_with_max_npmts_in_timebin<<<number_of_kernel_blocks,number_of_threads_per_block>>>(device_n_pmts_per_time_bin, device_max_number_of_pmts_in_time_bin, device_vertex_with_max_n_pmts);
  getLastCudaError("candidates_kernel execution failed\n");
  if( use_timing )
    elapsed_candidates_kernel = stop_cuda_clock();

  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- copy candidates from device to host \n");
  checkCudaErrors(cudaMemcpy(host_max_number_of_pmts_in_time_bin,
			     device_max_number_of_pmts_in_time_bin,
			     n_time_bins*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(host_vertex_with_max_n_pmts,
			     device_vertex_with_max_n_pmts,
			     n_time_bins*sizeof(unsigned int),
			     cudaMemcpyDeviceToHost));
  if( use_timing )
    elapsed_candidates_copy_host = stop_cuda_clock();

  for(unsigned int time_bin = 0; time_bin<n_time_bins - 1; time_bin++){ // loop over time bins
    // n_time_bins - 1 as we are checking the i and i+1 at the same time
    
    if(host_max_number_of_pmts_in_time_bin[time_bin] > threshold_number_of_pmts) {

      if( use_verbose )
	printf(" time %f vertex (%f, %f, %f) npmts %d \n", (time_bin + 2)*time_step_size - time_offset, vertex_x[host_vertex_with_max_n_pmts[time_bin]], vertex_y[host_vertex_with_max_n_pmts[time_bin]], vertex_z[host_vertex_with_max_n_pmts[time_bin]], host_max_number_of_pmts_in_time_bin[time_bin]);

      candidate_trigger_pair_vertex_time.push_back(std::make_pair(host_vertex_with_max_n_pmts[time_bin],time_bin+2));
      candidate_trigger_npmts_in_time_bin.push_back(host_max_number_of_pmts_in_time_bin[time_bin]);
    }

  }

  if( use_verbose )
    printf(" n candidates: %d \n", candidate_trigger_pair_vertex_time.size());





  ///////////////////////
  // coalesce triggers //
  ///////////////////////
  if( use_timing )
    start_cuda_clock();
  coalesce_triggers();
  if( use_timing )
    elapsed_coalesce = stop_cuda_clock();




  //////////////////////////////////
  // separate triggers into gates //
  //////////////////////////////////
  if( use_timing )
    start_cuda_clock();
  separate_triggers_into_gates();
  if( use_timing )
    elapsed_gates = stop_cuda_clock();




  // deallocate all memory 
  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- deallocate memory \n");
  free_event_memories();
  if( use_timing )
    elapsed_free = stop_cuda_clock();

  elapsed_total = stop_total_cuda_clock();

  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- deallocate tofs memory \n");
  free_global_memories();
  if( use_timing )
    elapsed_tofs_free = stop_cuda_clock();


  // CUDA exit -- needed to flush the buffer which holds printf from each thread
  if( use_timing )
    start_cuda_clock();
  if( use_verbose )
    printf(" --- reset device \n");
  //  cudaDeviceReset();
  if( use_timing )
    elapsed_reset = stop_cuda_clock();


  if( use_timing ){
    printf(" user parameters time : %f ms \n", elapsed_parameters);
    printf(" read pmts execution time : %f ms \n", elapsed_pmts);
    printf(" read detector execution time : %f ms \n", elapsed_detector);
    printf(" make test vertices execution time : %f ms \n", elapsed_vertices);
    printf(" setup threads execution time : %f ms \n", elapsed_threads);
    printf(" setup threads candidates execution time : %f ms \n", elapsed_threads_candidates);
    printf(" make table of tofs execution time : %f ms \n", elapsed_tof);
    printf(" read input execution time : %f ms \n", elapsed_input);
    printf(" allocate tofs memory on device execution time : %f ms \n", elapsed_memory_tofs_dev);
    printf(" fill tofs memory on device execution time : %f ms \n", elapsed_tofs_copy_dev);
    printf(" deallocate tofs memory execution time : %f ms \n", elapsed_tofs_free);
    printf(" device reset execution time : %f ms \n", elapsed_reset);
    printf(" allocate memory on device execution time : %f ms (%f) \n", elapsed_memory_dev, elapsed_memory_dev/elapsed_total);
    printf(" fill memory on device execution time : %f ms (%f) \n", elapsed_copy_dev, elapsed_copy_dev/elapsed_total);
    printf(" correct kernel execution time : %f ms (%f) \n", elapsed_kernel, elapsed_kernel/elapsed_total);
    printf(" allocate candidates memory on device execution time : %f ms (%f) \n", elapsed_candidates_memory_dev, elapsed_candidates_memory_dev/elapsed_total);
    printf(" copy candidates to host execution time : %f ms (%f) \n", elapsed_candidates_copy_host, elapsed_candidates_copy_host/elapsed_total);
    printf(" candidates kernel execution time : %f ms (%f) \n", elapsed_candidates_kernel, elapsed_candidates_kernel/elapsed_total);
    printf(" coalesce triggers execution time : %f ms (%f) \n", elapsed_coalesce, elapsed_coalesce/elapsed_total);
    printf(" separate triggers into gates execution time : %f ms (%f) \n", elapsed_gates, elapsed_gates/elapsed_total);
    printf(" deallocate memory execution time : %f ms (%f) \n", elapsed_free, elapsed_free/elapsed_total);
  }
  printf(" total execution time : %f ms \n", elapsed_total);


  return 1;
}




//
// kernel routine
// 

// __global__ identifier says it's a kernel function
__global__ void kernel_correct_times(unsigned int *ct){


  // get unique id for each thread in each block == vertex index
  unsigned int vertex_index = threadIdx.x + blockDim.x*blockIdx.x;

  // skip if thread is assigned to nonexistent vertex
  if( vertex_index >= constant_n_test_vertices ) return;

  //    printf( " vertex_index %d threadidx %d blockdim %d blockid %d \n",
  //	    vertex_index, threadIdx.x, blockDim.x, blockIdx.x);


  unsigned int vertex_block = constant_n_time_bins*vertex_index;

  unsigned int vertex_block2 = constant_n_PMTs*vertex_index;
  for(unsigned int hit_index=0; hit_index<constant_n_hits; hit_index++){
    atomicAdd(&
	      ct[
		 vertex_block 
		 + int(floor(
			     (tex1Dfetch(tex_times,hit_index)
			      - tex1Dfetch(tex_times_of_flight,
					   vertex_block2 
					   + tex1Dfetch(tex_ids,hit_index) - 1
					   )
			      + constant_time_offset)/constant_time_step_size
			     )
		       )
		 ],1);
  }


  return;

}



