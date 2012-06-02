using namespace std;

typedef double Data_t;

__global__ void cuda_distance_kernel(double *a, double *b, double *out) {
	int i = threadIdx.x;
	out[i] = abs(a[i] + b[i]);
}

void cuda_distance(const Data_t* a, const Data_t* b, size_t dim, Data_t* out) {
	// device copies
	Data_t *d_a, *d_b, *d_out;

	double size = dim * sizeof(double);
	cudaMalloc ((void **) &d_a, size);
	cudaMalloc ((void **) &d_b, size);
	cudaMalloc ((void **) &d_out, size);

	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
	
	cuda_distance_kernel<<<(dim+32-1)/32,32>>>(d_a, d_b, d_out);

	cudaMemcpy(out, d_out, size, cudaMemcpyDeviceToHost);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
	cudaThreadExit();
		
	return;
}