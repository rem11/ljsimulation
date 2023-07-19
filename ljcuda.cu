#include "common.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

float *vector_alloc(int N)
{
    return (float *)malloc(N * sizeof(float));
}

__global__ void ode(
    int N,
    paramsStruct params,
    float *f,
    float *x)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int M = N / 4;
    for (int i = index; i < M * M; i += stride)
    {
        int j = 4 * (i / M);
        int k = 4 * (i % M);
        if (j == k)
        {
            f[j + R_X] = x[j + V_X];
            f[j + R_Y] = x[j + V_Y];
        }
        else
        {
            float r_kj_x = x[k + R_X] - x[j + R_X];
            float r_kj_y = x[k + R_Y] - x[j + R_Y];
            float r_kj_sq = r_kj_x * r_kj_x + r_kj_y * r_kj_y;
            if (r_kj_sq < params.r_z_sq)
            {
                float r_rel_sq = params.r_m_sq / r_kj_sq;
                float coeff = (params.b / r_kj_sq) * (r_rel_sq * r_rel_sq * r_rel_sq -
                                                      r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq);
                atomicAdd(&f[j + V_X], r_kj_x * coeff);
                atomicAdd(&f[j + V_Y], r_kj_y * coeff);
            }
        }
    }
}

__global__ void vector_reset(int N, float *k_1, float *k_2, float *k_3, float *k_4)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride)
    {
        k_1[i] = 0;
        k_2[i] = 0;
        k_3[i] = 0;
        k_4[i] = 0;
    }
}

__global__ void vector_add(int N, float *x, float *k, float mul, float *r)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride)
    {
        r[i] = x[i] + k[i] * mul;
    }
}

__global__ void vector_rk(int N, float dt, float *x, float *k_1, float *k_2, float *k_3, float *k_4)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int j = index; j < N; j += stride)
    {
        x[j] += (dt / 6) * (k_1[j] + 2 * k_2[j] + 2 * k_3[j] + k_4[j]);
    }
}

void compute_iteration(int N, int M_stride, float dt, paramsStruct params, float *x, float *k_1, float *k_2, float *k_3, float *k_4, float *x_tmp)
{
    int blockSize = 64;
    int numBlocksOde = min(512, N * N / blockSize + 1);
    int numBlocksVector = min(512, N / blockSize + 1);    
    for (int i = 0; i < M_stride; i++)
    {
        vector_reset<<<numBlocksVector, blockSize>>>(N, k_1, k_2, k_3, k_4);
        ode<<<numBlocksOde, blockSize>>>(N, params, k_1, x);
        vector_add<<<numBlocksVector, blockSize>>>(N, x, k_1, dt / 2, x_tmp);
        ode<<<numBlocksOde, blockSize>>>(N, params, k_2, x_tmp);
        vector_add<<<numBlocksVector, blockSize>>>(N, x, k_2, dt / 2, x_tmp);
        ode<<<numBlocksOde, blockSize>>>(N, params, k_3, x_tmp);
        vector_add<<<numBlocksVector, blockSize>>>(N, x, k_3, dt, x_tmp);
        ode<<<numBlocksOde, blockSize>>>(N, params, k_4, x_tmp);
        vector_rk<<<numBlocksVector, blockSize>>>(N, dt, x, k_1, k_2, k_3, k_4);
    }
}

void compute(
    int N,
    int M_out,
    int M_stride,
    float dt,
    paramsStruct params,
    float *x_0,
    float *x_out)
{
    float *x, *k_1, *k_2, *k_3, *k_4, *x_tmp;
    cudaMalloc(&x, N * sizeof(float));
    cudaMalloc(&k_1, N * sizeof(float));
    cudaMalloc(&k_2, N * sizeof(float));
    cudaMalloc(&k_3, N * sizeof(float));
    cudaMalloc(&k_4, N * sizeof(float));
    cudaMalloc(&x_tmp, N * sizeof(float));
    cudaMemcpy(x, x_0, N * sizeof(float), cudaMemcpyDefault);
    memcpy(x_out, x_0, N * sizeof(float));
    for (int i = 1; i < M_out; i++)
    {
        compute_iteration(N, M_stride, dt, params, x, k_1, k_2, k_3, k_4, x_tmp);
        float *x_out_ptr = x_out + i * N;
        cudaMemcpy(x_out_ptr, x, N * sizeof(float), cudaMemcpyDefault);
        printf("Step: %d/%d, Progress: %g\n", i * M_stride, (M_out - 1) * M_stride, 100.0 * i / (M_out - 1));
        energy(N, params, x_out_ptr);
    }
}

int main(int argc, const char *argv[])
{
    return entry_point(argc, argv, &compute);
}
