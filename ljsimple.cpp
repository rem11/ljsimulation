#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ode(
    int N,
    paramsStruct params,
    float *f,
    float *x)
{
#pragma omp parallel for
    for (int j = 0; j < N; j += 4)
    {
        float a_x = 0, a_y = 0;
        for (int k = 0; k < N; k += 4)
        {
            if (k != j)
            {
                float r_kj_x = x[k + R_X] - x[j + R_X];
                float r_kj_y = x[k + R_Y] - x[j + R_Y];
                float r_kj_sq = r_kj_x * r_kj_x + r_kj_y * r_kj_y;
                if (r_kj_sq < params.r_z_sq)
                {
                    float r_rel_sq = params.r_m_sq / r_kj_sq;
                    float coeff = (params.b / r_kj_sq) * (r_rel_sq * r_rel_sq * r_rel_sq -
                                                          r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq);
                    a_x += r_kj_x * coeff;
                    a_y += r_kj_y * coeff;
                }
            }
        }
        f[j + R_X] = x[j + V_X];
        f[j + R_Y] = x[j + V_Y];
        f[j + V_X] = a_x;
        f[j + V_Y] = a_y;
    }
}

float *vector_alloc(int N)
{
    return (float *)malloc(N * sizeof(float));
}

void vector_add(int N, float *x, float *k, float mul, float *r)
{
    for (int i = 0; i < N; i++)
    {
        r[i] = x[i] + k[i] * mul;
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
    float *x = vector_alloc(N);
    float *k_1 = vector_alloc(N);
    float *k_2 = vector_alloc(N);
    float *k_3 = vector_alloc(N);
    float *k_4 = vector_alloc(N);
    float *x_tmp = vector_alloc(N);
    memcpy(x, x_0, N * sizeof(float));
    int steps = (M_out - 1) * M_stride;
    for (int i = 0; i <= steps; i++)
    {
        if (i % M_stride == 0)
        {
            memcpy(x_out + (i / M_stride) * N, x, N * sizeof(float));
            printf("Step: %d/%d, Progress: %g\n", i, steps, 100.0 * i / steps);
            energy(N, params, x);
        }
        ode(N, params, k_1, x);
        vector_add(N, x, k_1, dt / 2, x_tmp);
        ode(N, params, k_2, x_tmp);
        vector_add(N, x, k_2, dt / 2, x_tmp);
        ode(N, params, k_3, x_tmp);
        vector_add(N, x, k_3, dt, x_tmp);
        ode(N, params, k_4, x_tmp);
        for (int j = 0; j < N; j++)
        {
            x[j] += (dt / 6) * (k_1[j] + 2 * k_2[j] + 2 * k_3[j] + k_4[j]);
        }
    }
    free(x);
    free(k_1);
    free(k_2);
    free(k_3);
    free(k_4);
    free(x_tmp);
}

int main(int argc, const char *argv[])
{
    return entry_point(argc, argv, &compute);
}
