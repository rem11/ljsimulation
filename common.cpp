#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void energy(
    int N,
    paramsStruct params,
    float *x)
{
    float E_k = 0;
    float E_p = 0;
    float r_rel_z_sq = params.r_m_sq / params.r_z_sq;
    float E_p_z = (params.b / 12) * (r_rel_z_sq * r_rel_z_sq * r_rel_z_sq * r_rel_z_sq * r_rel_z_sq * r_rel_z_sq -
                                     2 * r_rel_z_sq * r_rel_z_sq * r_rel_z_sq);
    for (int j = 0; j < N; j += 4)
    {
        E_k += (x[j + V_X] * x[j + V_X] + x[j + V_Y] * x[j + V_Y]) / 2;
        for (int k = 0; k < j; k += 4)
        {
            float r_kj_x = x[k + R_X] - x[j + R_X];
            float r_kj_y = x[k + R_Y] - x[j + R_Y];
            float r_kj_sq = r_kj_x * r_kj_x + r_kj_y * r_kj_y;
            if (r_kj_sq < params.r_z_sq)
            {
                float r_rel_sq = params.r_m_sq / r_kj_sq;
                E_p += ((params.b / 12) * (r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq * r_rel_sq -
                                           2 * r_rel_sq * r_rel_sq * r_rel_sq) -
                        E_p_z);
            }
        }
    }
    printf("E_k=%g E_p=%g E=%g\n", E_k, E_p, E_k + E_p);
}

int entry_point(
    int argc,
    const char *argv[],
    void (*compute)(int, int, int, float, paramsStruct, float *, float *))
{
    if (argc < 3)
    {
        printf("Insufficient arguments, please specify input and output files");
        return 1;
    }
    time_t time1;
    time(&time1);
    printf("Reading input file...\n");
    FILE *f_input = fopen(argv[1], "r");
    if (f_input == NULL)
    {
        printf("Can't open input file");
        return 1;
    }
    int N;
    float T, dt_out, dt, b, r_m, r_z;
    if (fscanf(f_input, "%d %g %g %g %g %g %g", &N, &T, &dt_out, &dt, &b, &r_m, &r_z) < 6)
    {
        printf("Invalid input file format");
        return 1;
    }
    float *x_0 = (float *)malloc(N * 4 * sizeof(float));
    for (int i = 0; i < N * 4; i += 4)
    {
        if (fscanf(f_input, "%g %g %g %g", &x_0[i + R_X], &x_0[i + R_Y], &x_0[i + V_X], &x_0[i + V_Y]) < 4)
        {
            printf("Invalid input file format");
            return 1;
        }
    }
    fclose(f_input);
    time_t time2;
    time(&time2);
    printf("Done in %ld s\n", time2 - time1);
    time(&time1);
    printf("Computing...\n");
    int M_out = trunc(T / dt_out);
    int M_stride = trunc(dt_out / dt);
    float *x = (float *)malloc(M_out * N * 4 * sizeof(float));
    paramsStruct params;
    params.b = b;
    params.r_m_sq = r_m * r_m;
    params.r_z_sq = r_z * r_z;
    compute(N * 4, M_out, M_stride, dt, params, x_0, x);
    time(&time2);
    printf("Done in %ld s\n", time2 - time1);
    time(&time1);
    printf("Writing output file...\n");
    FILE *f_output = fopen(argv[2], "w");
    if (f_output == NULL)
    {
        printf("Can't open output file");
        return 1;
    }
    fprintf(f_output, "%d %d\n", N, M_out);
    for (int i = 0; i < M_out; i++)
    {
        fprintf(f_output, "%g\n", dt * M_stride * i);
        int i_o = i * N * 4;
        for (int j = 0; j < N * 4; j += 4)
        {
            fprintf(f_output, "%g %g %g %g\n", x[i_o + j + R_X], x[i_o + j + R_Y], x[i_o + j + V_X], x[i_o + j + V_Y]);
        }
    }
    fclose(f_output);
    time(&time2);    
    printf("Done in %ld s\n", time2 - time1);
    free(x_0);
    free(x);
    return 0;
}