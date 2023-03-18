#define R_X 0
#define R_Y 1
#define V_X 2
#define V_Y 3

typedef struct paramsStruct {
    float b;
    float r_m_sq;
    float r_z_sq;
} paramsStruct;

void energy(
    int N,
    paramsStruct params,
    float *x);

int entry_point(
    int argc,
    const char *argv[],
    void (*compute)(int, int, int, float, paramsStruct, float *, float *));
