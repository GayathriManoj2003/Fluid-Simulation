#include <omp.h>
#include<cmath>
#include "./const.h"

class Physics {
    private:
        void LinSolve(int b, float x[][SIZE], float x0[][SIZE], float a, float c, int iter, int N);
        void SetBnd(int b, float x[][SIZE], int N);
    public:
        void Diffuse(int b, float x[][SIZE], float x0[][SIZE], float diff, float dt, int iter, int N);
        void Project(float vx[][SIZE], float vy[][SIZE], float p[][SIZE], float div[][SIZE], int iter, int N);
        void Advect(int b, float d[][SIZE], float d0[][SIZE], float vx[][SIZE], float vy[][SIZE], float dt, int N);
};