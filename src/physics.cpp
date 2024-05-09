#include "./physics.h"



void Physics::SetBnd(int b, float x[][SIZE], int N) {
	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 1; i < N - 1; i++) {
				x[i][0] = b == 2 ? -x[i][1]: x[i][1];
				x[i][N-1] = b == 2 ? -x[i][N-2]: x[i][N-2];
		}

		#pragma omp for
		for(int j = 1; j < N - 1; j++) {
				x[0][j] = b == 1 ? -x[1][j] : x[1][j];
				x[N-1][j] = b == 1 ? -x[N-2][j] : x[N-2][j];
		}
	}

	x[0][0] = 0.33f * (x[1][0]
								+ x[0][1]
								+ x[0][0]);
	x[0][N-1] = 0.33f * (x[1][N-1]
								+ x[0][N-2]
								+ x[0][N-1]);
	x[N-1][0] = 0.33f * (x[N-2][0]
								+ x[N-1][1]
								+ x[N-1][0]);
	x[N-1][N-1] = 0.33f * (x[N-2][N-1]
								+ x[N-1][N-2]
								+ x[N-1][N-1]);
}

void Physics::LinSolve(int b, float x[][SIZE], float x0[][SIZE], float a, float c, int iter, int N) {
	float cRecip = 1.0 / c;
	for (int k = 0; k < iter; k++) {
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
					x[i][j] = (x0[i][j]+ a*(x[i+1][j]+x[i-1][j]+x[i][j+1]+x[i][j-1]+x[i][j]+x[i][j]
														)) * cRecip;
			}
		}
	}
	this->SetBnd(b, x, N);
}

void Physics::Diffuse(int b, float x[][SIZE], float x0[][SIZE], float diff, float dt, int iter, int N) {
	float a = dt * diff * (N - 2) * (N - 2);
	this->LinSolve(b, x, x0, a, 1 + 6 * a, iter, N);
}

void Physics::Project(float vx[][SIZE], float vy[][SIZE], float p[][SIZE], float div[][SIZE], int iter, int N) {
	#pragma omp parallel
	{
		#pragma omp for
		for (int j = 1; j < N - 1; j++) {
			#pragma omp for
			for (int i = 1; i < N - 1; i++) {
				div[i] [j]= -0.5f*(
												vx[i+1][j]
												-vx[i-1][j]
												+vy[i][j+1]
												-vy[i][j-1]
											)/N;
				p[i][j] = 0;
			}
		}
	}

	this->SetBnd(0, div, N);
	this->SetBnd(0, p, N);
	this->LinSolve(0, p, div, 1, 6, iter, N);

	#pragma omp parallel
	{
		#pragma omp for
		for (int j = 1; j < N - 1; j++) {
			#pragma omp for
			for (int i = 1; i < N - 1; i++) {
				vx[i][j]-= 0.5f * (p[i+1][j] - p[i-1][j]) * N;
				vy[i][j] -= 0.5f * (p[i][j+1] -p[i][j-1]) * N;
			}
		}
	}
	this->SetBnd(1, vx, N);
	this->SetBnd(2, vy, N);
}

void Physics::Advect(int b, float d[][SIZE], float d0[][SIZE], float vx[][SIZE], float vy[][SIZE], float dt, int N) {
	float i0, i1, j0, j1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);

	float s0, s1, t0, t1;
	float tmp1, tmp2, x, y;

	float Nfloat = N;
	float ifloat, jfloat;

	int i, j;


	for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * vx[i][j];
			tmp2 = dty * vy[i][j];
			x = ifloat - tmp1;
			y = jfloat - tmp2;

			if(x < 0.5f) x = 0.5f;
			if(x > Nfloat + 0.5f) x = Nfloat + 0.5f;
			i0 = ::floorf(x);
			i1 = i0 + 1.0f;
			if(y < 0.5f) y = 0.5f;
			if(y > Nfloat + 0.5f) y = Nfloat + 0.5f;
			j0 = ::floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0;
			s0 = 1.0f - s1;
			t1 = y - j0;
			t0 = 1.0f - t1;

			int i0i = i0;
			int i1i = i1;
			int j0i = j0;
			int j1i = j1;

			d[i][j] =
				s0 * (t0 * d0[i0i][j0i] + t1 * d0[i0i][j1i]) +
				s1 * (t0 * d0[i1i][j0i] + t1 * d0[i1i][j1i]);
        }
    }
	this->SetBnd(b, d, N);
}