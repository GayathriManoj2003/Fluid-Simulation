#include<vector>
#include <SFML/Graphics.hpp>

#include <omp.h>
#include "./physics.h"


class Container {
    private:
        Physics physics;

        int size;

        float dt;
        float diff;
        float visc;

        float px[SIZE][SIZE];
        float py[SIZE][SIZE];

        float x[SIZE][SIZE];
        float y[SIZE][SIZE];

        float previousDensity[SIZE][SIZE];
        float density[SIZE][SIZE];

        void InitArr(float arr[][SIZE], int size);
        float MapToRange(float value, float minIn, float maxIn, float minOut, float maxOut);
    public:
        Container();
        Container(float dt, float diff, float visc);

        void AddDensity(float x, float y, float amount);
        void AddVelocity(float x, float y, float px, float py);
        void Step();
        void Render(sf::RenderWindow& win, int color);
        void FadeDensity(int size);
};