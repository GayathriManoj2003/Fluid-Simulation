#include "./container.h"
#include <iostream>
#include <cmath>

int get1DIndex(int x, int y, int N);

Container::Container() : physics(Physics()) {}

Container::Container(float dt, float diff, float visc) : physics(Physics()){
	this->size = SIZE;
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;

	this->InitArr(this->px, SIZE*SIZE);
	this->InitArr(this->py, SIZE*SIZE);
	this->InitArr(this->x, SIZE*SIZE);
	this->InitArr(this->y, SIZE*SIZE);
	this->InitArr(this->previousDensity, SIZE*SIZE);
	this->InitArr(this->density, SIZE*SIZE);
}

void Container::InitArr(float arr[], int size) {
	#pragma omp parallel
		#pragma omp for
		for (int i = 0; i < size; i++) {
			arr[i] = 0;
		}1
}

void Container::AddDensity(float x, float y, float amount) {
	this->density[get1DIndex(x,y,this->size)] += amount;
}

void Container::AddVelocity(float x, float y, float px, float py) {
	int index = get1DIndex(x,y,this->size);

	this->x[index] += px;
	this->y[index] += py;
}

void Container::Step() {
	// Diffuse Velocities
	this->physics.Diffuse(1, this->px, this->x, this->visc, this->dt, 16, this->size);
	this->physics.Diffuse(2, this->py, this->y, this->visc, this->dt, 16, this->size);

	// Fix Up Velocities
	this->physics.Project(this->px, this->py, this->x, this->y, 16, this->size);

	// Move Velocities
	this->physics.Advect(1, this->x, this->px, this->px, this->py, this->dt, this->size);
	this->physics.Advect(2, this->y, this->py, this->px, this->py, this->dt, this->size);

	// Fix Up Velocities
	this->physics.Project(this->x, this->y, this->px, this->py, 16, this->size);

	// Diffuse and Move Densities
	this->physics.Diffuse(0, this->previousDensity, this->density, this->diff, this->dt, 16, this->size);
	this->physics.Advect(0, this->density, this->previousDensity, this->x, this->y, this->dt, this->size);
}

float Container::MapToRange(float val, float minIn, float maxIn, float minOut, float maxOut) {
	float x = (val - minIn) / (maxIn - minIn);
	float result = minOut + (maxOut - minOut) * x;
	return (result < minOut) ? minOut : (result > maxOut) ? maxOut : result;
}

void Container::Render(sf::RenderWindow& win, int color) {
	win.clear();
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < this->size; i++) {
			#pragma omp for
			for(int j = 0; j < this->size; j++) {
				sf::RectangleShape rect;
				rect.setSize(sf::Vector2f(SCALE, SCALE));
				rect.setPosition(j * SCALE, i * SCALE);

				switch (color) {
					case 0:
						rect.setFillColor(sf::Color(255, 255, 255, (this->density[get1DIndex(i,j,this->size)] > 255) ? 255 : this->density[get1DIndex(i,j,this->size)]));
						break;
					case 1: {
							int r = (int)this->MapToRange(this->x[get1DIndex(i,j,this->size)], -0.01f, 0.01f, 0, 255);
							int g = (int)this->MapToRange(this->y[get1DIndex(i,j,this->size)], -0.01f, 0.01f, 0, 255);
							rect.setFillColor(sf::Color(r, g, 200));
							break;
						}
				};

				win.draw(rect);
			}
		}
	}
}

void Container::FadeDensity(int size) {
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < size; i++) {
			float d = this->density[i];
			density[i] = (d - 0.05f < 0) ? 0 : d - 0.05f;
		}
	}
}
