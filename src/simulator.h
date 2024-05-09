#include<vector>

#include "./container.h"

class Simulator {
private:
	sf::RenderWindow window;
	Container container;
	int color = 0;
public:
	Simulator() : container(Container(0.2f, 0, 0.0000001f)) {
		window.create(sf::VideoMode(SIZE*SCALE, SIZE*SCALE),
		 		  "Fluid Simulation",
		 		   sf::Style::Titlebar | sf::Style::Close);
		// this->win.create(sf::VideoMode(SIZE * SCALE, SIZE * SCALE), "Fluid Simulation", sf::Style::Titlebar | sf::Style::Close, sf::ContextSettings(0, 0, 10));
	};
	void Run();
};
