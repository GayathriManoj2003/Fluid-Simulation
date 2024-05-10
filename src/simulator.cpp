#include <iostream>
#include "./simulator.h"

void Simulator::Run() {
	sf::Vector2i previousMouse = sf::Mouse::getPosition(window);
	sf::Vector2i currentMouse = sf::Mouse::getPosition(window);

	while (window.isOpen()) {
		sf::Event e;
		while (window.pollEvent(e)) {
			switch (e.type) {
				case sf::Event::Closed:
					window.close();
					break;
			}
		}
		int density = 1000;

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
			this->container.AddDensity(currentMouse.y/SCALE, currentMouse.x/SCALE, density);

		currentMouse = sf::Mouse::getPosition(window);

		float amountX = currentMouse.x - previousMouse.x;
		float amountY = currentMouse.y - previousMouse.y;
		float velX = amountX/10;
		float velY = amountY/10;

		this->container.AddVelocity(currentMouse.y/SCALE, currentMouse.x/SCALE, velX, velY);

		previousMouse = currentMouse;

		this->container.Step();
		this->container.Render(window, color);
		this->container.FadeDensity(SIZE);

		window.display();
	}
}

