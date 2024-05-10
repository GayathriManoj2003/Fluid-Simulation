build: ./src/*.cpp
	g++ -std=c++17 -o main ./src/*.cpp -lsfml-graphics -lsfml-window -lsfml-system
	./main