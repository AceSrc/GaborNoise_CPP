main: main.cpp gabor.hpp
	g++ main.cpp fft.cpp -o main `pkg-config opencv --cflags --libs` -g -lpthread

run: main
	./main

