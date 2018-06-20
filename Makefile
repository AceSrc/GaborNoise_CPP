tst: tst.cpp fft.cpp gabor.hpp
	g++ tst.cpp fft.cpp -o tst -std=c++11 -pthread -Wall -g

single: tst
	time ./tst 0


main: main.cpp gabor.hpp
	g++ main.cpp fft.cpp -o main `pkg-config opencv --cflags --libs` -g -lpthread

run: main
	./main

run1: tst
	./tst 0
	./tst 1
	./tst 2
	./tst 3
	convert noise.png fourier.png -append 1.png
	convert gabor.png gabor_fourier.png -append 2.png
	convert 1.png 2.png +append result.png
