make: four.cpp
	g++ four.cpp -o four -lm -fopenmp
simd: simd.cpp
	g++ simd.cpp -o five -lm -fopenmp
