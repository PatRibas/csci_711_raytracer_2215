CC=g++
CFLAGS=-std=c++17 -Wall -Werror 
FILES=main.cpp ppm.cpp ppm.hpp geometry.hpp  
FILENAME=-o raytracer

raytracer: $(FILES) 
	$(CC) $(CFLAGS) $(FILENAME) $(FILES) 

clean:
	rm -f $(FILENAME)
