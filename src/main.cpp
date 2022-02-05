// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// personal includes
#include "ppm.h"
#include <glm/vec3.hpp>


int main(int argc, char **argv)
{
	// process cmd input
	if ( argc < 3 )
	{
		std::cout << "error: expected " << std::endl 
				  << "\t./prog <width> <height> <options>" << std::endl 
				  << "but did not receive image dimensions!" << std::endl;
		std::exit(1);
	}
	int width, height;
	for ( auto i = 1; i < argc; i++ )
	{
		if ( i == 1 )
		{
			width = std::stoi(argv[i]);
		}
		else if ( i == 2 )
		{
			height = std::stoi(argv[i]);
		}
		else
		{
			std::string command =  argv[i] ;
			std::cout << command << std::endl;
		}
	}

	// init pixel array
	std::vector<std::vector<std::vector<double>>> pixels;
	std::vector<double> color{ 0.0f, 0.0f, 0.0f };
	std::vector<std::vector<double>> column;
	for ( auto i = 0; i < width; i++ )
	{
		column.push_back( color );
	}

	for ( auto i = 0; i < height; i++ )
	{
		pixels.push_back( column );
	}

	// write background color to pixels
	// TODO: ray tracing!
	for ( uint64_t y = 0; y < height; y++ )
	{
		for ( uint64_t x = 0; x < width; x++ )
		{
			pixels[y][x][0] = (135.0 + ((255.0 - 135.0) * ((double)y / (double)height))) / 255.0;
			pixels[y][x][1] = (206.0 + ((255.0 - 206.0) * ((double)y / (double)height))) / 255.0;
			pixels[y][x][2] = (235.0 + ((255.0 - 235.0) * ((double)y / (double)height))) / 255.0;
		}
	}

	// write to image
	write_to_ppm( "img.ppm", pixels );
	return 0;
}