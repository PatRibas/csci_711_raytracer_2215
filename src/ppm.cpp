#include "ppm.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int color_to_int( double color )
{
	double max = 255.0f;
	return static_cast<int>(max * color);
}

void write_to_ppm( std::string filename, std::vector<std::vector<std::vector<double>>> pixels )
{
	if ( filename.substr(filename.length() - 4, filename.length()) != ".ppm" )
	{
		filename.append(".ppm");
	}

	double r,g,b;
	std::ofstream file; 
	file.open(filename);

	file << "P3";
	file << std::endl;
	file << pixels.at(0).size();
	file << " ";
	file << pixels.size();
	file << std::endl;
	file << "255";
	file << std::endl;
	for ( size_t j = 0; j < pixels.size(); j++ )
	{
		for ( size_t i = 0; i < pixels[0].size(); i++ )
		{
			r = pixels[j][i][0];
			g = pixels[j][i][1];
			b = pixels[j][i][2];
			file << color_to_int( r ) << " ";
			file << color_to_int( g ) << " ";
			file << color_to_int( b ) << " ";
			file << std::endl;
		}
	}
	file.close();
}
