#ifndef PPM
#define PPM

#include <string>
#include <vector>

int color_to_int( double color );
void write_to_ppm( std::string filename, std::vector<std::vector<std::vector<double>>> pixels );

#endif