#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <glm/glm.hpp>

enum Mode // tone reproduction mode
{
	NONE, 
	WARD,
	REINHARD,
    ADAPTIVE
};

std::vector<int> color_to_int( glm::vec3 in_color, float lbar, float lwmax, Mode mode )
{
	std::vector<int> out_color;
    double ldmax = 255.0f;

	if ( mode == NONE ) 
	{
		
		for ( uint32_t i = 0; i < 3; i++ )
		{
            int num = static_cast<int>(ldmax * in_color[i]);
            if ( num > 255 )
            {
                num = 255;
            }
            else if ( num < 0 )
            {
                num = 0;
            }
			out_color.push_back( num );
		}
	}
    else if ( mode == WARD )
    {
        ldmax = 2048.0f;
        float sf = std::pow( ((1.219f + std::pow(ldmax/2.0f, 0.4))/(1.219 + std::pow(lbar, 0.4))), 2.5 );// * (1.0f/ldmax);
        for ( uint32_t i = 0; i < 3; i++ )
		{
            out_color.push_back( static_cast<int>(sf * in_color[i]));
        }
    }
    else if ( mode == REINHARD )
    {
        ldmax = 1024.0f;
        float a = 0.18f;
        for ( int i = 0; i < 3; i++ )
        {
            float scaled_luminance = a/lbar * in_color[i];
            float reflectance = scaled_luminance / ( 1.0f + scaled_luminance );
            out_color.push_back( reflectance * ldmax );
        }
    }
    else if ( mode == ADAPTIVE )
    {
        ldmax = (1024.0f + 512.0f)/2;
        float lw = (0.27 * in_color.r) + (0.67 * in_color.g) + (0.06 * in_color.b);
        float ld = (1.0f/std::log10(lwmax + 1.0f)) * (std::log(lw + 1.0f) / std::log(2.0f+((std::pow(lw/lwmax, std::log(0.85)/std::log(0.5)))*8.0f))) * ldmax;
        std::cout << ld << std::endl;
        for ( int i = 0; i < 3; i++ )
        {
            out_color.push_back( in_color[i] * ld );
        }
    }

	return out_color;
}

void write_to_ppm( std::string filename, std::vector<std::vector<std::vector<double>>> pixels, Mode mode )
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
    float lbar = 0.0;
    float lwmax = 0.0f;
    for ( size_t j = 0; j < pixels.size(); j++ )
	{
		for ( size_t i = 0; i < pixels[0].size(); i++ )
		{
			r = pixels[j][i][0];
			g = pixels[j][i][1];
			b = pixels[j][i][2];
            float luminance = (0.27 * r) + (0.67 * g) + (0.06 * b);
			lbar += std::log(0.0001 + luminance);
            lwmax = std::max(lwmax, luminance);
		}
    }
    /*
    std::cout << lbar << std::endl;
    lbar /= pixels.size();
    std::cout << lbar << std::endl;
    lbar /= pixels[0].size();
    std::cout << lbar << std::endl;
    lbar = std::exp(lbar);
    std::cout << lbar << std::endl;
    */

	for ( size_t j = 0; j < pixels.size(); j++ )
	{
		for ( size_t i = 0; i < pixels[0].size(); i++ )
		{
			r = pixels[j][i][0];
			g = pixels[j][i][1];
			b = pixels[j][i][2];
			glm::vec3 color(r, g, b);
			auto pixel_color = color_to_int( color, lbar, lwmax, mode );
			file << ( pixel_color[0] ) << " ";
			file << ( pixel_color[1] ) << " ";
			file << ( pixel_color[2] ) << " ";
			file << std::endl;
		}
	}
	file.close();
}
