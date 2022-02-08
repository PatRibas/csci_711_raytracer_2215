// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <optional>
#include <utility>

// GLM includes
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp> 

// personal includes
#include "ppm.hpp"
#include "geometry.hpp"



int main(int argc, char **argv)
{
	// process cmd input
	if ( argc < 3 )
	{
		std::cout << "error: expected " << std::endl 
				  << "\t./prog <width> <height> [filename] [options]" << std::endl 
				  << "but did not receive image dimensions!" << std::endl;
		std::exit(1);
	}
	int width, height;
	std::string filename = "img.ppm";
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
		else if ( i == 3 )
		{
			filename = argv[i];
		}
		else
		{
			std::string command =  argv[i] ;
			std::cout << "command \"" << command << "\": does not exist" << std::endl;
		}
	}

	// init pixel array
	std::vector<std::vector<std::vector<double>>> pixels;
	std::vector<std::vector<double>> column;
	for ( auto i = 0; i < width; i++ )
	{
		column.push_back( std::vector<double>{ 0.0f, 0.0f, 0.0f } );
	}

	for ( auto i = 0; i < height; i++ )
	{
		pixels.push_back( column );
	}

	// scene setup
	Scene scene( glm::vec3(135.0, 206.0, 235.0) );
	Sphere glass_sphere( glm::vec3(0.0, 0.2, 0.0), 0.5, glm::vec3(0.0, 0.0, 1.0) );
	Sphere solid_sphere( glm::vec3(0.7, -0.2, -1.0), 0.5, glm::vec3(1.0, 1.0, 0.0) );
	Floor floor( glm::vec3(0.75 - 1.5, -0.7, -0.5 + 1.5), 
				 glm::vec3(0.75 + 1.5, -0.7, -0.5 + 1.5), 
				 glm::vec3(0.75 - 1.5, -0.7, -0.5 - 1.5), 
				 glm::vec3(0.75 + 1.5, -0.7, -0.5 - 1.5), glm::vec3(1.0, 0.0, 0.0));

	glm::vec3 cam_pos = glm::vec3(0.0, 0.2, 2.0);
	glm::vec3 look_at = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);
	Camera camera( cam_pos, look_at, up ); 
	scene.add_object( &glass_sphere );
	scene.add_object( &solid_sphere );
	scene.add_object( &floor );

	// variable setup
	double aspect_ratio = width / height;
	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio; 
	double focal_length = 1.0;
	glm::vec3 horizontal(viewport_width, 0, 0);
	glm::vec3 vertical(0, viewport_width, 0);
	glm::vec3 upper_left_corner = camera.position - horizontal*0.5f + vertical*0.5f - glm::vec3(0, 0, focal_length);
	std::optional<std::pair<glm::vec3, glm::vec3>> intersection, closest_intersection;
	glm::vec3 intersect_point, intersect_normal;
	std::pair<glm::vec3, glm::vec3> container;
	double closest_intersect_dist = 0;
	Object *closest_object = nullptr;

	// ray tracing!
	for ( auto y = 0; y < height; y++ )
	{
		for ( auto x = 0; x < width; x++ )
		{
			// generate ray for the pixel we are on
			float u = double(x) / (width-1);
            float v = double(y) / (height-1);
            Ray test_ray(camera.position, upper_left_corner + horizontal*u - vertical*v - camera.position);
			
			// intersection tests
			intersection = std::nullopt;
			closest_intersection = std::nullopt;
			closest_object = nullptr;
			for ( Object *obj : scene.get_objects() )
			{
				intersection = obj->intersect( test_ray );
				if ( intersection.has_value() )
				{
					container = intersection.value();
					intersect_point = std::get<0>(container); // TODO: record this for just closest object
					intersect_normal = std::get<1>(container);
					if ( closest_object == nullptr || glm::distance(camera.position, intersect_point) < closest_intersect_dist  )
					{
						closest_intersection = intersection;
						closest_intersect_dist = glm::distance(camera.position, intersect_point);
						closest_object = obj;
					}
				}
			}

			// draw to the pixel!
			if ( closest_intersection.has_value() )
			{
				pixels[y][x][0] = closest_object->material.x;
				pixels[y][x][1] = closest_object->material.y;
				pixels[y][x][2] = closest_object->material.z;
			}
			else
			{
				pixels[y][x][0] = (135.0 + ((255.0 - 135.0) * ((double)y / (double)height))) / 255.0;
				pixels[y][x][1] = (206.0 + ((255.0 - 206.0) * ((double)y / (double)height))) / 255.0;
				pixels[y][x][2] = (235.0 + ((255.0 - 235.0) * ((double)y / (double)height))) / 255.0;
			}
		}
	}

	// write to image
	write_to_ppm( filename, pixels );
	return 0;
}