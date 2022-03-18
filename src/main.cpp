// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <cstdlib>
#include <ctime>

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
			std::string command = argv[i] ;
			std::cout << "command \"" << command << "\": does not exist" << std::endl;
		}
	}

	// init pixel array
	double aspect_ratio = (double)width / (double)height;
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


	PhongMaterial glass_material( glm::vec3(0.0, 0.0, 0.1), glm::vec3(0.0, 0.0, 1.0), glm::vec3(1.0, 1.0, 1.0), &scene );
	Sphere glass_sphere( glm::vec3(0.0, 0.2, 0.0), 0.5, &glass_material, "glass" );

	PhongMaterial solid_material( glm::vec3(0.1, 0.1, 0.0), glm::vec3(1.0, 1.0, 0.0), glm::vec3(1.0, 1.0, 1.0), &scene );
	Sphere solid_sphere( glm::vec3(0.7, -0.2, -1.0), 0.5, &solid_material, "solid" );

	CheckeredPhongMaterial floor_material( glm::vec3(0.1, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0), glm::vec3(1.0, 1.0, 0.0), 5, glm::vec3(1.0, 1.0, 1.0), &scene );
	CheckeredPhongMaterial floor_material_1( glm::vec3(0.1, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0), glm::vec3(1.0, 1.0, 0.0), 5, glm::vec3(1.0, 1.0, 1.0), &scene );
	CheckeredPhongMaterial floor_material_2( glm::vec3(0.1, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0), glm::vec3(1.0, 1.0, 0.0), 5, glm::vec3(1.0, 1.0, 1.0), &scene );
	Floor floor( glm::vec3(0.75 - 1.5, -0.7, -0.5 - 1.5), // top left : ACTUALLY TOP LEFT
				 glm::vec3(0.75 + 1.5, -0.7, -0.5 - 1.5), // top right : ACTUALLY TOP RIGHT 
				 glm::vec3(0.75 - 1.5, -0.7, -0.5 + 1.5), // bottom left : ACTUALLY BOTTOM LEFT
				 glm::vec3(0.75 + 1.5, -0.7, -0.5 + 1.5), // bottom right : ACTUALLY BOTTOM RIGHT
				 &floor_material); // color
				 

	glm::vec3 cam_pos = glm::vec3(0.0, 0.2, 2.0);
	glm::vec3 look_at = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);
	Camera camera( cam_pos, look_at, up, aspect_ratio ); 

	Light light1( glm::vec3(1.0, 1.0, 1.0), glm::vec3(0.0, 20.0, 20.0) );
	Light light2( glm::vec3(1.0, 1.0, 1.0), glm::vec3(0.0, 2.0, 0.0) );

	scene.add_object( &glass_sphere );
	scene.add_object( &solid_sphere );
	scene.add_object( &floor );
	scene.add_light( &light1 );
	//scene.add_light( &light2 );
/*
	SoftPhongMaterial ground( glm::vec3(0.1, 0.1, 0.1), glm::vec3(0.0, 1.0, 0.0), glm::vec3(1.0, 1.0, 1.0), &scene );
	SoftPhongMaterial trimat( glm::vec3(0.1, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0), glm::vec3(1.0, 1.0, 1.0), &scene );

	Triangle shadow_caster( glm::vec3(0.0, -0.3, 0.0 - 0.5), 
							glm::vec3(-0.5,-0.3, 0.5 - 0.5),
							glm::vec3(0.5, -0.3, 0.5 - 0.5),
							&trimat );
	Floor floor( glm::vec3(0.0 - 1.5, -0.7, -0.5 + 1.5), // top left
				 glm::vec3(0.0 + 1.5, -0.7, -0.5 + 1.5), // top right
				 glm::vec3(0.0 - 1.5, -0.7, -0.5 - 1.5), // bottom left
				 glm::vec3(0.0 + 1.5, -0.7, -0.5 - 1.5), // bottom right
				 &ground); // color

	SphereLight lightsource( glm::vec3(0.0, 0.5, 0.0), 0.1, glm::vec3(1.0, 1.0, 1.0) );
	
	glm::vec3 cam_pos = glm::vec3(0.0, 0.2, 2.0);
	glm::vec3 look_at = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);
	Camera camera( cam_pos, look_at, up, aspect_ratio ); 

	scene.add_geometry_light( &lightsource );
	scene.add_object( &floor );
	//scene.add_object( &fake_light );
	scene.add_object( &shadow_caster );
*/

	/*
	 * Begin Ray casting main!
	 */
	// variable setup
	Intersection closest_intersection;
	double closest_intersect_dist = 0;
	Object *closest_object = nullptr;
	uint16_t num_samples = 10;

	std::srand(std::time(nullptr)); 

	// ray tracing!
	for ( auto y = 0; y < height; y++ )
	{
		for ( auto x = 0; x < width; x++ )
		{
			// generate ray for the pixel we are on

			for ( auto sample = 0; sample < num_samples; sample++ )
			{
				auto sample_u = (x + std::rand()/(RAND_MAX + 1.0))/(width + 1.0);
				auto sample_v = (y + std::rand()/(RAND_MAX + 1.0))/(height + 1.0);
				Ray test_ray = camera.spawn_ray(sample_u, sample_v);
				
				// intersection tests
				closest_intersection = Intersection();
				closest_object = nullptr;

				auto intersections = scene.intersect_geometry( test_ray );

				for ( uint32_t i = 0; i < intersections.size(); i++ )
				{
					auto temp_container = intersections[i];
					auto object = std::get<0>(temp_container);
					auto intersection = std::get<1>(temp_container);
					auto intersect_point = intersection.position;
					//auto intersect_normal = intersection.normal; 
					if ( !closest_intersection.exists() || glm::distance(camera.position, intersect_point) < closest_intersect_dist  )
					{
						closest_intersection = intersection;
						closest_intersect_dist = glm::distance(camera.position, intersect_point);
						closest_object = object;
					}
				}

				// draw to the pixel!
				glm::vec3 color;
				if ( closest_intersection.exists() )
				{
					color = closest_object->get_color(test_ray, closest_intersection, closest_object);
				}
				else
				{
					color = scene.get_sky_color(x, y, width, height);					
				}
				pixels[y][x][0] += color.x / num_samples;
				pixels[y][x][1] += color.y / num_samples;
				pixels[y][x][2] += color.z / num_samples;
			}
		}
	}

	// write to image
	write_to_ppm( filename, pixels );
	return 0;
}