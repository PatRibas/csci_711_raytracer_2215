#ifndef GEOMETRY
#define GEOMETRY

#include <optional>
#include <vector>
#include <utility>
#include <iostream>
#include <random>
#include <chrono> 
#include <cmath>

#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtx/normal.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/gtx/vec_swizzle.hpp> 
#include <glm/gtx/vector_angle.hpp>

class Ray;
class Intersection;
class Light;
class Object;
class Camera;
class Scene;
class Sphere;
class Triangle;
class Floor;

class Material;
class BasicMaterial;
class PhongMaterial;


void vec3_print( glm::vec3 vector )
{
    std::cout << "[" << vector.x << ", " << vector.y << ", " << vector.z << "]" << std::endl;
}

class Ray
{
public:
    glm::vec3 origin;
    glm::vec3 direction;

    Ray(glm::vec3 origin, glm::vec3 direction)
    {
        this->origin = origin;
        this->direction = glm::normalize(direction);
    }

    Ray()
    {
        this->origin = glm::vec3();
        this->direction = glm::vec3();
    }

    ~Ray(){}

};

class Intersection
{
private:
    bool has_value;
public:
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;
    Object *object_intersected;

    bool exists()
    {
        return has_value;
    }

    Intersection()
    {
        has_value = false;
    }

    Intersection(glm::vec3 position, glm::vec3 normal, glm::vec2 uv, Object *object_intersected)
    {
        has_value = true;
        this->position = position;
        this->normal = normal;
        this->uv = uv;
        this->object_intersected = object_intersected;
    }
};

class Light
{
public:
    glm::vec3 color;
    glm::vec3 position;
    float intensity;
    
    Light( glm::vec3 color, glm::vec3 position, float intensity )
    {
        this->color = color;
        this->position = position;
        this->intensity = intensity;
    }

    float intensity_to( glm::vec3 pos )
    {
        float distance = glm::distance(position, pos);
        return intensity * 1.0f/( distance * distance );
    }
};


class Object
{
public:
    Material *material;
    std::string name;

    // Intersect should return vectors of the form:
    // (intersection_position, intersection_normal)
    virtual Intersection intersect(const Ray &ray) = 0;
    virtual glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) = 0;
    virtual void transform( glm::mat4 transform ){}
    virtual std::vector<glm::vec3> get_points()
    {
        return std::vector<glm::vec3>();
    }
};

class GeometryLight : public Object
{
public:
    glm::vec3 position;

    virtual glm::vec3 sample_geometry() = 0;
    virtual std::vector<glm::vec3> get_edge_points() = 0;
    virtual Intersection intersect(const Ray &ray) = 0;
    virtual glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) = 0;
    
};

class SphereLight : public GeometryLight
{
private:
    glm::vec3 center;
    glm::vec3 color;
    std::vector<glm::vec3> edge_points;
    double radius;
    
public:
    SphereLight( glm::vec3 center, double radius, glm::vec3 color )
        : center(center), color(color), radius(radius)
    {
        this->position = center;
        edge_points.push_back(center);
    }

    glm::vec3 sample_geometry() override
    {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> sampler(0.0, 1.0);

        double x, y, z;
        double theta = 2 * M_PI * sampler(generator);
        double phi = acos(1 - 2 * sampler(generator));
        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);
        glm::vec3 sample(x, y, z); 
        sample = glm::normalize(sample) * (float)radius;

        return sample;
    }

    Intersection intersect(const Ray &ray) override
    {
        glm::vec3 intersect_position(0, 0, 0); 
        glm::vec3 intersect_normal(0, 0, 0); 

        bool intersects = glm::intersectRaySphere( ray.origin, 
                                                   ray.direction, 
                                                   center, 
                                                   radius,  
                                                   intersect_position, 
                                                   intersect_normal);
        auto intersection_dir = glm::normalize(intersect_position - ray.origin);
        auto ray_dir = glm::normalize(ray.direction - ray.origin);
        if ( intersects && glm::dot(intersection_dir, ray_dir) > 0 )
        {
            // TODO: (u,v)
            return Intersection(intersect_position, intersect_normal, glm::vec2(0, 0), this);
        }
        return Intersection();
    }

    std::vector<glm::vec3> get_edge_points()
    {
        return edge_points;
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return color;
    }
    
};


class Camera
{
public:
    glm::vec3 position;
    glm::vec3 look_at;
    glm::vec3 up;
    glm::vec3 horizontal, vertical, upper_left_corner;
    float aspect_ratio, viewport_height, viewport_width;

    Camera(glm::vec3 pos, glm::vec3 look, glm::vec3 up_dir, float ratio)
    {
        position = pos;
        look_at = look;
        up = glm::normalize(up_dir);
        aspect_ratio = ratio;
        viewport_height = 1.0; // might make this a parameter at some point...
        viewport_width = aspect_ratio * viewport_height;
        // temp variables
        glm::vec3 w = glm::normalize(position - look_at);
        glm::vec3 u = glm::normalize(glm::cross(up, w));
        glm::vec3 v = glm::cross(w, u);
        // arbitrary horizontal and vertical axes for our viewport
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        upper_left_corner = position - (0.5f * horizontal) + (0.5f * vertical) - w;
    }

    Ray spawn_ray(float x, float y)
    {
        return Ray(position, upper_left_corner + x*horizontal - y*vertical - position);
    }

    ~Camera(){}
};

class Scene
{
private: 
    std::vector<Object*> object_list;
    std::vector<Light*> light_list;
    std::vector<GeometryLight*> geometry_light_list;
    glm::vec3 sky_color;

public:
    Scene(glm::vec3 c)
    {
        sky_color = c;
    }

    ~Scene(){}

    void add_object(Object *object)
    {
        object_list.push_back( object );
    }

    void add_light(Light *light)
    {
        light_list.push_back( light );
    }

    void add_geometry_light(GeometryLight *light)
    {
        geometry_light_list.push_back( light );
    }

    void apply_transform( glm::mat4 transform )
    {
        for ( Object *o : object_list )
        {
            o->transform( transform );
        }
    }

    void construct_penumbra_wedges()
    {
        
        // for light in scene
        //   for object in scene
        //     for edge of light
        //       for edge of object
        //         cast ray from points on edge of light to points on edge of object
        //         see if they intersect another object
        //         if yes, construct penumbra triangle 
        /*
        for ( GeometryLight* light : geometry_light_list )
        {
            for ( Object* obj : object_list )
            {
                uint32_t i = 0, j = 1;
                auto light_points = light->get_edge_points();
                while ( i < light_points.size() )
                {
                    auto light_p1 = light_points[i];
                    auto light_p2 = light_points[j];
                    
                    uint32_t k = 0, l = 1;
                    auto obj_points = obj->get_points(); 
                    while ( k < obj_points.size() )
                    {
                        auto obj_p1 = obj_points[i];
                        auto obj_p2 = obj_points[j];
                        
                        
                        
                        k++;
                        l = (l + 1) % obj_points.size();
                    }
                    
                    i++;
                    j = (j + 1) % light_points.size();
                }
            }
        }*/
    }

    std::vector<std::pair<Object*, Intersection>> intersect_objects(Ray ray)
    {
        std::vector<std::pair<Object*, Intersection>> intersections;
        
        // TODO: spatial data structure
        for ( Object *obj : object_list )
        {
            auto intersection = obj->intersect( ray );
            if ( intersection.exists() )
            {
                intersections.push_back( std::make_pair(obj, intersection) );
            }
        }

        return intersections;
    }

    std::vector<std::pair<Object*, Intersection>> intersect_geometry(Ray ray)
    {
        std::vector<std::pair<Object*, Intersection>> intersections;
        
        // TODO: spatial data structure
        for ( Object *obj : object_list )
        {
            auto intersection = obj->intersect( ray );
            if ( intersection.exists() && glm::dot( intersection.position - ray.origin, ray.direction ) > 0 )
            {
                intersections.push_back( std::make_pair(obj, intersection) );
            }
        }

        for ( GeometryLight *g : geometry_light_list )
        {
            auto intersection = g->intersect( ray );
            if ( intersection.exists() )
            {
                intersections.push_back( std::make_pair((Object*)g, intersection) );
            }
        }
    
        return intersections;
    }

    std::vector<Object*> get_objects() 
    {
        return object_list;
    }

    std::vector<Light*> get_lights() 
    {
        return light_list;
    }

    std::vector<GeometryLight*> get_geometry_lights() 
    {
        return geometry_light_list;
    }
    
    glm::vec3 get_sky_color( int x, int y, int width, int height )
    {
        glm::vec3 color;
        color.x = (sky_color.x + ((255.0 - sky_color.x) * ((double)y / (double)height))) / 255.0;
        color.y = (sky_color.y + ((255.0 - sky_color.y) * ((double)y / (double)height))) / 255.0;
        color.z = (sky_color.z + ((255.0 - sky_color.z) * ((double)y / (double)height))) / 255.0;
        return color;
    }

};

class Material
{
public:
    virtual glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) = 0;
};

class BasicMaterial : public Material
{
private:
    glm::vec3 color;
public: 
    BasicMaterial( glm::vec3 color )
    {
        this->color = color;
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return color;
    }
};

class PhongMaterial : public Material
{
private:
    glm::vec3 ambient_color;
    glm::vec3 diffuse_color;
    glm::vec3 specular_color;
    Scene *scene;
    float ka = 1.0;
    float kd = 1.0;
    float ks = 1.0;
    float ke = 20.0;
    float kr, kt, refractive_index;
public: 
    PhongMaterial( glm::vec3 ambient_color, glm::vec3 diffuse_color, glm::vec3 specular_color, Scene *scene, float kr, float kt, float refractive_index ):
        ambient_color(ambient_color), diffuse_color(diffuse_color), specular_color(specular_color), scene(scene), kr(kr), kt(kt), refractive_index(refractive_index)
    {}

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        glm::vec3 color(0.0, 0.0, 0.0);
        auto lights = scene->get_lights();

        for ( Light *light : lights )
        {
            // spawn shadow rays 
            Ray shadow_ray(intersection.position, light->position - intersection.position);
            auto intersections = scene->intersect_objects( shadow_ray );
            bool lit = true;
            for ( auto inter : intersections )
            {
                auto obj = std::get<0>(inter);
                glm::vec3 obj_intersection = std::get<1>(inter).position;
                auto intersection_dir = obj_intersection - intersection.position;
                if ( obj != current_object && glm::dot(intersection_dir, shadow_ray.direction) > 0 && glm::distance(intersection.position, obj_intersection) < glm::distance(intersection.position, light->position) )
                {
                    lit = false;
                }
            }

            // color calculations
            if ( lit && inside == 0 )
            {
                glm::vec3 S = glm::normalize( light->position - intersection.position );
                glm::vec3 N = glm::normalize( intersection.normal );
                glm::vec3 R = glm::normalize( S - 2*(glm::dot(S,N))*N );
                glm::vec3 V = -glm::normalize( incoming_ray.origin - incoming_ray.direction );
                float lambertian = std::max(glm::dot(N, S), 0.0f);
                color += diffuse_color * lambertian * kd * light->intensity_to( intersection.position );
                if ( lambertian > 0.0f )
                {
                    color += specular_color * std::max((float)pow(glm::dot(R,V), ke), 0.0f) * ks * light->intensity_to( intersection.position );
                }
            }

            //color += ambient_color * ka;

            // color is direct illumination
            // so we need to calculate indirect illumination

            // REFLECTION
            glm::vec3 ind_color;
            if ( depth > 0 && kr > 0 )
            {
                glm::vec3 reflection_ray = glm::reflect(incoming_ray.direction, glm::faceforward(intersection.normal, intersection.normal, incoming_ray.direction));
                Ray outgoing_ray(intersection.position + (reflection_ray * 0.01f), reflection_ray);
                auto intersections = scene->intersect_geometry( outgoing_ray );
                if ( intersections.size() > 0 ) 
                {
                    Intersection closest_intersection;
                    double closest_intersect_dist = 0;
                    Object *closest_object = nullptr;
                    for ( uint32_t i = 0; i < intersections.size(); i++ )
                    {
                        auto temp_container = intersections[i];
                        auto object = std::get<0>(temp_container);
                        auto next_intersection = std::get<1>(temp_container);
                        if ( !closest_intersection.exists() || glm::distance(intersection.position, next_intersection.position) < closest_intersect_dist )
                        {
                            closest_intersection = next_intersection;
                            closest_intersect_dist = glm::distance(intersection.position, next_intersection.position);
                            closest_object = object;
                        }
                    }
                    if ( closest_object != nullptr )
                    { 
                        ind_color = closest_object->get_color(outgoing_ray, closest_intersection, closest_object, depth - 1, 0);
                    }
                    else 
                    {
                        ind_color = scene->get_sky_color(1, 1, 800, 600);
                    }
                }
                else 
                {
                    ind_color = scene->get_sky_color(1, 1, 800, 600);
                }
                color = color + kr*ind_color;
            }
            
            // TRANSMISSION
            
            if ( depth > 0 && kt > 0 )
            {
                glm::vec3 ind_color;
                glm::vec3 refraction_ray;
                if ( inside )
                {
                    refraction_ray = glm::refract(glm::normalize(incoming_ray.direction), glm::faceforward(intersection.normal, intersection.normal, incoming_ray.direction), refractive_index);
                }
                else 
                {
                    refraction_ray = glm::refract(glm::normalize(incoming_ray.direction), glm::faceforward(intersection.normal, intersection.normal, incoming_ray.direction), 1.0f/refractive_index);
                }
                /*if ( glm::length(refraction_ray) == 0 )
                {
                    // total internal reflection!
                    refraction_ray = glm::reflect(glm::normalize(incoming_ray.direction), glm::faceforward(intersection.normal, intersection.normal, incoming_ray.direction));
                }*/
                Ray outgoing_ray(intersection.position + (refraction_ray * 0.01f), glm::normalize(refraction_ray));
                auto intersections = scene->intersect_geometry( outgoing_ray );
                if ( intersections.size() > 0 ) 
                {
                    Intersection closest_intersection;
                    double closest_intersect_dist = 0;
                    Object *closest_object = nullptr;
                    for ( uint32_t i = 0; i < intersections.size(); i++ )
                    {
                        auto temp_container = intersections[i];
                        auto object = std::get<0>(temp_container);
                        auto next_intersection = std::get<1>(temp_container);
                        if ( !closest_intersection.exists() || glm::distance(intersection.position, next_intersection.position) < closest_intersect_dist )
                        {
                            closest_intersection = next_intersection;
                            closest_intersect_dist = glm::distance(intersection.position, next_intersection.position);
                            closest_object = object;
                        }
                    }
                    if ( closest_object == current_object )
                    {
                        ind_color = closest_object->get_color(outgoing_ray, closest_intersection, closest_object, depth - 1, 1);
                    }
                    else if ( closest_object != nullptr )
                    { 
                        ind_color = closest_object->get_color(outgoing_ray, closest_intersection, closest_object, depth - 1, 0);
                    }
                    else 
                    {
                        ind_color = scene->get_sky_color(1, 1, 800, 600);
                    }
                }
                else 
                {
                    ind_color = scene->get_sky_color(1, 1, 800, 600);
                }
                color = color + kt*ind_color;
            }


            return color;
        }

        // """ tone reproduction """
        // make sure our values are clamped to [0, 255]
        color.x = std::min(1.0f, color.x);
        color.y = std::min(1.0f, color.y);
        color.z = std::min(1.0f, color.z);

        color.x = std::max(0.0f, color.x);
        color.y = std::max(0.0f, color.y);
        color.z = std::max(0.0f, color.z);

        return color;
    }
}; 

class CheckeredPhongMaterial : public Material
{
    Scene *scene;
    glm::vec3 ambient_color;
    glm::vec3 diffuse_color_1;
    glm::vec3 diffuse_color_2;
    glm::vec3 specular_color;
    float ka = 1.0;
    float kd = 1.0;
    float ks = 1.0;
    float ke = 20.0;
    int checkernum;
public: 
    CheckeredPhongMaterial( glm::vec3 ambient, glm::vec3 diffuse_1, glm::vec3 diffuse_2, int checker_num, glm::vec3 specular, Scene *world )
    {
        ambient_color = ambient;
        this->diffuse_color_1 = diffuse_1;
        this->diffuse_color_2 = diffuse_2;
        specular_color = specular;
        scene = world;
        this->checkernum = checker_num;
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        glm::vec3 color(0.0, 0.0, 0.0);
        auto lights = scene->get_lights();

        for ( Light *light : lights )
        {
            // spawn shadow rays 
            Ray shadow_ray(intersection.position, light->position - intersection.position);
            auto intersections = scene->intersect_objects( shadow_ray );
            bool lit = true;
            for ( auto inter : intersections )
            {
                auto obj = std::get<0>(inter);
                glm::vec3 obj_intersection = std::get<1>(inter).position;
                auto intersection_dir = obj_intersection - intersection.position;
                if ( obj != current_object && glm::dot(intersection_dir, shadow_ray.direction) > 0 && glm::distance(intersection.position, obj_intersection) < glm::distance(intersection.position, light->position) )
                {
                    lit = false;
                }
            }

            // color calculations
            if ( lit )
            {
                glm::vec3 S = glm::normalize( light->position - intersection.position );
                glm::vec3 N = glm::normalize( intersection.normal );
                glm::vec3 R = glm::normalize( S - 2*(glm::dot(S,N))*N );
                glm::vec3 V = -glm::normalize( incoming_ray.origin - incoming_ray.direction );
                float lambertian = std::max(glm::dot(N, S), 0.0f);
                auto u = intersection.uv.x;
                auto v = intersection.uv.y;
                if ( (fmod(u * checkernum, 1.0)) < 0.5 )
                {
                    if ( (fmod(v * checkernum, 1.0)) < 0.5 )
                    {
                        color += diffuse_color_1 * lambertian * kd * light->intensity_to( intersection.position );
                    }
                    else
                    {
                        color += diffuse_color_2 * lambertian * kd * light->intensity_to( intersection.position );
                    }
                }
                else
                {
                    if ( (fmod(v * checkernum, 1.0)) < 0.5 )
                    {
                        color += diffuse_color_2 * lambertian * kd * light->intensity_to( intersection.position );
                    }
                    else
                    {
                        color += diffuse_color_1 * lambertian * kd * light->intensity_to( intersection.position );
                    }
                }
                if ( lambertian > 0.0f )
                {
                    color += specular_color * std::max((float)pow(glm::dot(R,V), ke), 0.0f) * ks * light->intensity_to( intersection.position );
                }
            }
        }

        // """ tone reproduction """
        // make sure our values are clamped to [0, 255]
        color.x = std::min(1.0f, color.x);
        color.y = std::min(1.0f, color.y);
        color.z = std::min(1.0f, color.z);

        color.x = std::max(0.0f, color.x);
        color.y = std::max(0.0f, color.y);
        color.z = std::max(0.0f, color.z);

        return color;
    }
};



class Triangle : public Object
{
public:
    glm::vec3 a, b, c, normal;
    glm::vec2 a_uv, b_uv, c_uv;

    Triangle(){}

    Triangle(glm::vec3 A, glm::vec3 B, glm::vec3 C, Material *mat, glm::vec2 a_uv, glm::vec2 b_uv, glm::vec2 c_uv )
    {
        a = A;
        b = B;
        c = C;
        normal = glm::triangleNormal(a, b, c); // TODO: make sure this normal is correct
        material = mat;
        this->a_uv = a_uv;
        this->b_uv = b_uv;
        this->c_uv = c_uv;
    }

    Intersection intersect(const Ray &ray) override
    {
        glm::vec2 bary_coords(0, 0); 
        float distance = 0;

        bool intersects = glm::intersectRayTriangle( ray.origin, 
                                                     ray.direction, 
                                                     a, 
                                                     b, 
                                                     c, 
                                                     bary_coords, // if triangle is vecX, this must be vec(X-1)
                                                     distance // this parameter is undocumented in the default Google result for GLM! 
                                                     );

        if ( intersects )
        {
            // compute barycentric interpolation
            float u,v;
            u = bary_coords.x;
            v = bary_coords.y;

            glm::vec3 intersect_position( (1-u-v)*a + u*b + v*c );
            glm::vec2 uv_interp( (1-u-v)*a_uv + u*b_uv + v*c_uv  );
            if ( uv_interp.x > 1 || uv_interp.y > 1 ) 
            {
                std::cout << "(" << uv_interp.x << ", " << uv_interp.y << ")"  << std::endl;
            }

            return Intersection(intersect_position, normal, uv_interp, this);
        }
        return Intersection();
    }

    void transform( glm::mat4 transform ) override
    {
        a = glm::vec3( transform * glm::vec4(a, 1.0) );
        b = glm::vec3( transform * glm::vec4(b, 1.0) );
        c = glm::vec3( transform * glm::vec4(c, 1.0) );
        normal = glm::vec3( transform * glm::vec4(normal, 1.0) );
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return material->get_color( incoming_ray, intersection, current_object, depth, 0 );
    }

    std::vector<glm::vec3> get_points()
    {
        std::vector<glm::vec3> points;
        points.push_back(a);
        points.push_back(b);
        points.push_back(c);
        return points;
    }
};

class BoxLight : public GeometryLight
{
private:
    std::vector<Triangle> polygons;
    glm::vec3 color;
    glm::vec3 top_left, top_right, bottom_left, bottom_right;
    std::vector<glm::vec3> edge_points;
    Triangle t1, t2;
    
public:
    BoxLight( glm::vec3 top_left, glm::vec3 top_right, glm::vec3 bottom_left, glm::vec3 bottom_right, glm::vec3 color )
        :top_left(top_left), top_right(top_right), bottom_left(bottom_left), bottom_right(bottom_right)
    {
        BasicMaterial temp(color);
        t1 = Triangle(top_left, bottom_right, bottom_left, &temp, glm::vec2(), glm::vec2(), glm::vec2());
        t2 = Triangle(top_left, top_right, bottom_right, &temp, glm::vec2(), glm::vec2(), glm::vec2());
        polygons.push_back(t1);
        polygons.push_back(t2);
        /*edge_points.push_back(top_left);
        edge_points.push_back(top_right);
        edge_points.push_back(bottom_left);
        edge_points.push_back(bottom_right);
        */
        float num = 32.0f;
        for ( float i = 0.0f; i <= num; i += 1.0f )
        {
            glm::vec3 cur_left = lerp(top_left, bottom_left, i / num);
            glm::vec3 cur_right = lerp(top_right, bottom_right, i / num);;
            for ( float j = 0.0f; j <= num; j += 1.0f )
            {
                edge_points.push_back( lerp(cur_left, cur_right, j / num) );
            }
        }
        this->color = color;
    }

    glm::vec3 lerp(glm::vec3 a, glm::vec3 b, float t)
    {
        return a + ((b - a) * t);
    }

    glm::vec3 sample_geometry() override
    {
        float x, y;
        x = ((float) rand() / (RAND_MAX));
        y = ((float) rand() / (RAND_MAX));
        glm::vec3 y_vec_1 = lerp(bottom_left, top_left, y); // lerp left side vertically by y sample
        glm::vec3 y_vec_2 = lerp(bottom_right, top_right, y); // lerp right side vertically by y sample
        return lerp(y_vec_1, y_vec_2, x); // lerp between two vertical lerps by x sample
    }

    Intersection intersect(const Ray &ray) override
    {
        auto intersection = t1.intersect(ray);
        if ( intersection.exists() )
        {
            return intersection;
        }
        intersection = t2.intersect(ray);    
        if ( intersection.exists() )
        {
            return intersection;
        }        
        return Intersection();
    }

    std::vector<glm::vec3> get_edge_points()
    {
        return edge_points;
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return color;
    }
    
};


class SoftPhongMaterial : public Material
{
private:
    Scene *scene;
    glm::vec3 ambient_color;
    glm::vec3 diffuse_color;
    glm::vec3 specular_color;
    float ka = 1.0;
    float kd = 1.0;
    float ks = 1.0;
    float ke = 20.0;
    uint16_t num_samples = 10;

public: 
    SoftPhongMaterial( glm::vec3 ambient, glm::vec3 diffuse, glm::vec3 specular, Scene *world )
    {
        ambient_color = ambient;
        diffuse_color = diffuse;
        specular_color = specular;
        scene = world;
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        glm::vec3 color(0.0, 0.0, 0.0);
        auto lights = scene->get_geometry_lights();
        for ( GeometryLight *light : lights )
        {
            // MONTE CARLO
            /*
            for ( uint16_t i = 0; i < num_samples; i++ )
            {
                // spawn shadow rays 
                Ray shadow_ray(intersection.position, light->sample_geometry() - intersection.position);
                auto intersections = scene->intersect_objects( shadow_ray );
                bool lit = true;
                for ( auto inter : intersections )
                {
                    auto obj = std::get<0>(inter);
                    glm::vec3 obj_intersection = std::get<1>(inter).position;
                    auto intersection_dir = obj_intersection - intersection.position;
                    if ( obj != current_object && glm::dot(intersection_dir, shadow_ray.direction) > 0 && glm::distance(intersection.position, obj_intersection) < glm::distance(intersection.position, light->position) )
                    {
                        lit = false;
                    }
                }

                // color calculations
                if ( true )
                {
                    glm::vec3 S = glm::normalize( light->position - intersection.position );
                    glm::vec3 N = glm::normalize( intersection.normal );
                    glm::vec3 R = glm::normalize( S - 2*(glm::dot(S,N))*N );
                    glm::vec3 V = -glm::normalize( incoming_ray.origin - incoming_ray.direction );
                    float lambertian = std::max(glm::dot(N, S), 0.0f);
                    color += diffuse_color * lambertian * kd * (1.0f/num_samples);// * ((float)num_in_shade/(float)points.size());
                    if ( lambertian > 0.0f )
                    {
                        color += specular_color * std::max((float)pow(glm::dot(R,V), ke), 0.0f) * ks * (1.0f/num_samples);// * ((float)num_in_shade/(float)points.size());
                    }
                }
                color += ambient_color * ka * (1.0f/num_samples);
            }
            */
            // NOT MONTE CARLO
            uint64_t num_in_shade = 0;
            auto points = light->get_edge_points();
            for ( glm::vec3 point : points )
            {
                Ray shadow_ray(intersection.position, point - intersection.position);
                auto intersections = scene->intersect_objects( shadow_ray );
                bool lit = true;
                for ( auto inter : intersections )
                {
                    auto obj = std::get<0>(inter);
                    glm::vec3 obj_intersection = std::get<1>(inter).position;
                    auto intersection_dir = obj_intersection - intersection.position;
                    if ( obj != current_object && glm::dot(intersection_dir, shadow_ray.direction) > 0 && glm::distance(intersection.position, obj_intersection) < glm::distance(intersection.position, light->position) )
                    {
                        lit = false;
                    }
                }
                if ( lit == false )
                {
                    num_in_shade++;
                }
            }

            // color calculations
            if ( num_in_shade < points.size() )
            {
                glm::vec3 S = glm::normalize( light->position - intersection.position );
                glm::vec3 N = glm::normalize( intersection.normal );
                glm::vec3 R = glm::normalize( S - 2*(glm::dot(S,N))*N );
                glm::vec3 V = -glm::normalize( incoming_ray.origin - incoming_ray.direction );
                float lambertian = std::max(glm::dot(N, S), 0.0f);
                color += diffuse_color * lambertian * kd * ((float)(points.size() - num_in_shade)/(float)points.size());
                if ( lambertian > 0.0f )
                {
                    color += specular_color * std::max((float)pow(glm::dot(R,V), ke), 0.0f) * ks * ((float)(points.size() - num_in_shade)/(float)points.size());
                }
            }
                color += ambient_color * ka * (1.0f/num_samples);
        }

        // """ tone reproduction """
        // make sure our values are clamped to [0, 255]
        color.x = std::min(1.0f, color.x);
        color.y = std::min(1.0f, color.y);
        color.z = std::min(1.0f, color.z);

        color.x = std::max(0.0f, color.x);
        color.y = std::max(0.0f, color.y);
        color.z = std::max(0.0f, color.z);

        return color;
    }
}; 

class Sphere : public Object 
{
private:
    glm::vec3 center;
    double radius;
    // material is inherited from Object

public:
    Sphere(glm::vec3 point, double r, Material *mat)
    {
        center = point;
        radius = r;
        material = mat;
    }

    Sphere(glm::vec3 point, double r, Material *mat, std::string name)
    {
        center = point;
        radius = r;
        material = mat;
        this->name = name;
    }

    Intersection intersect(const Ray &ray) override
    {
        glm::vec3 intersect_position(0, 0, 0); 
        glm::vec3 intersect_normal(0, 0, 0); 

        bool intersects = glm::intersectRaySphere( ray.origin, 
                                                   ray.direction, 
                                                   center, 
                                                   radius,  
                                                   intersect_position, 
                                                   intersect_normal);
        auto intersection_dir = glm::normalize(intersect_position - ray.origin);
        auto ray_dir = glm::normalize(ray.direction - ray.origin);
        if ( intersects && glm::dot(intersection_dir, ray_dir) > 0 )
        {
            // TODO: (u,v)
            return Intersection(intersect_position, intersect_normal, glm::vec2(0,0), this);
        }
        return Intersection();
    }   

    void transform( glm::mat4 transform ) override
    {
        center = glm::vec3( transform * glm::vec4(center, 1.0) );
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return material->get_color( incoming_ray, intersection, current_object, depth, 0 );
    }

    ~Sphere(){}
};


class Floor : public Object
{
private:
    std::vector<Triangle> polygons;

public:
    Floor(glm::vec3 top_left, glm::vec3 top_right, glm::vec3 bottom_left, glm::vec3 bottom_right, Material *mat)
    {
        material = mat;
        glm::vec2 tl_uv(0, 1);
        glm::vec2 tr_uv(1, 1);
        glm::vec2 bl_uv(0, 0);
        glm::vec2 br_uv(1, 0);
        Triangle t1( bottom_left, top_right, top_left, mat, bl_uv, tr_uv, tl_uv );
	    Triangle t2( bottom_left, bottom_right, top_right, mat, bl_uv, br_uv, tr_uv );
        polygons.push_back( t1 );
        polygons.push_back( t2 );
        this->name = "floor";
    }

    Intersection intersect(const Ray &ray) override
    {
        for ( Triangle t : polygons )
        {
            auto intersection = t.intersect(ray);
            if ( intersection.exists() )
            {
                return intersection;
            }
        }        
        return Intersection();
    }

    void transform( glm::mat4 transform ) override
    {
        for ( Triangle t : polygons )
        {
            t.transform( transform );
        }
    }

    glm::vec3 get_color( Ray incoming_ray, Intersection intersection, Object *current_object, int depth, int inside ) override
    {
        return material->get_color( incoming_ray, intersection, current_object, depth, 0 );
    }

    ~Floor(){}
    
};

#endif