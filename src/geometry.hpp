#ifndef GEOMETRY
#define GEOMETRY

#include <optional>
#include <vector>
#include <utility>
#include <iostream>
#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtx/normal.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/gtx/vec_swizzle.hpp> 


void vec3_print( glm::vec3 vector )
{
    std::cout << "[" << vector.x << ", " << vector.y << ", " << vector.z << "]" << std::endl;
}

class Ray
{
public:
    glm::vec3 origin;
    glm::vec3 direction;

    Ray(glm::vec3 o, glm::vec3 d)
    {
        origin = o;
        direction = glm::normalize(d);
    }

    ~Ray(){}

};


class Object
{
public:
    glm::vec3 material;

    // Intersect should return vectors of the form:
    // (intersection_position, intersection_normal)
    virtual std::optional<std::pair<glm::vec3, glm::vec3>> intersect(const Ray &ray)
    {
        return std::nullopt;
    }

    virtual void transform( glm::mat4 transform ){}
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
    glm::vec3 sky_color;

public:
    Scene(glm::vec3 c)
    {
        sky_color = c;
    }

    ~Scene(){}

    void add_object(Object* object)
    {
        object_list.push_back( object );
    }

    void apply_transform( glm::mat4 transform )
    {
        for ( Object *o : object_list )
        {
            o->transform( transform );
        }
    }

    std::vector<Object*> get_objects() // TODO: spatial data structure
    {
        return object_list;
    }

};

class Sphere : public Object 
{
private:
    glm::vec3 center;
    double radius;
    // material is inherited from Object

public:
    Sphere(glm::vec3 point, double r, glm::vec3 color)
    {
        center = point;
        radius = r;
        material = color;
    }

    std::optional<std::pair<glm::vec3, glm::vec3>> intersect(const Ray &ray) override
    {
        glm::vec3 intersect_position(0, 0, 0); 
        glm::vec3 intersect_normal(0, 0, 0); 

        bool intersects = glm::intersectRaySphere( ray.origin, 
                                                   ray.direction, 
                                                   center, 
                                                   radius,  
                                                   intersect_position, 
                                                   intersect_normal);
        
        if ( intersects )
        {
            return std::optional<std::pair<glm::vec3, glm::vec3>>{ std::pair<glm::vec3, glm::vec3>{intersect_position, intersect_normal} };
        }
        return std::nullopt;
    }   

    void transform( glm::mat4 transform ) override
    {
        center = glm::vec3( transform * glm::vec4(center, 1.0) );
    }

    ~Sphere(){}
};

class Triangle : public Object
{
public:
    glm::vec3 a, b, c, normal;

    Triangle(glm::vec3 A, glm::vec3 B, glm::vec3 C, glm::vec3 mat)
    {
        a = A;
        b = B;
        c = C;
        normal = glm::triangleNormal(a, b, c); // TODO: make sure this normal is correct
        material = mat;
    }

    std::optional<std::pair<glm::vec3, glm::vec3>> intersect(const Ray &ray) override
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

            return std::optional<std::pair<glm::vec3, glm::vec3>>{ std::pair<glm::vec3, glm::vec3>{intersect_position, normal} };
        }
        return std::nullopt;
    }

    void transform( glm::mat4 transform ) override
    {
        a = glm::vec3( transform * glm::vec4(a, 1.0) );
        b = glm::vec3( transform * glm::vec4(b, 1.0) );
        c = glm::vec3( transform * glm::vec4(c, 1.0) );
        normal = glm::vec3( transform * glm::vec4(normal, 1.0) );
    }

    ~Triangle(){}
};

class Floor : public Object
{
private:
    std::vector<Triangle> polygons;

public:
    Floor(glm::vec3 top_left, glm::vec3 top_right, glm::vec3 bottom_left, glm::vec3 bottom_right, glm::vec3 color)
    {
        material = color;
        Triangle t1(top_left, top_right, bottom_right, material);
        Triangle t2(top_left, bottom_right, bottom_left,  material);
        polygons.push_back( t1 );
        polygons.push_back( t2 );
    }

    std::optional<std::pair<glm::vec3, glm::vec3>> intersect(const Ray &ray) override
    {
        std::optional<std::pair<glm::vec3, glm::vec3>> intersects;
        for ( Triangle t : polygons )
        {
            intersects = t.intersect(ray);
            if ( intersects.has_value() )
            {
                return intersects;
            }
        }        
        return std::nullopt;
    }

    void transform( glm::mat4 transform ) override
    {
        for ( Triangle t : polygons )
        {
            t.transform( transform );
        }
    }

    ~Floor(){}
    
};

#endif