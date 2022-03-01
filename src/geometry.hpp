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

class Scene;
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

    ~Ray(){}

};

class Intersection
{
private:
    bool has_value;
public:
    glm::vec3 position;
    glm::vec3 normal;

    bool exists()
    {
        return has_value;
    }

    Intersection()
    {
        has_value = false;
    }

    Intersection(glm::vec3 position, glm::vec3 normal)
    {
        has_value = true;
        this->position = position;
        this->normal = normal;
    }
};

class Light
{
public:
    glm::vec3 color;
    glm::vec3 position;
    
    Light( glm::vec3 color, glm::vec3 position )
    {
        this->color = color;
        this-> position = position;
    }
};


class Object
{
public:
    Material *material;
    std::string name;

    // Intersect should return vectors of the form:
    // (intersection_position, intersection_normal)
    virtual Intersection intersect(const Ray &ray)
    {
        std::runtime_error("Should never call the Intersect function of the base material class!\n");
        return Intersection();
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
    std::vector<Light*> light_list;
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

    void apply_transform( glm::mat4 transform )
    {
        for ( Object *o : object_list )
        {
            o->transform( transform );
        }
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

    std::vector<Object*> get_objects() 
    {
        return object_list;
    }

    std::vector<Light*> get_lights() 
    {
        return light_list;
    }

};

class Material
{
public:
    virtual glm::vec3 get_color( Ray incoming_ray, glm::vec3 intersection, glm::vec3 normal, Object *current_object )
    {
        return glm::vec3(0.0, 0.0, 0.0);
    }
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

    glm::vec3 get_color( Ray incoming_ray, glm::vec3 intersection, glm::vec3 normal, Object *current_object ) override
    {
        return color;
    }
};

class PhongMaterial : public Material
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
public: 
    PhongMaterial( glm::vec3 ambient, glm::vec3 diffuse, glm::vec3 specular, Scene *world )
    {
        ambient_color = ambient;
        diffuse_color = diffuse;
        specular_color = specular;
        scene = world;
    }

    glm::vec3 get_color( Ray incoming_ray, glm::vec3 intersection, glm::vec3 normal, Object *current_object ) override
    {
        glm::vec3 color(0.0, 0.0, 0.0);
        auto lights = scene->get_lights();

        for ( Light *light : lights )
        {
            // spawn shadow rays 
            Ray shadow_ray(intersection, light->position - intersection);
            auto intersections = scene->intersect_objects( shadow_ray );
            bool lit = true;
            for ( auto inter : intersections )
            {
                auto obj = std::get<0>(inter);
                glm::vec3 obj_intersection = std::get<1>(inter).position;
                auto intersection_dir = obj_intersection - intersection;
                if ( obj != current_object && glm::dot(intersection_dir, shadow_ray.direction) > 0 )
                {
                    lit = false;
                }
            }
            // color calculations
            if ( lit )
            {
                glm::vec3 S = glm::normalize( light->position - intersection );
                glm::vec3 N = glm::normalize( normal );
                glm::vec3 R = glm::normalize( S - 2*(glm::dot(S,N))*N );
                glm::vec3 V = -glm::normalize( incoming_ray.origin - incoming_ray.direction );
                float lambertian = std::max(glm::dot(N, S), 0.0f);
                color += diffuse_color * lambertian * kd;
                if ( lambertian > 0.0f )
                {
                    color += specular_color * std::max((float)pow(glm::dot(R,V), ke), 0.0f) * ks;
                }
            }
            color += ambient_color * ka;
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
        
        if ( intersects )
        {
            return Intersection(intersect_position, intersect_normal);
        }
        return Intersection();
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

    Triangle(glm::vec3 A, glm::vec3 B, glm::vec3 C, Material *mat)
    {
        a = A;
        b = B;
        c = C;
        normal = glm::triangleNormal(a, b, c); // TODO: make sure this normal is correct
        material = mat;
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

            return Intersection(intersect_position, normal);
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

    ~Triangle(){}
};

class Floor : public Object
{
private:
    std::vector<Triangle> polygons;

public:
    Floor(glm::vec3 top_left, glm::vec3 top_right, glm::vec3 bottom_left, glm::vec3 bottom_right, Material *mat)
    {
        material = mat;
        Triangle t1(top_left, top_right, bottom_right, material);
        Triangle t2(top_left, bottom_right, bottom_left,  material);
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

    ~Floor(){}
    
};

#endif