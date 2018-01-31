//TODO:
    //ray-plane shading isn't working properly - weird backwards inputs
    //shadow rays
//OPT-TODO:
    //cornell box (ray-quad)
    //ray-triangle & .obj files
    //reflections
    //refractions
    //textures
    //area light
    //anti-aliasing

#include "OpenGP/Image/Image.h"
#include "bmpwrite.h"

using namespace OpenGP;

using Colour = Vec3; // RGB Value
Colour red() { return Colour(1.0f, 0.0f, 0.0f); }
Colour white() { return Colour(1.0f, 1.0f, 1.0f); }
Colour black() { return Colour(0.0f, 0.0f, 0.0f); }

///Class Surface used for all objects, members = position, ambient/diffuse/specular lighting coefficients, Phong exponent
class Surface {
public:
    Vec3 pos;
    Colour ambient;
    Colour diffuse;
    Colour specular;
    float phongexp;
    Surface(Vec3 p,Colour a,Colour d,Colour s,float x) : pos(p), ambient(a), diffuse(d), specular(s), phongexp(x) {}
};

///Class Plane inherits Surface, members = all Surface members plus plane normal
class Plane : public Surface {
public:
    Vec3 normal;
    float distance;
    Vec3 intersect;
    Plane(Vec3 p,Colour a,Colour d,Colour s,float x, Vec3 n) : Surface(p,a,d,s,x), normal(n.normalized()) {}
    bool intersects(Vec3 ray, Vec3 E) {
        float d = -1.0f;
        float denom = ray.dot(normal);  //denominator of ray-plane intersection equation
        if(denom > 0.00000000001f){     //we can't calculate if denom is zero
            Vec3 emp = E - pos;         //vector from camera to intersection point
            d = -emp.dot(normal);
        }
        if(d >= 0){
            distance = d;
            intersect = E + distance*ray;
            return true;
        }
        return false;
    }
};

///Class Sphere inherits Surface, members = all Surface members plus sphere radius
class Sphere : public Surface {
public:
    float radius;
    float distance;
    Vec3 intersect;
    Sphere(Vec3 p,Colour a,Colour d,Colour s,float x, float r) : Surface(p,a,d,s,x), radius(r) {}
    bool intersects(Vec3 ray, Vec3 E){
        Vec3 ems = E - pos;
        float disc = std::powf(ray.dot(ems),2) - ems.dot(ems) + radius*radius;
        if(disc >= 0){
            distance = -(ray.dot(ems)) - std::sqrtf(disc);   //find distance where ray hits sphere
            intersect = E + distance*ray;
            return true;
        }
        return false;
    }
    Vec3 findNormal(){
        return (intersect - pos)/radius;
    }

};

///Returns a ray from a light source to a given position
Vec3 lightDirection(Vec3 lightpos, Vec3 pos){
    return lightpos - pos;
}

///Computes the ambient light for a pixel based on ambient coefficient of a surface & ambient light intensity
Colour ambientShading(Colour ka, float ia){
    return ka*ia;
}

///Computes the colour of a pixel based on the given ray, normal to surface, vector pointing to light source, and surface/light coefficients/Phong exponent
Colour phongShading(Vec3 v, Vec3 n, Vec3 l, Colour ka, float ia, Colour kd, float i, Colour ks, float p){
    ///normalize all vectors just in case
    v = -v;
    v = v.normalized();
    n = n.normalized();
    l = l.normalized();

    ///compute ambient term
    Colour ambient = ambientShading(ka,ia);

    ///compute diffuse term
    Colour diffuse = kd*i*std::fmaxf(l.dot(n),0.0f);

    ///compute specular term
    Vec3 h = (v+l).normalized();
    Colour specular = ks*i*std::powf(std::fmaxf(h.dot(n),0.0f),p);

    ///return complete shade
    return ambient + diffuse + specular;
}

///Vec3 point = (left->right(x), up->down(y), far->near(z))
///Vec3 vector = (right(x),up(x),near(x))
int main(int, char**){

    /// #rows = hResolution, #cols = wResolution
    int wresolution = 1080;
    int hresolution = 720;
    float aratio = float(wresolution)/float(hresolution);
    Image<Colour> image(hresolution, wresolution);

    ///view bounds:
    float left = -1.0f*aratio;      //leftmost view plane boundary
    float right = 1.0f*aratio;      //rightmost
    float top = 1.0f;               //top
    float bottom = -1.0f;           //bottom
    ///without aratio, image will be distorted: our bounds define a square but the image aspect ratio is 4:3

    ///camera position definition:
    Vec3 W = Vec3(0.0f,0.0f,-1.0f); //out direction
    Vec3 V = Vec3(0.0f,1.0f,0.0f);  //up direction
    Vec3 U = Vec3(1.0f,0.0f,0.0f);  //side direction
    float d = 1.0f;                 //focal length
    Vec3 E = -d*W;                  //camera position

    ///light source & intensity:
    Vec3 lightpos = Vec3(0.0f, 0.0f, 0.0f);
    float lightintensity = 1.0f;

    ///Ambient light intensity:
    float ambientlight = 0.5f;

    ///floor position & material definition:
    Vec3 floorp = Vec3(0.0f,-50.0f,0.0f);
    Colour floora = Colour(0.5f,0.5f,0.5f);
    Colour floord = Colour(0.0f,0.2f,0.0f);
    Colour floors = white();
    float floorx = 10.0f;
    Vec3 floorn = Vec3(0.0f,-1.0f,0.0f);
    ///floor construction:
    Plane floor(floorp,floora,floord,floors,floorx,floorn);

    ///sphere1 position & material definition:
    Vec3 sphere1p = Vec3(3.0f, 1.0f, -5.0f);
    Colour sphere1a = Colour(0.5f,0.5f,0.5f);
    Colour sphere1d = red();
    Colour sphere1s = white();
    float sphere1x = 10.0f;
    float sphere1r = 3.0f;
    ///sphere1 construction:
    Sphere sphere1(sphere1p,sphere1a,sphere1d,sphere1s,sphere1x,sphere1r);

    ///For each row, compute rays, intersections, & shading
    for (int row = 0; row < image.rows(); ++row) {
        for (int col = 0; col < image.cols(); ++col) {

            /// (row,column) = current point on grid
            Vec3 pixel = left*U + (col*(right-left)/image.cols())*U;    //x-coordinate
            pixel += bottom*V + (row*(top-bottom)/image.rows())*V;      //y-coordinate

            ///build primary rays
            Vec3 ray = pixel - E;       //coordinate minus camera position; standard vector equation
            ray = ray.normalized();     //unit vector of ray

            /// structure:
            /// bool hit = false
            /// for each object in scene from closest to farthest:
            ///     if ray intersects object:
            ///         hit = true
            ///         intersection normal = normal from intersection with object
            ///         break out of for loop
            /// if hit:
            ///     shoot shadow ray at light source
            ///     if shadow ray hits any object:
            ///         image(row,col) = ambientShading
            ///     else:
            ///         image(row,col) = phongShading
            /// else:
            ///     image(row,col) = white();

            bool hit = false;
            Vec3 normal,lightdir,ambient,diffuse,specular;
            float phongexp;
            if(sphere1.intersects(ray,E)){
                hit = true;
                normal = sphere1.findNormal();
                lightdir = lightDirection(lightpos,sphere1.intersect);
                ambient = sphere1.ambient;
                diffuse = sphere1.diffuse;
                specular = sphere1.specular;
                phongexp = sphere1.phongexp;
                //image(row,col) = phongShading(ray, normal, lightdir, sphere1.ambient, ambientlight, sphere1.diffuse, lightintensity, sphere1.specular, sphere1.phongexp);
            }else if(floor.intersects(ray,E)){
                hit = true;
                normal = floor.normal;
                lightdir = lightDirection(floor.intersect,lightpos);   //TODO: this is backwards for some reason (floor.intersect is negative?)
                ambient = floor.ambient;
                diffuse = floor.diffuse;
                specular = floor.specular;
                phongexp = floor.phongexp;
                //image(row,col) = phongShading(ray, normal, lightdir, floor.ambient, ambientlight, floor.diffuse, lightintensity, floor.specular, floor.phongexp);
            }
            if(hit){
                /// shoot ray from (row,column) back towards light source
                Vec3 shadowray;
                //if obstructed:
                    //image(row,col) = ambientShading(ambient, ambientlight);
                //else:
                image(row,col) = phongShading(ray, normal, lightdir, ambient, ambientlight, diffuse, lightintensity, specular, phongexp);
            }else{
                image(row,col) = white();
            }

            /// TODO shadow rays (not here probably):
                /// shoot ray from (row,column) back towards light source
                /// if it's obstructed, only use ambient light
                /// else use full Phong shading

       }
    }

    bmpwrite("../../out.bmp", image);
    imshow(image);

    return EXIT_SUCCESS;
}
