//TODO:
    //ray-plane shading isn't working properly - weird backwards inputs
    //shadow rays from sphere always intersect sphere, with or without epsilon of any size
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
Colour green() { return Colour(0.0f, 1.0f, 0.0f); }
Colour blue() { return Colour(0.0f, 0.0f, 1.0f); }
Colour yellow() { return Colour(1.0f, 1.0f, 0.0f); }
Colour turquoise() { return Colour(0.0f, 1.0f, 1.0f); }
Colour purple() { return Colour(1.0f, 0.0f, 1.0f); }
Colour white() { return Colour(1.0f, 1.0f, 1.0f); }
Colour grey() { return Colour(0.5f, 0.5f, 0.5f); }
Colour black() { return Colour(0.0f, 0.0f, 0.0f); }

float epsilon = 0.0001f;
/*
///returns area of triangle given by 3 points
float triangleArea(Vec3 v1, Vec3 v2, Vec3 v3){
    return 0.5f*(v1(0)*v2(1) + v2(0)*v3(1) + v3(0)*v1(1) - v1(0)*v3(1) - v2(0)*v1(1) - v3(0)*v2(1));
}

Vec3 findCentroid(Vec3 v1, Vec3 v2, Vec3 v3){
    float x = (v1(0)+v2(0)+v3(0))/3;
    float y = (v1(1)+v2(1)+v3(1))/3;
    float z = (v1(2)+v2(2)+v3(2))/3;
    return Vec3(x,y,z);
}
*/

///Class Surface used for all objects, members = position, ambient/diffuse/specular lighting coefficients, Phong exponent
class Surface {
public:
    Vec3 pos;
    Colour ambient;
    Colour diffuse;
    Colour specular;
    float phongexp;
    Surface(Vec3 p,Colour a,Colour d,Colour s,float x) : pos(p), ambient(a), diffuse(d), specular(s), phongexp(x) {}
    void intersects() {}
    void findNormal() {}
};

///Class Plane inherits Surface, members = all Surface members plus plane normal
class Plane : public Surface {
public:
    Vec3 normal;
    float distance;
    Vec3 intersect;
    Plane(Vec3 p,Colour a,Colour d,Colour s,float x, Vec3 n) : Surface(p,a,d,s,x), normal(n.normalized()) {}
    bool intersects(Vec3 ray, Vec3 origin) {
        float d;
        float denom = normal.dot(ray);  //denominator of ray-plane intersection equation
        if(std::abs(denom) > epsilon){     //we can't calculate if denom is zero
            Vec3 pmo = (pos-origin);         //vector from camera to intersection point
            d = pmo.dot(normal)/denom;
        }
        if(d > epsilon){
            distance = d;
            intersect = origin + distance*ray;
            return true;
        }
        return false;
    }
    Vec3 findNormal(){
        return normal;
    }
};

///Class Sphere inherits Surface, members = all Surface members plus sphere radius
class Sphere : public Surface {
public:
    float radius;
    float distance;
    Vec3 intersect;
    Sphere(Vec3 p,Colour a,Colour d,Colour s,float x,float r) : Surface(p,a,d,s,x), radius(r) {}
    bool intersects(Vec3 ray, Vec3 origin){
        Vec3 pmo = pos - origin;  //fix these names
        float tca = pmo.dot(ray);
        float disc = pmo.dot(pmo) - tca*tca;
        if (disc > radius*radius) return false;
        float thc = std::sqrtf(radius*radius - disc);
        float d0 = tca - thc;
        float d1 = tca + thc;
        if(d0 > d1) std::swap(d0,d1);
        if(d0 < 0){
            d0 = d1;
            if(d0 < 0) return false;
        }

        distance = d0;
        intersect = origin + distance*ray;
        return true;
    }
    Vec3 findNormal(){
        return (intersect - pos)/radius;
    }
};
/*
///Class Triangle inherits Surface, members = all Surface members plus 3 triangle vertices
class Triangle : public Surface {
public:
    Vec3 vertex1;
    Vec3 vertex2;
    Vec3 vertex3;
    Vec3 normal;
    float distance;
    float area;
    Vec3 intersect;
    Triangle(Vec3 p,Colour a,Colour d,Colour s,float x,Vec3 v1,Vec3 v2,Vec3 v3) : Surface(p,a,d,s,x), vertex1(v1), vertex2(v2), vertex3(v3), normal((v2-v1).cross(v3-v1)), area(triangleArea(v1,v2,v3)) {}
    bool intersects(Vec3 ray, Vec3 origin){
        //check if it's in the plane first
        float d;
        float denom = normal.dot(ray);  //denominator of ray-plane intersection equation
        if(std::abs(denom) > epsilon){     //we can't calculate if denom is zero
            Vec3 pminuso = (pos-origin);         //vector from camera to intersection point
            d = pminuso.dot(normal)/denom;
        }
        if(d > epsilon){
            //check if it's in the triangle specifically
            Vec3 point = origin + d*ray;
            float alpha = triangleArea(vertex2,point,vertex3)/area;
            float beta = triangleArea(vertex1,point,vertex3)/area;
            float gamma = triangleArea(vertex1,point,vertex2)/area;
            if (alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1){
                distance = d;
                intersect = point;
                return true;
            }
        }
        return false;   //TODO: compute ray-triangle intersection
    }
    Vec3 findNormal(){
        return normal;
    }
};
*/

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
    Vec3 lightpos = Vec3(0.0f, 5.0f, -5.0f);
    float lightintensity = 1.0f;

    ///Ambient light intensity:
    float ambientlight = 0.5f;

    ///floor position & material definition:
    Vec3 floorp = Vec3(0.0f,-5.0f,0.0f);
    Colour floora = purple();
    Colour floord = red();
    Colour floors = yellow();
    float floorx = 100000.0f;            //this produces an inverse effect: higher values = less specularity
    Vec3 floorn = Vec3(0.0f,1.0f,0.0f);
    ///floor construction:
    Plane floor(floorp,floora,floord,floors,floorx,floorn);

    ///sphere1 position & material definition:
    Vec3 sphere1p = Vec3(0.0f, 0.0f, -5.0f);
    Colour sphere1a = blue();
    Colour sphere1d = green();
    Colour sphere1s = yellow();
    float sphere1x = 10.0f;
    float sphere1r = 1.0f;
    ///sphere1 construction:
    Sphere sphere1(sphere1p,sphere1a,sphere1d,sphere1s,sphere1x,sphere1r);

    /*
    ///triangle1 position & material definition:
    Vec3 tri1v1 = Vec3(2.0f, 2.0f, -1.0f);
    Vec3 tri1v2 = Vec3(2.0f, 0.0f, -1.0f);
    Vec3 tri1v3 = Vec3(0.0f, 0.0f, 0.0f);
    Colour tri1a = black();
    Colour tri1d = grey();
    Colour tri1s = white();
    float tri1x = 10.0f;
    ///tri1 construction:
    Triangle tri1(findCentroid(tri1v1,tri1v2,tri1v3),tri1a,tri1d,tri1s,tri1x,tri1v1,tri1v2,tri1v3);
    */

    ///For each row, compute rays, intersections, & shading
    for (int row = 0; row < image.rows(); ++row) {
        for (int col = 0; col < image.cols(); ++col) {

            /// (row,column) = current point on grid
            Vec3 pixel = left*U + (col*(right-left)/image.cols())*U;    //x-coordinate
            pixel += bottom*V + (row*(top-bottom)/image.rows())*V;      //y-coordinate

            ///build primary rays
            Vec3 ray = pixel - E;       //coordinate minus camera position; standard vector equation
            ray = ray.normalized();     //unit vector of ray

            ///hit detection
            bool hit = false;
            Vec3 normal,lightdir,ambient,diffuse,specular,intersection;
            float phongexp;

            if(sphere1.intersects(ray,E)){
                hit = true;
                normal = sphere1.findNormal();
                intersection = sphere1.intersect;
                lightdir = lightDirection(lightpos,intersection);
                ambient = sphere1.ambient;
                diffuse = sphere1.diffuse;
                specular = sphere1.specular;
                phongexp = sphere1.phongexp;
            }/*else if(tri1.intersects(ray,E)){
                hit = true;
                normal = tri1.findNormal();
                intersection = tri1.intersect;
                lightdir = lightDirection(lightpos,intersection);
                ambient = tri1.ambient;
                diffuse = tri1.diffuse;
                specular = tri1.specular;
                phongexp = tri1.phongexp;
            }*/else if(floor.intersects(ray,E)){
                hit = true;
                normal = floor.findNormal();
                intersection = floor.intersect;
                lightdir = lightDirection(lightpos,intersection);
                ambient = floor.ambient;
                diffuse = floor.diffuse;
                specular = floor.specular;
                phongexp = floor.phongexp;
            }

            ///compute shading (ambient only or Phong)
            if(hit){
                /// shoot ray from intersection point back towards light source
                Vec3 shadowray = lightpos - intersection;
                shadowray = shadowray.normalized();
                Vec3 epsvec = shadowray*epsilon;
                ///check shadow ray intersection
                if(sphere1.intersects(shadowray,intersection+epsvec) || floor.intersects(shadowray,intersection+epsvec)){
                    image(row,col) = ambientShading(ambient, ambientlight);
                }else{
                    image(row,col) = phongShading(ray, normal, lightdir, ambient, ambientlight, diffuse, lightintensity, specular, phongexp);
                }
            }else{
                image(row,col) = white();
            }
       }
    }

    bmpwrite("../../out.bmp", image);
    imshow(image);

    return EXIT_SUCCESS;
}
