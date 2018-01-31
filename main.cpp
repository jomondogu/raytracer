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

///Returns true if the given ray intersects with a sphere at point ems from camera with radius r, false otherwise
bool raySphereIntersect(float disc){
    if(disc >= 0){  //then ray intersects sphere
        return true;
    }
    return false;
}

///Returns true if the given ray from intersects with a plane with point p and normal n, false otherwise
bool rayPlaneIntersect(float t){
    if(t >= 0){
        return true;
    }
    return false;
}

///Computes the ambient light for a pixel based on ambient coefficient of a surface & ambient light intensity
Colour ambientShading(Colour ka, float ia){
    return ka*ia;
}

///Computes the colour of a pixel based on the given ray, normal to surface, vector pointing to light source, and surface/light coefficients/Phong exponent
Colour phongShading(Vec3 v, Vec3 n, Vec3 l, Colour ka, float ia, Colour kd, float i, Colour ks, float p){
    ///normalize all vectors
    v = -v;
    v = v.normalized();
    n = n.normalized();
    l = l.normalized();

    Colour ambient = ambientShading(ka,ia);

    Colour diffuse = kd*i*std::fmaxf(l.dot(n),0.0f);

    Vec3 h = (v+l).normalized();
    Colour specular = ks*i*std::powf(std::fmaxf(h.dot(n),0.0f),p);

    return ambient + diffuse + specular;
}

//Vec3 point = (left->right(x), up->down(y), far->near(z))
//Vec3 vector = (right(x),up(x),near(x))
int main(int, char**){

    int wResolution = 1080;
    int hResolution = 720;
    float aratio = float(wResolution)/float(hResolution);
    // #rows = hResolution, #cols = wResolution
    Image<Colour> image(hResolution, wResolution);

    ///Camera position definition:
    Vec3 W = Vec3(0.0f,0.0f,-1.0f); //out direction
    Vec3 V = Vec3(0.0f,1.0f,0.0f);  //up direction
    Vec3 U = Vec3(1.0f,0.0f,0.0f);  //side direction
    float d = 1.0f;                 //focal length
    Vec3 E = -d*W;

    ///Ambient light intensity:
    float ambientLight = 0.2f;

    ///Sphere position & material definition:
    Vec3 spherePos = Vec3(3.0f, 1.0f, -5.0f);   //the intersection equation defines the sphere's boundaries later; we only need pos & radius here
    float sphereRadius = 3.0f;
    Colour sphereAmbient = Colour(0.5f,0.5f,0.5f);
    Colour sphereDiffuse = red();
    Colour sphereSpecular = white();
    float spherePhong = 10.0f;

    ///plane position definition:
    Vec3 floorPos = Vec3(0.0f,-50.0f,0.0f);  //positioned at "floor height"
    Vec3 floorNormal = Vec3(0.0f,-1.0f,0.0f);  //directed straight up
    Colour floorAmbient = Colour(0.5f,0.5f,0.5f);
    Colour floorDiffuse = Colour(0.0f,0.2f,0.0f);
    Colour floorSpecular = white();
    float floorPhong = 10.0f;

    ///view bounds:
    float left = -1.0f*aratio;     //leftmost view plane boundary
    float right = 1.0f*aratio;     //rightmost
    float top = 1.0f;       //top
    float bottom = -1.0f;   //bottom
    ///without aratio, image will be distorted: our bounds define a square but the image aspect ratio is 4:3

    ///light source & intensity:
    Vec3 lightPos = Vec3(0.0f, 0.0f, 0.0f);
    float lightIntensity = 1.0f;

    for (int row = 0; row < image.rows(); ++row) {
        for (int col = 0; col < image.cols(); ++col) {

            /// (row,column) = current point on grid
            Vec3 pixel = left*U + (col*(right-left)/image.cols())*U;    //x-coordinate
            pixel += bottom*V + (row*(top-bottom)/image.rows())*V;      //y-coordinate

            ///build primary rays
            Vec3 ray = pixel - E;       //coordinate minus camera position; standard vector equation
            ray = ray.normalized();     //unit vector of ray
            
            /// ray-sphere intersection
            Vec3 ems = E - spherePos; //camera position minus sphere position
            float disc = std::powf(ray.dot(ems),2) - ems.dot(ems) + sphereRadius*sphereRadius; //discriminant: https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

            ///ray-plane intersection: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
            floorNormal = floorNormal.normalized();
            float denom = ray.dot(floorNormal);              //denominator of ray-plane intersection equation
            float t = -1;
            if(denom > 0.00000000001f){                    //we can't calculate if denom is zero
                Vec3 omp = E - floorPos;               //plane point minus origin pixel
                t = -omp.dot(floorNormal);
            }

            if(raySphereIntersect(disc)){
                float d = -(ray.dot(ems)) - std::sqrtf(disc);   //find distance where ray hits sphere
                Vec3 pos = E + d*ray;                           //find actual position
                Vec3 sphereNormal = (pos - spherePos)/sphereRadius;   //find normal to sphere at that position
                Vec3 lightDir = lightPos - pos;
                image(row,col) = phongShading(ray, sphereNormal, lightDir, sphereAmbient, ambientLight, sphereDiffuse, lightIntensity, sphereSpecular, spherePhong);
            }else if(rayPlaneIntersect(t)){
                //find distance, position, & normal to intersection point
                Vec3 pos = E + t*ray;
                Vec3 lightDir = pos - lightPos;
                image(row,col) = phongShading(ray, floorNormal, lightDir, floorAmbient, ambientLight, floorDiffuse, lightIntensity, floorSpecular, floorPhong);
            }else{
                image(row,col) = white();       //and here
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
