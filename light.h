#ifndef LIGHT_H
#define LIGHT_H
#include "vector3d.h"
#include "color.h"

class Light {
public:
    Vector3D position;
    double intensity;

    Light(const Vector3D &p, double i) : position(p), intensity(i){
    }
};

#endif // LIGHT_H
