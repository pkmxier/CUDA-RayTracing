#ifndef LIGHT_H
#define LIGHT_H
#include "vector3d.h"
#include "color.h"

class Light {
public:
    Vector3D position;
    Color color;
    double intensity;

    Light(const Vector3D &pos, double intensity);
};

#endif // LIGHT_H
