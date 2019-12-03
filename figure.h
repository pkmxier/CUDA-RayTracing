#ifndef FIGURE_H
#define FIGURE_H
#include "color.h"
#include "vector3d.h"
#include <vector>
#include <iostream>

using namespace std;

class Ray {
public:
    Vector3D source;
    Vector3D direction;

    Ray() {}
    Ray(const Vector3D &s, const Vector3D &d) : source(s), direction(d) {
    }
};

class Material {
public:
    double ambCoeff, difCoeff, specCoeff;

    Vector3D color;

    double p;
    double refraction;

    Material() {}
    Material(const Vector3D &color, double ambCoeff, double difCoeff,
             double specCoeff, double p, double refraction) :
        refraction(refraction), p(p), color(color), ambCoeff(ambCoeff),
        difCoeff(difCoeff), specCoeff(specCoeff) {
    }
};

struct IntersectPoint {
    Vector3D point;
    Vector3D normal;
    double t;

    Material material;
};

class Figure {
public:
    Material material;

    Figure(const Material &material) : material(material) {
    }

    virtual bool Intersection(const Ray &ray, IntersectPoint &intersection) const = 0;
    virtual Vector3D Normal(const Vector3D &point) const = 0;
};

class Sphere : public Figure {
public:
    Vector3D center;
    double R;

    Sphere(const Material &material,
           const Vector3D &c, double R);

    bool Intersection(const Ray &ray, IntersectPoint &intersection) const;
    Vector3D Normal(const Vector3D &point) const;
};

class Triangle : public Figure {
public:
    Vector3D v[3];

    Triangle(const Material &material, const Vector3D p[3]);
    Triangle(const Material &material,
             const Vector3D &v0, const Vector3D &v1, const Vector3D &v2);

    bool Intersection(const Ray &ray, IntersectPoint &intersection) const;
    Vector3D Normal(const Vector3D &point) const;
};

#endif // FIGURE_H
