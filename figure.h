#ifndef FIGURE_H
#define FIGURE_H
#include "color.h"
#include "vector3d.h"
#include <vector>

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
    double ambient, diffuse, specular;

    Vector3D color;

    double transparency;
    double power;
    double refraction;

    std::vector<double> albedo;

    Material() {}
    Material(double refraction,
             const std::vector<double> &albedo,
             const Vector3D &color,
             double power) :
        refraction(refraction), albedo(albedo), power(power), color(color) {
    }
};

class Figure {
public:
    Material material;

    Figure(const Material &material) : material(material) {
    }

    virtual bool Intersection(const Ray &ray, Vector3D &point) const = 0;
    virtual Vector3D Normal(const Vector3D &point) const = 0;
};

class Sphere : public Figure {
public:
    Vector3D center;
    double R;

    Sphere(const Material &material,
           const Vector3D &c, double R);

    bool Intersection(const Ray &ray, Vector3D &point) const;
    Vector3D Normal(const Vector3D &point) const;
};

class Triangle : public Figure {
public:
    Vector3D v[3];

    Triangle(const Material &material, const Vector3D p[3]);
    Triangle(const Material &material,
             const Vector3D &v0, const Vector3D &v1, const Vector3D &v2);

    bool Intersection(const Ray &ray, Vector3D &point) const;
    Vector3D Normal(const Vector3D &point) const;
};

#endif // FIGURE_H
