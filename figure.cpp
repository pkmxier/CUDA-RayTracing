#include "figure.h"
#include <cmath>

Sphere::Sphere(const Material &material,
               const Vector3D &center, double R) :
    Figure(material), center(center), R(R) {
}

bool Sphere::Intersection(const Ray &ray, Vector3D &point) const {
    Vector3D fromCenter = ray.source - center;
    double b = fromCenter * ray.direction;
    double c = fromCenter * fromCenter - R * R;
    double discriminant = b * b - c;

    double t = 0;
    if (discriminant >= 0) {
        double sqrtDiscriminant = std::sqrt(discriminant);
        double t1 = -b - sqrtDiscriminant;
        double t2 = -b + sqrtDiscriminant;

        double minT = fmin(t1, t2);
        double maxT = fmax(t1, t2);

        t = (minT >= 0) ? minT : maxT;
        point = ray.source + t * ray.direction;
    }

    return t > 0;
}

Vector3D Sphere::Normal(const Vector3D &point) const {
    Vector3D normal = point - center;
    normal.Normalize();
    return normal;
}


Triangle::Triangle(const Material &material,
                   const Vector3D p[]) : Figure (material) {
    for (int i = 0; i < 3; ++i) {
        v[i] = p[i];
    }
}

Triangle::Triangle(const Material &material,
                   const Vector3D &v0, const Vector3D &v1, const Vector3D &v2) :
    Figure (material) {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
}

#include <iostream>
bool Triangle::Intersection(const Ray &ray, Vector3D &point) const {
    Vector3D v0 = v[0];
    Vector3D v1 = v[1];
    Vector3D v2 = v[2];

    Vector3D v0v1 = v1 - v0;
    Vector3D v0v2 = v2 - v0;
    // no need to normalize
    Vector3D N = v0v1 ^ v0v2; // N
    float area2 = N.Norm();

    // Step 1: finding P

    double kEpsilon = 0.01;
    // check if ray and plane are parallel ?
    float NdotRayDirection = N * ray.direction;
    if (fabs(NdotRayDirection) < kEpsilon) // almost 0
        return false; // they are parallel so they don't intersect !

    // compute d parameter using equation 2
    float d = N * v0;

    // compute t (equation 3)
    double t = (N *(ray.source) + d) / NdotRayDirection;
    // check if the triangle is in behind the ray
    if (t < 0) return false; // the triangle is behind

    // compute the intersection point using equation 1
    Vector3D P = ray.source + t * ray.direction;
    point = P;

    // Step 2: inside-outside test
    Vector3D C; // vector perpendicular to triangle's plane

    // edge 0
    Vector3D edge0 = v1 - v0;
    Vector3D vp0 = P - v0;
    C = edge0 ^(vp0);
    if (N *(C) < 0) return false; // P is on the right side

    // edge 1
    Vector3D edge1 = v2 - v1;
    Vector3D vp1 = P - v1;
    C = edge1 ^ (vp1);
    if (N *(C) < 0)  return false; // P is on the right side

    // edge 2
    Vector3D edge2 = v0 - v2;
    Vector3D vp2 = P - v2;
    C = edge2 ^ (vp2);
    if (N *(C) < 0) return false; // P is on the right side;

    return true; // this ray hits the triangle
}

Vector3D Triangle::Normal(const Vector3D &point) const {
    return ((v[2] - v[0]) ^ (v[1] - v[0])).Normalized();
}
