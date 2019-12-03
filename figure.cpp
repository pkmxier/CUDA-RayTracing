#include "figure.h"
#include <cmath>

Sphere::Sphere(const Material &material,
               const Vector3D &center, double R) :
    Figure(material), center(center), R(R) {
}

bool Sphere::Intersection(const Ray &ray, IntersectPoint &intersection) const {
    Vector3D fromCenter = ray.source - center;
    double b = fromCenter * ray.direction;
    double c = fromCenter * fromCenter - R * R;
    double discriminant = b * b - c;

    intersection.t = 0;
    if (discriminant >= 0) {
        double sqrtDiscriminant = std::sqrt(discriminant);
        double t1 = -b - sqrtDiscriminant;
        double t2 = -b + sqrtDiscriminant;

        double minT = fmin(t1, t2);
        double maxT = fmax(t1, t2);

        intersection.t = (minT >= 0) ? minT : maxT;
        intersection.point = ray.source + intersection.t * ray.direction;
    }

    intersection.normal = Normal(intersection.point);
    intersection.material = material;

    return intersection.t > 0;
}

Vector3D Sphere::Normal(const Vector3D &point) const {
    Vector3D normal = (point - center).Normalized();
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

bool Triangle::Intersection(const Ray &ray, IntersectPoint &intersection) const {
    Vector3D v0v1 = v[1] - v[0];
    Vector3D v0v2 = v[2] - v[0];
    Vector3D pvec = ray.direction ^ v0v2;
    double det = v0v1 * pvec;

    double BIAS = 1e-8;
    if (det < BIAS) return false;
    if (std::fabs(det) < BIAS) return false;

    double invDet = 1 / det;

    Vector3D tvec = ray.source - v[0];
    double u = tvec * pvec * invDet;
    if (u < 0 || u > 1) return false;

    Vector3D qvec = tvec ^ v0v1;
    double v = ray.direction * qvec * invDet;
    if (v < 0 || u + v > 1) return false;

    intersection.t = v0v2 * qvec * invDet;
    intersection.normal = Normal(v0v1);
    intersection.point = ray.source + intersection.t * ray.direction;
    intersection.material = material;

    return intersection.t > BIAS;
}

Vector3D Triangle::Normal(const Vector3D &point) const {
    return ((v[1] - v[0]) ^ (v[2] - v[0])).Normalized();
}
