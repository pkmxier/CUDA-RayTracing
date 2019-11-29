#include <iostream>
#include "camera.h"
#include "figure.h"
#include "light.h"
#include "vector3d.h"
#include "color.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>

#define M_PI 3.14159265359

using namespace std;

bool SceneIntersect(const Ray &ray, const vector<Figure *> &spheres,
                     Vector3D &hit, Vector3D &N, Material &material) {
    double spheresDistance = std::numeric_limits<double>::max();
    for (int i = 0; i < spheres.size(); ++i) {
        Vector3D point;

        if (spheres[i]->Intersection(ray, point)) {
            double distance = (ray.source - point).Norm();
            if (distance < spheresDistance) {
                spheresDistance = distance;
                hit = point;
                Vector3D normal = spheres[i]->Normal(hit);
                N = normal.Normalized();

                if (N * ray.direction > 0) {
                    N = N * (-1);
                }

                material = spheres[i]->material;
            }
        }
    }

    double checkerboardDist = std::numeric_limits<double>::max();

//    if (fabs(ray.direction.y) > 1e-3) {
//        double d = -(ray.source.y + 4) / ray.direction.y;
//        Vector3D pt = ray.source + ray.direction * d;
//        if (d > 0 && fabs(pt.x) < 10 &&
//                pt.z < -10 && pt.z > -30 && d < spheresDistance) {
//            checkerboardDist = d;
//            hit = pt;
//            N = Vector3D(0, 1, 0);
//            material.albedo = vector<double>(4, 0);
//            material.albedo[0] = 1;
//            material.color =
//                    (int(.5 * hit.x + 1000) + int(.5 * hit.z)) &
//                    1 ? Vector3D(.3, .3, .3) : Vector3D(.3, .2, .1);
//        }
//    }

    return min(spheresDistance, checkerboardDist) < 1000;
}

Ray Reflect(const Ray &I, const Vector3D &N, const Vector3D &point) {
    Ray ray;
    ray.direction = I.direction - N * 2.0 * (I.direction * N);
    ray.source = ray.direction * N < 0 ? point - N * (1e-3): point + N * 1e-3;
    return ray;
}

Ray Refract(const Ray &I, const Vector3D &N, const Vector3D &point,
            const double eta_t, const double eta_i = 1.0) {
    double cosi = - max(-1.0, min(1.0, I.direction * N));
    if (cosi < 0) {
        return Refract(I, N * (-1), point, eta_i, eta_t);
    }
    double eta = eta_i / eta_t;
    double k = 1 - eta * eta * (1 - cosi * cosi);

    Ray ray;
    ray.direction = k < 0 ?
             Vector3D(1, 0, 0) : I.direction * eta + N * (eta * cosi - sqrt(k));
    ray.source = ray.direction * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    return ray;
}

Vector3D CastRay(const Ray &ray,
                 const vector<Figure *> &spheres, const vector<Light> &lights,
                 int depth = 0) {
    Vector3D point, N;
    Material material;

    if (depth > 4 || !SceneIntersect(ray, spheres, point, N, material)) {
        return Vector3D(0.4, 0.3, 1);
    }

    Ray reflectRay = Reflect(ray, N, point);
    Ray refractRay = Refract(ray, N, point, material.refraction);

    Vector3D reflectColor = CastRay(reflectRay, spheres, lights, depth + 1);
    Vector3D refractColor = CastRay(refractRay, spheres, lights, depth + 1);

    double diffuseIntensity = 0, specularIntensity = 0;
    for (int i = 0; i < lights.size(); ++i) {
        Ray lightRay;
        lightRay.direction = lights[i].position - point;
        double lightDistance = lightRay.direction.Norm();
        lightRay.direction.Normalize();

        lightRay.source = (lightRay.direction * N < 0) ?
                    point - N * 1e-3 : point + N * 1e-3;
        Vector3D shadowPt, shadowN;
        Material tmpMaterial;

        if (SceneIntersect(lightRay, spheres, shadowPt, shadowN, tmpMaterial)) {
            if (shadowPt.Norm() < lightDistance) {
                continue;
            }
        }

        diffuseIntensity +=
                lights[i].intensity * max(0.0, lightRay.direction * N);
        specularIntensity +=
                powf(max(0.0,(-1) *
                         Reflect(lightRay, N, point).direction * ray.direction),
                     material.power * lights[i].intensity);
    }

    Vector3D color = material.color * diffuseIntensity * material.albedo[0] +
            Vector3D(1, 1, 1) * specularIntensity * material.albedo[1] +
            reflectColor * material.albedo[2] +
            refractColor * material.albedo[3];

    return color;
}

void Render(const vector<Figure *> &spheres, const vector<Light> &lights) {
    int width = 800;
    int height = 600;
    double fov = M_PI / 3;

    vector<Vector3D> out(width * height);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            Vector3D start(0, 0, 0);
            Vector3D direction;
            direction.x = (i + 0.5) - width / 2.0;
            direction.y = -(j + 0.5) + height / 2.0;
            direction.z = -height / (2.0 * tan(fov / 2.0));
            direction.Normalize();

            Ray ray(start, direction);
            out[i + j * width] = CastRay(ray, spheres, lights);
        }
    }

    ofstream ofs;
    ofs.open("../out.ppm", std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";

    for (int i = 0; i < height * width; ++i) {
        Vector3D &c = out[i];

        double m = max(c.x, max(c.y, c.z));
        if (m > 1) {
            c = c / m;
        }

        ofs << (char)(255 * max(0.0, min(1.0, out[i].x)));
        ofs << (char)(255 * max(0.0, std::min(1.0, out[i].y)));
        ofs << (char)(255 * max(0.0, std::min(1.0, out[i].z)));
    }

    ofs.close();
}

void PushS(const Material &material, std::vector<Figure *> &spheres,
           Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
    spheres.push_back(new Triangle(material, a, b, c));
    spheres.push_back(new Triangle(material, a, d, c));
//    cout << a << b << c << d << endl;
}

void PushSquare(const Material &material, std::vector<Figure *> &spheres,
                const Vector3D &center, double size, std::vector<Light> &lights) {
    vector<Vector3D> v(8);
    v[0] = center - size * Vector3D(1, 1, 1);
    v[1] = center - size * Vector3D(1, -1, 1);
    v[2] = center - size * Vector3D(-1, -1, 1);
    v[3] = center - size * Vector3D(-1, 1, 1);

    v[4] = center - size * Vector3D(1, 1, -1);
    v[5] = center - size * Vector3D(1, -1, -1);
    v[6] = center - size * Vector3D(-1, -1, -1);
    v[7] = center - size * Vector3D(-1, 1, -1);

    PushS(material, spheres, v[0], v[1], v[2], v[3]);
    PushS(material, spheres, v[0], v[1], v[5], v[4]);
    PushS(material, spheres, v[1], v[5], v[6], v[2]);
    PushS(material, spheres, v[2], v[6], v[7], v[3]);
    PushS(material, spheres, v[0], v[4], v[7], v[3]);
    PushS(material, spheres, v[4], v[5], v[6], v[7]);

    for (int i = 0; i < v.size(); ++i) {
        lights.push_back(Light(v[i] + ((center - v[i]).Normalized() * 0.1), 2));
    }
}

using v4 = vector<double>;
int main() {
    Material ivory(1.0, v4{0.6,  0.3, 0.1, 0.0}, Vector3D(0.4, 0.4, 0.3),   50.);
    Material glass(1, v4{0.0,  0.5, 0.1, 0.8}, Vector3D(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, v4{0.9,  0.1, 0.0, 0.0}, Vector3D(0.3, 0.1, 0.1),   10.);
    Material mirror(1.0, v4{0.0, 10.0, 0.8, .0}, Vector3D(1.0, 1.0, 1.0), 1425.);

    std::vector<Figure *> spheres;
//    spheres.push_back(new Sphere(ivory, Vector3D(-3,    0,   -16), 2));
//    spheres.push_back(new Sphere(glass, Vector3D(0, 0, -10), 3));
//    spheres.push_back(new Sphere(red_rubber, Vector3D(20, 0, 5), 2));
//    spheres.push_back(new Sphere(red_rubber, Vector3D( 1.5, -0.5, -18), 3));
//    spheres.push_back(new Sphere(mirror, Vector3D( 7,    5,   -18), 4));

//    spheres.push_back(new Triangle(glass,
//                                   Vector3D(-5, -2, -10),
//                                   Vector3D(-5, 5, -10),
//                                   Vector3D(5, 5, -10)));
//    spheres.push_back(new Triangle(glass,
//                                   Vector3D(-5, -2, -10),
//                                   Vector3D(5, -2, -10),
//                                   Vector3D(5, 5, -10)));

//    spheres.push_back(new Triangle(glass,
//                                   Vector3D(-5, -2, 2),
//                                   Vector3D(5, -2, 2),
//                                   Vector3D(0, 5, 2)));

    std::vector<Light> lights;
    lights.push_back(Light(Vector3D(0, 3, 3), 1.5));
//    lights.push_back(Light(Vector3D( 30, 50, -25), 1.8));
//    lights.push_back(Light(Vector3D( 30, 20,  30), 1.7));
    PushSquare(glass, spheres, Vector3D(3, -3, -12), 2, lights);


    cout << "Rendering started\n";
    Render(spheres, lights);
    cout << "Done\n";

    return 0;
}
