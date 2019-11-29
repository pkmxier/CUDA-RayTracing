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
#include <string>
#include "omp.h"

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

    return spheresDistance < 1000;
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
        return Vector3D(0.3, 0.4, 0.8);
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

void PushS(const Material &material, std::vector<Figure *> &spheres,
           Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
    spheres.push_back(new Triangle(material, a, b, c));
    spheres.push_back(new Triangle(material, a, d, c));
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

        lights.push_back(Light(0.5*(v[0] + v[2]), 2));
}

void Render(const vector<Figure *> &figures, const vector<Light> &lights,
            Camera &camera, int frames, const string &fileName) {
    int w = camera.width;
    int h = camera.height;
    double fov = camera.fov;

    camera.position = CylinderPoint(10, M_PI, 0);
    camera.lookAt = CylinderPoint(0, 0, 0);

    vector<Vector3D> image(w * h);
    CameraParameters params1;
    params1.amplitude = CylinderPoint(0, 0, 0);
    params1.frequency = CylinderPoint(M_PI / 6, 1, M_PI / 6);
    params1.shift = CylinderPoint(0, 0, 0);

    CameraParameters params2;
    params2.amplitude = CylinderPoint(0, 0, 0);
    params2.frequency = CylinderPoint(0, 0, 0);
    params2.shift = CylinderPoint(0, 0, 0);

    CylinderPoint pos = camera.position;
    CylinderPoint lookAt = camera.lookAt;

    for (double t = 0; t < 2 * M_PI; t += M_PI / 6.0) {
        camera.Move(pos, lookAt, t, params1, params2);

        Vector3D toCenter = (camera.lookAt.ToCartesian() - camera.position.ToCartesian()).Normalized();
        Vector3D up(0, 1, 0);
        Vector3D u = up ^ toCenter;
        Vector3D v = toCenter ^ u;
        Vector3D direction;
        double scale = 0.5 / tan(fov / 2);
        cout << scale <<endl;
        double ratio = h / double(w);
        Ray ray;

        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                double x = double(i) / w  - 0.5;
                double y = -double(j) / h + 0.5;
                direction = (toCenter * scale + x * u + y * v * ratio).Normalized();

                ray.source = camera.position.ToCartesian();
                ray.direction = direction.Normalized();
                image[i + j * w] = CastRay(ray, figures, lights);
            }
        }
        cout << "t = " << t * 180 / M_PI << " done\n";

        ofstream ofs;
        ofs.open(string("../out") + to_string(t * 180 / M_PI) + string(".ppm"), std::ios::binary);
        ofs << "P6\n" << w << " " << h << "\n255\n";

        for (int i = 0; i < w * h; ++i) {
            Vector3D &c = image[i];

            double m = max(c.x, max(c.y, c.z));
            if (m > 1) {
                c = c / m;
            }

            ofs << (char)(255 * max(0.0, min(1.0, image[i].x)));
            ofs << (char)(255 * max(0.0, min(1.0, image[i].y)));
            ofs << (char)(255 * max(0.0, min(1.0, image[i].z)));
        }

        ofs.close();
    }
}

using v4 = vector<double>;
int main() {
    Material ivory(1.0, v4{0.6,  0.3, 0.1, 0.0}, Vector3D(0.4, 0.4, 0.3),   50.);
    Material glass(1, v4{0.0,  0.5, 0.1, 0.8}, Vector3D(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, v4{0.9,  0.1, 0.0, 0.0}, Vector3D(0.3, 0.1, 0.1),   10.);
    Material mirror(1.0, v4{0.0, 10.0, 0.8, .0}, Vector3D(1.0, 1.0, 1.0), 1425.);

    std::vector<Figure *> figures;
    std::vector<Light> lights;
    lights.push_back(Light(Vector3D(0, 1, 10), 1.5));
//    lights.push_back(Light(Vector3D(0, 1, -10), 1.5));

    figures.push_back(new Triangle(red_rubber, Vector3D(-1, 0, 0),
                                   Vector3D(0, 1, 0),
                                   Vector3D(1, 0, 0)));

    Camera camera;
    camera.fov = M_PI / 2.0;
    camera.width = 800;
    camera.height = 600;

    cout << "Rendering started\n";
    Render(figures, lights, camera, 1, "");
    cout << "Done\n";

    return 0;
}
