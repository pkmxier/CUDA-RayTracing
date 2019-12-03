#include <iostream>
#include "camera.h"
#include "figure.h"
#include "light.h"
#include "vector3d.h"
#include "scene.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include "omp.h"
#include "time.h"
#include <cstdlib>

#define M_PI 3.14159265359

using namespace std;

Vector3D CastRay(const Ray &ray, const Scene &scene) {
    IntersectPoint intersection;
    if (scene.Intersects(ray, intersection)) {
        Vector3D ambient = scene.AmbientColor(ray, intersection);
        Vector3D diffuseSpecular = scene.DiffuseSpecularColor(ray, intersection);
        Vector3D reflected = scene.ReflectedColor(ray, intersection);

        Vector3D color = ambient + diffuseSpecular;

        return color;// + reflected;
    } else {                                                          // default color
        Vector3D unit_direction = ray.direction.Normalized();
        float t = 0.5 * (unit_direction.y + 1.0);
        return (1.0 - t) * Vector3D(1, 1, 1) + t * Vector3D(1, 0, 0);
    }
}

float clamp(const float &lo, const float &hi, const float &v) {
    return std::max(lo, std::min(hi, v));
}

void ImageToFile(vector<Vector3D> &image, int w, int h, const string &fileName) {
    ofstream ofs;
    ofs.open(fileName, ios::binary);
    ofs << "P6\n" << w << " " << h << "\n255\n";

    for (int i = 0; i < image.size(); ++i) {
        Vector3D &c = image[i];

        char r = (char)(255 * clamp(0, 1, c.x));
                char g = (char)(255 * clamp(0, 1, c.y));
                char b = (char)(255 * clamp(0, 1, c.z));
        ofs << r << g << b;
//        ofs << char(255.99 * c.x) << char(255.99 * c.y) << char(255.99 * c.z);
    }

    ofs.close();
}

vector<Vector3D> Antialiasing(vector<Vector3D> &image, int n, int width) {
    if (n == 1) {
        return image;
    }

    int cnt = 0;
    vector<Vector3D> result(image.size() / (n * n), Vector3D(0, 0, 0));
    for (int k = 0; k < image.size() && cnt < result.size(); k += n) {
        int x = k % width;
        int y = k / width;
        if (y % n != 0) {
            k += width * (n - 1);
            y = k / width;
        }


        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int index = x + i + (y + j) * width;
                result[cnt] = result[cnt] + image[index];
            }
        }

        result[cnt] = result[cnt] / (n * n);
        ++cnt;
    }

    return result;
}

//#define ITERATIVE
#define NON_ITERATIVE

void Render(Scene &scene, int frames, int aa = 1) {
    Camera &camera = scene.camera;
    camera.height *= aa;
    camera.width *= aa;
    double scale = 0.5 / tan(camera.fov * M_PI / (2 * 180.0));
    double ratio = camera.height / double(camera.width);
    Vector3D up(0, 1, 0);

    for (double t = 0; t < 2 * M_PI; t += (2 * M_PI) / frames) {
        #ifdef NON_ITERATIVE
            vector<Vector3D> image(camera.width * camera.height);
        #endif

        #ifdef ITERATIVE
            ofstream ofs;
            ofs.open(string("../images/out") + to_string(t * 180 / M_PI) + string(".ppm"), ios::binary);
            ofs << "P6\n" << camera.width << " " << camera.height << "\n255\n";
        #endif

        cout << "t = " << t * 180 / M_PI << ": ";
        CameraPoints cameraPoints = camera.Move(t);

        Vector3D toCenter = (cameraPoints.lookAt.ToCartesian() -
                             cameraPoints.position.ToCartesian()).Normalized();
        Vector3D u = up ^ toCenter;
        Vector3D v = toCenter ^ u;
        Ray ray;

        for (int j = 0; j < camera.height; ++j) {
            for (int i = 0; i < camera.width; ++i) {
                double x = double(i) / camera.width  - 0.5;
                double y = -(double(j) / camera.height - 0.5);

                ray.source = cameraPoints.position.ToCartesian();
                ray.direction = (toCenter * scale + x * u + y * v * ratio).Normalized();


                #ifdef NON_ITERATIVE
                    image[i + j * camera.width] = CastRay(ray, scene);
                #endif

                #ifdef ITERATIVE
                    Vector3D c = CastRay(ray, scene);
                    ofs << char(255.99 * c.x) << char(255.99 * c.y) << char(255.99 * c.z);
                #endif
            }
        }
        cout << "done\n";

        #ifdef NON_ITERATIVE
            image = Antialiasing(image, aa, camera.width);
            ImageToFile(image, camera.width / aa, camera.height / aa,
                        string("../images/out") + to_string(t * 180 / M_PI) + string(".ppm"));
        #endif

        #ifdef ITERATIVE
            ofs.close();
        #endif
    }
}

int main() {
    int frames = 1;
    Scene scene;
    CameraParameters positionParams(CylinderPoint(0, 0, 0),
                                    CylinderPoint(0, 1, 0),
                                    CylinderPoint(0, 0, 0));
    CameraParameters lookAtParams(CylinderPoint(0, 0, 0),
                                  CylinderPoint(0, 0, 0),
                                  CylinderPoint(0, 0, 0));
    scene.camera = Camera(CylinderPoint(1, M_PI, 0), CylinderPoint(0, 0, 0), // position, lookAtPosition
                          positionParams, lookAtParams,                   // movement parameters
                          800, 600, 90);                                  // width, height, fov


    Material green(Vector3D(0, 1, 0), 0.2, 0.4, 0.4, 10, 0);
    Material blue(Vector3D(0, 0, 1), 0.2, 0.6, 0.2, 3, 0);
    scene.lights.push_back(Light(Vector3D(1, 5, 0), 1.5));
    scene.figures.push_back(new Sphere(green, Vector3D(0, 0, 1), 0.5));
    scene.figures.push_back(new Sphere(blue, Vector3D(0, -100.5, 1), 100));

    cout << "Rendering started\n";
    Render(scene, frames, 1);
    cout << "Done\n";

    return 0;
}

//void PushS(const Material &material, std::vector<Figure *> &spheres,
//           Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
//    spheres.push_back(new Triangle(material, a, c, b));
//    spheres.push_back(new Triangle(material, a, d, c));
//}

//void PushSquare(const Material &material, std::vector<Figure *> &spheres,
//                const Vector3D &center, double size, std::vector<Light> &lights) {
//    vector<Vector3D> v(8);
//    v[0] = center - size * Vector3D(1, 1, 1);
//    v[1] = center - size * Vector3D(1, -1, 1);
//    v[2] = center - size * Vector3D(-1, -1, 1);
//    v[3] = center - size * Vector3D(-1, 1, 1);

//    v[4] = center - size * Vector3D(1, 1, -1);
//    v[5] = center - size * Vector3D(1, -1, -1);
//    v[6] = center - size * Vector3D(-1, -1, -1);
//    v[7] = center - size * Vector3D(-1, 1, -1);

//    PushS(red, spheres, v[0], v[1], v[2], v[3]);
//    PushS(green, spheres, v[0], v[4], v[5], v[1]);
//    PushS(yellow, spheres, v[1], v[5], v[6], v[2]);
//    PushS(blue, spheres, v[2], v[6], v[7], v[3]);
//    PushS(white, spheres, v[0], v[3], v[7], v[4]);
//    PushS(grey, spheres, v[4], v[7], v[6], v[5]);
//}
