#ifndef SCENE_H
#define SCENE_H
#include "figure.h"
#include "light.h"
#include "camera.h"
#include <vector>

class Scene {
public:
    std::vector<Figure *> figures;
    std::vector<Light> lights;
    Camera camera;

    Scene() {}

    bool Intersects(const Ray &ray, IntersectPoint &intersection) const;

    Vector3D ReflectVector(const Vector3D &v, const Vector3D &normal) const;

    Vector3D AmbientColor(const Ray &ray, const IntersectPoint &intersection) const;
    Vector3D DiffuseSpecularColor(const Ray &ray, const IntersectPoint &intersection) const;
    Vector3D ReflectedColor(const Ray &ray, const IntersectPoint &intersection) const;
};

#endif // SCENE_H
