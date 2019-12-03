#include "scene.h"


bool Scene::Intersects(const Ray &ray, IntersectPoint &intersection) const {
    intersection.t = 1e8;
    IntersectPoint currentIntersection;
    for (int i = 0; i < figures.size(); ++i) {
        if (figures[i]->Intersection(ray, currentIntersection)) {
            if (currentIntersection.t < intersection.t) {
                intersection = currentIntersection;
            }
        }
    }

    return intersection.t < 1e8;
}

Vector3D Scene::ReflectVector(const Vector3D &v, const Vector3D &normal) const {
    return v - 2 * normal * (v * normal);
}

Vector3D Scene::AmbientColor(const Ray &ray, const IntersectPoint &intersection) const {
    const Material &material = intersection.material;
    return material.ambCoeff * material.color;
}


Vector3D Scene::DiffuseSpecularColor(const Ray &ray, const IntersectPoint &intersection) const {
    double diffuse = 0;
    double specular = 0;
    const Material &material = intersection.material;

    for (int i = 0; i < lights.size(); ++i) {
        Vector3D fromLight = intersection.point - lights[i].position;
        double distanceToLight = fromLight.Norm();
        fromLight.Normalize();

        double cosTheta = intersection.normal * fromLight;

        if (cosTheta <= 0) {
            Ray lightRay(lights[i].position, fromLight);
            IntersectPoint lightIntersection;

            if (Intersects(lightRay, lightIntersection)) {
                double distance = (lightIntersection.point - lights[i].position).Norm(); // = lightIntersection.t;
                if (distance < distanceToLight - 1e-8) {
                    continue;
                }
            }

            diffuse += (-cosTheta) * lights[i].intensity;

            Vector3D reflected = ReflectVector(fromLight, intersection.normal);
            double cosPhi = ray.direction * reflected;

            if (cosPhi <= 0) {
                specular += std::pow(-cosPhi, material.p) * lights[i].intensity;
            }
        }
    }

    return material.difCoeff * diffuse * material.color +
           material.specCoeff * specular * Vector3D(1, 1, 1);
}

Vector3D Scene::ReflectedColor(const Ray &ray, const IntersectPoint &intersection) const {
    return Vector3D(0, 0, 0);
}
