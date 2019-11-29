#ifndef CAMERA_H
#define CAMERA_H
#include "vector3d.h"

class CylinderPoint {
public:
    double r, phi, z;

    CylinderPoint() {}
    CylinderPoint(double r, double phi, double z) : r(r), phi(phi), z(z) {}
    Vector3D ToCartesian() {
        return Vector3D(r * sin(phi), z, r * cos(phi));
    }

    static CylinderPoint FromCartesian(const Vector3D v) {
        return CylinderPoint(sqrt(v.x * v.x + v.z * v.z),
                             atan(v.x / v.z),
                             v.y);
    }
};

class CameraParameters {
public:
    CylinderPoint amplitude;
    CylinderPoint frequency;
    CylinderPoint shift;

    CameraParameters() {}
    CameraParameters(const CylinderPoint &amplitude,
                     const CylinderPoint &frequency,
                     const CylinderPoint &shift) :
        amplitude(amplitude), frequency(frequency), shift(shift) {
    }
};

class Camera {
private:
    CylinderPoint Move(const CylinderPoint &start, const CameraParameters &params,
                   double t);
public:
    CylinderPoint position;
    CylinderPoint lookAt;
    int width;
    int height;
    double fov;
    double viewRange = 1e9;

    Camera() {}
    Camera(const CylinderPoint &position, const CylinderPoint &lookAt,
           int width, int height, double fov, double viewRange) :
        position(position), lookAt(lookAt),
        width(width), height(height), fov(fov), viewRange(viewRange) {
    }

    void Move(const CylinderPoint &startPos, const CylinderPoint &startLookAt,
              double t,
              const CameraParameters &posParams,
              const CameraParameters &lookAtParams);
};

#endif // CAMERA_H
