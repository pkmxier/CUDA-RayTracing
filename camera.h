#ifndef CAMERA_H
#define CAMERA_H
#include "vector3d.h"

class CylinderPoint {
public:
    double r, phi, z;

    CylinderPoint() {}
    CylinderPoint(double r, double phi, double z) : r(r), phi(phi), z(z) {}
    Vector3D ToCartesian() const {
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

struct CameraPoints {
    CylinderPoint position;
    CylinderPoint lookAt;

    CameraPoints() {}
    CameraPoints(const CylinderPoint &pos, const CylinderPoint &lookAt) : position(pos), lookAt(lookAt) {}
};

class Camera {
private:
    CylinderPoint Move(const CylinderPoint &start, const CameraParameters &params, double t) const;
public:
    CameraPoints points;
    CameraParameters posParams;
    CameraParameters lookAtParams;
    int width;
    int height;
    double fov;

    Camera() {}
    Camera(const CylinderPoint &position, const CylinderPoint &lookAt,
           const CameraParameters &posParams,
           const CameraParameters &lookAtParams,
           int width, int height, double fov) :
        points(position, lookAt), posParams(posParams), lookAtParams(lookAtParams),
        width(width), height(height), fov(fov) {
    }

    CameraPoints Move(double t) const;
};

#endif // CAMERA_H
