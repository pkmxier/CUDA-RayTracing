#include "camera.h"
#include <cmath>

CylinderPoint Camera::Move(const CylinderPoint &start, const CameraParameters &params,
                      double t) {
    CylinderPoint newPoint = start;
    newPoint.r += params.amplitude.r * sin(params.frequency.r * t + params.shift.r);
    newPoint.z += params.amplitude.z * sin(params.frequency.z * t + params.shift.z);
    newPoint.phi += params.frequency.phi * t;

    return newPoint;
}

void Camera::Move(const CylinderPoint &startPosition, const CylinderPoint &startLookAt,
                  double t,
                  const CameraParameters &positionParams,
                  const CameraParameters &lookAtParams) {
    position = Move(startPosition, positionParams, t);
    lookAt = Move(startLookAt, lookAtParams, t);
}
