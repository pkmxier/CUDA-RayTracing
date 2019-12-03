#include "vector3d.h"
#include <cmath>

Vector3D::Vector3D(double x, double y, double z) : x(x), y(y), z(z) {
}

double Vector3D::Norm() const {
    return std::sqrt((*this) * (*this));
}

void Vector3D::Normalize() {
    (*this) = (*this) / this->Norm();
}

Vector3D Vector3D::Normalized() const {
    Vector3D v = *this;
    return v / v.Norm();
}

double operator *(const Vector3D &lhs, const Vector3D &rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vector3D operator *(const Vector3D &lhs, double scale) {
    return Vector3D(lhs.x * scale, lhs.y * scale, lhs.z * scale);
}

Vector3D operator *(double scale, const Vector3D &rhs) {
    return rhs * scale;
}

Vector3D Vector3D::operator /(double scale) {
    return (*this) * (1.0 / scale);
}

Vector3D operator ^(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.y * rhs.z - lhs.z * rhs.y,
                    - lhs.x * rhs.z + lhs.z * rhs.x,
                    lhs.x * rhs.y - lhs.y * rhs.x);
}

Vector3D operator +(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

Vector3D operator -(const Vector3D &lhs, const Vector3D &rhs) {
    return Vector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

