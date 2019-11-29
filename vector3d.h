#ifndef VECTOR3D_H
#define VECTOR3D_H
#include <iostream>

class Vector3D {
public:
    double x, y, z;

    Vector3D() {}
    Vector3D(double x, double y, double z);
    Vector3D(const Vector3D &a, const Vector3D &b);

    double Norm() const;
    double Length();
    void Normalize();
    Vector3D Normalized() const;

    friend double operator *(const Vector3D &lhs, const Vector3D &rhs);
    Vector3D operator /(double scale);
    Vector3D operator ^(const Vector3D &rhs);
    friend Vector3D operator +(const Vector3D &lhs, const Vector3D &rhs);
    friend Vector3D operator -(const Vector3D &lhs, const Vector3D &rhs);

    friend Vector3D operator *(const Vector3D &lhs, double scale);
    friend Vector3D operator *(double scale, const Vector3D &rhs);

    friend std::ostream & operator <<(std::ostream &os, const Vector3D &v) {
        os << v.x << " " << v.y << " " << v.z << std::endl;
        return os;
    }
};



#endif // VECTOR3D_H
