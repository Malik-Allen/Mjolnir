// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_TYPE_DEFS_H
#define MJOLNIR_TYPE_DEFS_H

#include "../src/Quaternion.h"
#include "../src/Vector2.h"

typedef Mjolnir::Math::Matrix3<float> Matrix3f;
typedef Mjolnir::Math::Matrix3<double> Matrix3d;

typedef Mjolnir::Math::Matrix4<float> Matrix4f;
typedef Mjolnir::Math::Matrix4<double> Matrix4d;

typedef Mjolnir::Math::Vector2<int> Vector2i;
typedef Mjolnir::Math::Vector2<float> Vector2f;
typedef Mjolnir::Math::Vector2<double> Vector2d;

typedef Mjolnir::Math::Vector3<int> Vector3i;
typedef Mjolnir::Math::Vector3<float> Vector3f;
typedef Mjolnir::Math::Vector3<double> Vector3d;

typedef Mjolnir::Math::Vector4<int> Vector4i;
typedef Mjolnir::Math::Vector4<float> Vector4f;
typedef Mjolnir::Math::Vector4<double> Vector4d;

typedef Mjolnir::Math::Quaternion<float> Quatf;
typedef Mjolnir::Math::Quaternion<double> Quatd;

#endif //MJOLNIR_TYPE_DEFS_H