// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_MATRIX4_H
#define MJOLNIR_MATRIX4_H

#include "Vector3.h"
#include "Vector4.h"

	// Shoutout and thank you to Scott
	/// These are the default vectors of the eye (camera) according to OpenGL and the right hand rule
	///								+Y   -Z
	///	                             |  /
	///   origin(0.0f,0.0f,0.0f);    | /
	///   up(0.0f,1.0f,0.0f);        |/_____+X
	///   forward(0.0f,0.0f,-1.0f);  /
	///                             /
	///                           +Z
	///

	/// Let's just make sure that all is clear about how this matrix is layed out. 

	/// How a matrix is really layed out is pretty much abitrary but we need to agree
	/// and the world has agreed (except for Microsoft) on the right-hand rule. 
	/// First, the 4x4 matrix is really just an array of 16 numbers.  
	/// We need to think of the array as a matrix in the following way
	/// 4x4 matrix - COLUMN MAJOR - The OpenGL, science, physics, mathematics and engineering way. 
	///	0  4  8  12        [0][0]  [1][0]  [2][0]  [3][0]   
	///	1  5  9  13  (or)  [0][1]  [1][1]  [2][1]  [3][1]
	///	2  6  10 14        [0][2]  [1][2]  [2][2]  [3][2]
	///	3  7  11 15        [0][3]  [1][3]  [2][3]  [3][3]

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Matrix4
		{
			~Matrix4() = default;

			explicit Matrix4(T s = 0)
			{
				if(s == 0)
				{
					SetIdentity();
				}
				else
				{
					m[0] = s;	m[1] = s;	m[2] = s;   m[3] = s;
					m[4] = s;	m[5] = s;   m[6] = s;	m[7] = s;
					m[8] = s;   m[9] = s;   m[10] = s;	m[11] = s;
					m[12] = s;   m[13] = s;   m[14] = s;	m[15] = s;
				}
			}

			Matrix4(const Matrix4<T>& _m)
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];   m[3] = _m[3];
				m[4] = _m[4];   m[5] = _m[5];   m[6] = _m[6];   m[7] = _m[7];
				m[8] = _m[8];   m[9] = _m[9];   m[10] = _m[10];   m[11] = _m[11];
				m[12] = _m[12];   m[13] = _m[13];   m[14] = _m[14]; m[15] = _m[15];
			}

			Matrix4(Matrix4<T>&& _m) noexcept
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];   m[3] = _m[3];
				m[4] = _m[4];   m[5] = _m[5];   m[6] = _m[6];   m[7] = _m[7];
				m[8] = _m[8];   m[9] = _m[9];   m[10] = _m[10];   m[11] = _m[11];
				m[12] = _m[12];   m[13] = _m[13];   m[14] = _m[14]; m[15] = _m[15];

				_m[0] = T();     _m[1] = T();     _m[2] = T();     _m[3] = T();
				_m[4] = T();    _m[5] = T();    _m[6] = T();   _m[7] = T();
				_m[8] = T();     _m[9] = T();     _m[10] = T();    _m[11] = T();
				_m[12] = T();    _m[13] = T();    _m[14] = T();    _m[15] = T();
			}

			Matrix4(
				T m00, T m01, T m02, T m03,
				T m10, T m11, T m12, T m13,
				T m20, T m21, T m22, T m23,
				T m30, T m31, T m32, T m33)
			{
				m[0] = m00; m[1] = m01; m[2] = m02; m[3] = m03;
				m[4] = m10; m[5] = m11; m[6] = m12; m[7] = m13;
				m[8] = m20; m[9] = m21; m[10] = m22; m[11] = m23;
				m[12] = m30; m[13] = m31; m[14] = m32; m[15] = m33;
			}

			Matrix4(const Vector4<T>& v)
			{
				m[0] = v.x; m[1] = 0.0; m[2] = 0.0; m[3] = 0.0;
				m[4] = 0.0; m[5] = v.y; m[6] = 0.0; m[7] = 0.0;
				m[8] = 0.0; m[9] = 0.0; m[10] = v.z; m[11] = 0.0;
				m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = v.w;
			}

			Matrix4& operator=(const Matrix4<T>& _m)
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];   m[3] = _m[3];
				m[4] = _m[4];   m[5] = _m[5];   m[6] = _m[6];   m[7] = _m[7];
				m[8] = _m[8];   m[9] = _m[9];   m[10] = _m[10];   m[11] = _m[11];
				m[12] = _m[12];   m[13] = _m[13];   m[14] = _m[14]; m[15] = _m[15];
				return *this;
			}

			Matrix4& operator=(Matrix4<T>&& _m) noexcept
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];   m[3] = _m[3];
				m[4] = _m[4];   m[5] = _m[5];   m[6] = _m[6];   m[7] = _m[7];
				m[8] = _m[8];   m[9] = _m[9];   m[10] = _m[10];   m[11] = _m[11];
				m[12] = _m[12];   m[13] = _m[13];   m[14] = _m[14]; m[15] = _m[15];

				_m[0] = T();     _m[1] = T();     _m[2] = T();     _m[3] = T();
				_m[4] = T();    _m[5] = T();    _m[6] = T();   _m[7] = T();
				_m[8] = T();     _m[9] = T();     _m[10] = T();    _m[11] = T();
				_m[12] = T();    _m[13] = T();    _m[14] = T();    _m[15] = T();
				return *this;
			}

			bool operator==(const Matrix4<T>& _m) const
			{
				return (m == _m.m);
			}

			bool operator!=(const Matrix4<T>& _m) const
			{
				return !(*this == _m);
			}

			T operator[](const unsigned int index) const
			{
				return *(m + index);
			}

			T& operator[](const unsigned int index)
			{
				return *(m + index);
			}

			Matrix4 operator*(const T& s) const
			{
				Matrix4 result = *this;
				result.m[0] *= s;   result.m[1] *= s;   result.m[2] *= s;   result.m[3] *= s;
				result.m[4] *= s;   result.m[5] *= s;   result.m[6] *= s;   result.m[7] *= s;
				result.m[8] *= s;   result.m[9] *= s;   result.m[10] *= s;  result.m[11] *= s;
				result.m[12] *= s;  result.m[13] *= s;  result.m[14] *= s;  result.m[15] *= s;
				return result;
			}

			Matrix4 operator*(const Matrix4& _m) const
			{
				Matrix4 result;
				for(int i = 0; i < 4; ++i)
				{
					for(int k = 0; k < 4; ++k)
					{
						result[i * 4 + k] = (m[0 * 4 + k] * _m[i * 4 + 0]) + (m[1 * 4 + k] * _m[i * 4 + 1]) + (m[2 * 4 + k] * _m[i * 4 + 2]) + (m[3 * 4 + k] * _m[i * 4 + 3]);
					}
				}
				return result;
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * matrix = product;
			static friend Matrix4 operator*(const T& s, const Matrix4& _m)
			{
				return _m * s;
			}

			Matrix4& operator*=(const T& s)
			{
				*this = *this * s;
				return *this;
			}

			Matrix4& operator*=(const Matrix4& _m)
			{
				*this = *this * _m;
				return *this;
			}

			Vector4<T> operator*(const Vector4<T>& v) const
			{
				const T x = v.x * m[0] + v.y * m[4] + v.z * m[8] + v.w * m[12];
				const T y = v.x * m[1] + v.y * m[5] + v.z * m[9] + v.w * m[13];
				const T z = v.x * m[2] + v.y * m[6] + v.z * m[10] + v.w * m[14];
				const T w = v.x * m[3] + v.y * m[7] + v.z * m[11] + v.w * m[15];
				return Vector4<T>(x, y, z, w);
			}

			void SetIdentity()
			{
				m[0] = 1;  m[1] = 0;  m[2] = 0;  m[3] = 0;
				m[4] = 0;  m[5] = 1;  m[6] = 0;  m[7] = 0;
				m[8] = 0;  m[9] = 0;  m[10] = 1; m[11] = 0;
				m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
			}

			void SetIdentity(T s)
			{
				m[0] = s;  m[1] = 0;  m[2] = 0;  m[3] = 0;
				m[4] = 0;  m[5] = s;  m[6] = 0;  m[7] = 0;
				m[8] = 0;  m[9] = 0;  m[10] = s; m[11] = 0;
				m[12] = 0; m[13] = 0; m[14] = 0; m[15] = s;
			}

			Vector4<T> GetRow(const unsigned int index) const
			{
				return Vector4<T>(m[4 * index + 0], m[4 * index + 1], m[4 * index + 2], m[4 * index + 3]);
			}

			Vector4<T> GetColumn(const unsigned int index) const
			{
				return Vector4<T>(m[0 + index], m[4 + index], m[8 + index], m[12 + index]);
			}

			T Determinant() const
			{
				return
					m[3] * m[6] * m[9] * m[12] - m[2] * m[7] * m[9] * m[12] - m[3] * m[5] * m[10] * m[12] + m[1] * m[7] * m[10] * m[12] +
					m[2] * m[5] * m[11] * m[12] - m[1] * m[6] * m[11] * m[12] - m[3] * m[6] * m[8] * m[13] + m[2] * m[7] * m[8] * m[13] +
					m[3] * m[4] * m[10] * m[13] - m[0] * m[7] * m[10] * m[13] - m[2] * m[4] * m[11] * m[13] + m[0] * m[6] * m[11] * m[13] +
					m[3] * m[5] * m[8] * m[14] - m[1] * m[7] * m[8] * m[14] - m[3] * m[4] * m[9] * m[14] + m[0] * m[7] * m[9] * m[14] +
					m[1] * m[4] * m[11] * m[14] - m[0] * m[5] * m[11] * m[14] - m[2] * m[5] * m[8] * m[15] + m[1] * m[6] * m[8] * m[15] +
					m[2] * m[4] * m[9] * m[15] - m[0] * m[6] * m[9] * m[15] - m[1] * m[4] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
			}

			Matrix4<T> Inverse() const
			{
				const T determinant = Determinant();
				// When the determinant is less than zero, it means the Matrix has been reflected
				// When the determinant is zero then this Matrix has been projected, meaning that a column or row is filled with zeroes
				// When the determinant is one than this Matrix is an identity Matrix
				// With any other determinant result, this inverse matrix will reverse a rotation, transformation, scale, etc impacting said matrix 
#ifdef _DEBUG  
				if(fabs(determinant) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				Matrix4 inverse;

				inverse[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
				inverse[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
				inverse[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
				inverse[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
				inverse[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
				inverse[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
				inverse[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
				inverse[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
				inverse[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
				inverse[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
				inverse[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
				inverse[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
				inverse[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
				inverse[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
				inverse[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
				inverse[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

				determinant = 1.0f / determinant;

				for(int i = 0; i < 16; i++)
				{
					inverse[i] *= determinant;
				}
				return inverse;
			}

			Matrix4<T> Transpose() const
			{
				return Matrix4(
					m[0], m[4], m[8], m[12],
					m[1], m[5], m[9], m[13],
					m[2], m[6], m[10], m[14],
					m[3], m[7], m[11], m[15]
				);
			}

			std::string ToString() const;

			static Matrix4<T> Translate(const Vector3<T>& translation);
			static Matrix4<T> Scale(const Vector3<T>& scale);
			static Matrix4<T> Rotate(T angle, const Vector3<T>& axis);

			static Matrix4<T> PerspectiveMatrix(T fov, T aspectRatio, T zNear, T zFar);
			static Matrix4<T> OrthographicMatrix(T left, T right, T bottom, T top, T zNear, T zFar);
			static Matrix4<T> UnorthoMatrix(const Matrix4& orthoMatrix, T left, T right, T bottom, T top, T zNear, T zFar);
			static Matrix4<T> LookAt(const Vector3<T>& _eye, const Vector3<T>& _at, const Vector3<T>& _up);

			static Matrix4 ViewportNormalizedDeviceCoordinates(T width, T height);

			static const Matrix4 Identity;

		private:
			T m[16];
		};

		template<typename T>
		std::string Matrix4<T>::ToString() const
		{
			std::ostringstream oss;
			for (int i = 0; i < 4; ++i) {
				oss << "(";
				for (int j = 0; j < 4; ++j) {
					oss << m[i * 4 + j];
					if (j < 3) oss << ", ";
				}
				oss << ")\n";
			}
			return oss.str();
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::Translate(const Vector3<T>& translation)
		{
			return Matrix4(
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				translation.x, translation.y, translation.z, 1
			);
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::Scale(const Vector3<T>& scale)
		{
			return Matrix4(
				scale.x, 0, 0, 0,
				0, scale.y, 0, 0,
				0, 0, scale.z, 0,
				0, 0, 0, 1
			);
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::Rotate(T angle, const Vector3<T>& axis)
		{
			Vector3<T> normalizedAxis = axis.Normal();
			T radians = static_cast<T>(angle * DEGREES_TO_RADIANS);    // Convert degrees to radians

			T cosAngle = static_cast<T>(cos(radians));
			T sinAngle = static_cast<T>(sin(radians));
			T oneMinusCos = static_cast<T>(1.0f - cosAngle);

			cosAngle = fabs(cosAngle) < NEARLY_ZERO ? 0 : cosAngle;
			sinAngle = fabs(sinAngle) < NEARLY_ZERO ? 0 : sinAngle;

			Matrix4<T> result;
			result[0] = (normalizedAxis.x * normalizedAxis.x * oneMinusCos) + cosAngle;
			result[1] = (normalizedAxis.x * normalizedAxis.y * oneMinusCos) - (normalizedAxis.z * sinAngle);
			result[2] = (normalizedAxis.x * normalizedAxis.z * oneMinusCos) + (normalizedAxis.y * sinAngle);
			result[3] = 0;
			result[4] = (normalizedAxis.y * normalizedAxis.x * oneMinusCos) + (normalizedAxis.z * sinAngle);
			result[5] = (normalizedAxis.y * normalizedAxis.y * oneMinusCos) + cosAngle;
			result[6] = (normalizedAxis.y * normalizedAxis.z * oneMinusCos) - (normalizedAxis.x * sinAngle);
			result[7] = 0;
			result[8] = (normalizedAxis.z * normalizedAxis.x * oneMinusCos) - (normalizedAxis.y * sinAngle);
			result[9] = (normalizedAxis.z * normalizedAxis.y * oneMinusCos) + (normalizedAxis.x * sinAngle);
			result[10] = (normalizedAxis.z * normalizedAxis.z * oneMinusCos) + cosAngle;
			result[11] = 0;
			result[12] = 0;
			result[13] = 0;
			result[14] = 0;
			result[15] = 1;
			return result;
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::PerspectiveMatrix(T fovY, T aspectRatio, T zNear, T zFar)
		{
			// Aspect ratio is calculated by dividing the height of the screen by its width in pixels
			// Fov along the y axis is usually 90 degrees for humans
			// zNear defines the position of out near plane and zFar the position of the far plane.

			T cot = static_cast<T>(1.0 / tan(fovY * 0.5 * DEGREES_TO_RADIANS));

			// Keep in mind we are using the right-hand rule
			// We are using the Matrix to create the View Frustum with a Perspective Projection
			// To be then projected on the Normalized Device Coordinates becase that is what is used by OpenGL
			Matrix4 result(
				cot / aspectRatio, 0, 0, 0,
				0, cot, 0, 0,
				0, 0, (zNear + zFar) / (zNear - zFar), -1,
				0, 0, (static_cast<T>(2) * zNear * zFar) / (zNear - zFar), 0
			);
			return result;
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::OrthographicMatrix(T left, T right, T bottom, T top, T zNear, T zFar)
		{
			// Derived-from: https://en.wikipedia.org/wiki/Orthographic_projection

			Matrix4 m1 = Matrix4::Scale(Vector3<T>(
				static_cast<T>(2) / (right - left),
				static_cast<T>(2) / (top - bottom),
				static_cast<T>(-2) / (zFar - zNear)));

			Matrix4 m2 = Matrix4::Translate(Vector3<T>(
				-(right + left) / (right - left),
				-(top + bottom) / (top - bottom),
				-(zFar + zNear) / (zFar - zNear)));

			Matrix4 orthographicMatrix = m2 * m1;
			return orthographicMatrix;
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::UnorthoMatrix(const Matrix4& orthoMatrix, T left, T right, T bottom, T top, T zNear, T zFar)
		{
			// Derived-from: https://en.wikipedia.org/wiki/Orthographic_projection
			// Inverting the projection matrix, will unortho

			Matrix4 m1 = Matrix4::Scale(Vector3<T>(
				(right - left) / static_cast<T>(2),
				(top - bottom) / static_cast<T>(2),
				(zFar - zNear) / static_cast<T>(-2)));

			Matrix4 m2 = Matrix4::Translate(Vector3<T>(
				(left + right) / static_cast<T>(2), 
				(top + bottom) / static_cast<T>(2), 
				-(zFar + zNear) / static_cast<T>(2)));

			Matrix4 unorthMatrix = m2 * m1;
			return unorthMatrix;
		}

		template <typename T>
		Matrix4<T> Matrix4<T>::LookAt(const Vector3<T>& _eye, const Vector3<T>& _at, const Vector3<T>& _up)
		{
			Vector3<T> at = _at;
			Vector3<T> up = _up;
			Vector3<T> eye = _eye;

			Matrix4 result;

			Vector3<T> forward = Vector3<T>::Normalize(at - eye);   // Creating the direction vector which is the -eyeDirec currently
			up = up.Normal();
			Vector3<T> side = Vector3<T>::Normalize(Vector3<T>::CrossProduct(forward, up));    // Creating the side of the camera
			up = Vector3<T>::CrossProduct(side, forward);   // Then using the side and forward vector's we can make the new up direction


			// Below is the Matrix to used represent the Eye(LookAt)
			result[0] = side.x;
			result[1] = side.y;
			result[2] = side.z;
			result[3] = static_cast<T>(0);

			result[4] = up.x;
			result[5] = up.y;
			result[6] = up.z;
			result[7] = static_cast<T>(0);

			result[8] = -forward.x;
			result[9] = -forward.y;
			result[10] = -forward.z;
			result[11] = static_cast<T>(0);

			result[12] = -Vector3<T>::DotProduct(side, eye);
			result[13] = -Vector3<T>::DotProduct(up, eye);
			result[14] = Vector3<T>::DotProduct(forward, eye);
			result[15] = static_cast<T>(1);

			return result;
		}

		/// This creates a transform from Normalized Device Coordinates (NDC) to 
		/// screen coordinates. OpenGL uses NDC as	(Left-Hand Rule)		 
		///	              ------------------------------
		///	             /|                           /|
		///	            / |                          / |
		///	           /  |                         /  |
		///	          /   |                        /   |
		///	         /    |                       /    |
		///	        /     |                      /     |
		///	       /      |                     /      |
		///	      /       |                    /       |
		///	     /        |                   /        |
		///	    /         |                  /         |
		///	   /----------------------------/ (1.0,1.0)|
		///	   |          |                 |          |
		///	   |          |                 |          |      +Y
		///	   |          |                 |          | 
		///	   |          |-----------------|----------|      ^
		///	   |         /                  |         /       |
		///	   |        /                   |        /        |       -Z
		///	   |       /                    |       /         |
		///	   |      /                     |      /          |     /
		///	   |     /                      |     /           |    /
		///	   |    /                       |    /            |   /
		///	   |   /                        |   /             |  /
		///	   |  /                         |  /              | /
		///	   | /                          | /               |/
		///	   |/ (-1.0,-1.0)               |/                ----------------> +X
		///	   ------------------------------
		template <typename T>
		Matrix4<T> Matrix4<T>::ViewportNormalizedDeviceCoordinates(T width, T height)
		{
			// When using the NDC the z that we use will not ever be negative, it will be between 0 and 1, anything outside that will not be drawn.

			const T minZ = 0;
			const T maxZ = 1;

			Matrix4 m1 = Matrix4::Scale(1, -1, 1);
			Matrix4 m2 = Matrix4::Scale(width / static_cast<T>(2), height / static_cast<T>(2), maxZ - minZ);
			Matrix4 m3 = Matrix4::Translate(width / static_cast<T>(2), height / static_cast<T>(2), minZ);
			Matrix4 result = m3 * m2 * m1;

			/***
			 * Above scale and translation accomplish the below code

				result[0] = float(width_)/2.0f;
				result[5] = -float(height_)/2.0f;
				result[10] =  maxZ - minZ;
				result[12] = float(width_)/2.0f;
				result[13] = float(height_)/2.0f;
				result[14] = minZ;
				result[15] = 1.0f;
			***/

			return result;
		}

		template<typename T>
		const Matrix4<T> Matrix4<T>::Identity = Matrix4<T>(
			1, 0, 0, 0, 
			0, 1, 0, 0, 
			0, 0, 1, 0, 
			0, 0, 0, 1);


		// ~~~~~~~~~~~~~~~~~ Unit Tests ~~~~~~~~~~~~~~~~~~~~~~~~

#if _DEBUG

	// Helper function to print matrix
		template<typename T>
		void PrintMatrix(const Matrix4<T>& matrix, const std::string& name)
		{
			std::cout << name << ":\n";
			std::cout << matrix.ToString().c_str() << std::endl;
		}

		// Helper function to print matrix and compare with expected values
		template<typename T>
		void PrintAndCompareMatrix(const Matrix4<T>& matrix, const std::string& name, const Matrix4<T>& expected)
		{
			std::cout << name << ":\n";
			bool pass = true;
			std::cout << "Result Matrix:\n";
			for(int i = 0; i < 4; ++i)
			{
				for(int j = 0; j < 4; ++j)
				{
					if(matrix[i * 4 + j] != expected[i * 4 + j])
					{
						pass = false;
					}
				}
			}
			std::cout << matrix.ToString();
			std::cout << "Expected Matrix:\n";
			std::cout << expected.ToString();
			std::cout << (pass ? "Pass" : "Fail") << "\n\n";
		}

		template<typename T>
		void TestMatrix4Type()
		{
			// Test default constructor
			Matrix4<T> mat1;
			Matrix4<T> expected1;
			expected1.SetIdentity();
			PrintAndCompareMatrix(mat1, "Default Constructor", expected1);

			// Test constructor with scalar
			Matrix4<T> mat2(static_cast<T>(5));
			Matrix4<T> expected2(
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5), static_cast<T>(5),
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5), static_cast<T>(5),
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5), static_cast<T>(5),
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5), static_cast<T>(5)
			);
			PrintAndCompareMatrix(mat2, "Constructor with Scalar", expected2);

			// Test copy constructor
			Matrix4<T> mat3 = mat2;
			PrintAndCompareMatrix(mat3, "Copy Constructor", expected2);

			// Test move constructor
			Matrix4<T> mat4 = std::move(mat3);
			PrintAndCompareMatrix(mat4, "Move Constructor", expected2);

			// Test parameterized constructor
			Matrix4<T> mat5(
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(4),
				static_cast<T>(5), static_cast<T>(6), static_cast<T>(7), static_cast<T>(8),
				static_cast<T>(9), static_cast<T>(10), static_cast<T>(11), static_cast<T>(12),
				static_cast<T>(13), static_cast<T>(14), static_cast<T>(15), static_cast<T>(16)
			);
			Matrix4<T> expected5(
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(4),
				static_cast<T>(5), static_cast<T>(6), static_cast<T>(7), static_cast<T>(8),
				static_cast<T>(9), static_cast<T>(10), static_cast<T>(11), static_cast<T>(12),
				static_cast<T>(13), static_cast<T>(14), static_cast<T>(15), static_cast<T>(16)
			);
			PrintAndCompareMatrix(mat5, "Parameterized Constructor", expected5);

			// Test assignment operator
			Matrix4<T> mat6;
			mat6 = mat5;
			PrintAndCompareMatrix(mat6, "Assignment Operator", expected5);

			// Test scalar multiplication
			Matrix4<T> mat7 = mat5 * static_cast<T>(2);
			Matrix4<T> expected7(
				static_cast<T>(2), static_cast<T>(4), static_cast<T>(6), static_cast<T>(8),
				static_cast<T>(10), static_cast<T>(12), static_cast<T>(14), static_cast<T>(16),
				static_cast<T>(18), static_cast<T>(20), static_cast<T>(22), static_cast<T>(24),
				static_cast<T>(26), static_cast<T>(28), static_cast<T>(30), static_cast<T>(32)
			);
			PrintAndCompareMatrix(mat7, "Scalar Multiplication", expected7);

			// Test matrix multiplication
			Matrix4<T> mat8(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(1), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(1)
			);
			Matrix4<T> mat9 = mat5 * mat8;
			PrintAndCompareMatrix(mat9, "Matrix Multiplication", expected5);

			// Test determinant
			Matrix4<T> mat10(
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(4),
				static_cast<T>(5), static_cast<T>(6), static_cast<T>(7), static_cast<T>(8),
				static_cast<T>(9), static_cast<T>(10), static_cast<T>(11), static_cast<T>(12),
				static_cast<T>(13), static_cast<T>(14), static_cast<T>(15), static_cast<T>(16)
			);
			T determinant = mat10.Determinant();
			std::cout << "Determinant: " << determinant << " | Expected: 0\n";
			std::cout << (determinant == static_cast<T>(0) ? "Pass" : "Fail") << std::endl;

			// Test transpose
			Matrix4<T> mat11 = mat10.Transpose();
			Matrix4<T> expected11(
				static_cast<T>(1), static_cast<T>(5), static_cast<T>(9), static_cast<T>(13),
				static_cast<T>(2), static_cast<T>(6), static_cast<T>(10), static_cast<T>(14),
				static_cast<T>(3), static_cast<T>(7), static_cast<T>(11), static_cast<T>(15),
				static_cast<T>(4), static_cast<T>(8), static_cast<T>(12), static_cast<T>(16)
			);
			PrintAndCompareMatrix(mat11, "Transpose", expected11);

			// Test scale
			Vector3<T> scaleFactors(static_cast<T>(2), static_cast<T>(3), static_cast<T>(4));
			Matrix4<T> scaleMatrix = Matrix4<T>::Scale(scaleFactors);
			Matrix4<T> expectedScaleMatrix(
				static_cast<T>(2), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(3), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(4), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(1)
			);
			PrintAndCompareMatrix(scaleMatrix, "Scale", expectedScaleMatrix);

			// Test translate
			Vector3<T> translationFactors(static_cast<T>(1), static_cast<T>(2), static_cast<T>(3));
			Matrix4<T> translateMatrix = Matrix4<T>::Translate(translationFactors);
			Matrix4<T> expectedTranslateMatrix(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(1), static_cast<T>(0),
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(1)
			);
			PrintAndCompareMatrix(translateMatrix, "Translate", expectedTranslateMatrix);

			// Test rotate
			Vector3<T> rotationAxis(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0));
			Matrix4<T> rotateMatrix = Matrix4<T>::Rotate(90, rotationAxis);
			Matrix4<T> expectedRotateMatrix(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-1), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(0), static_cast<T>(1)
			);
			PrintAndCompareMatrix(rotateMatrix, "Rotate", expectedRotateMatrix);

			// Test perspective matrix
			Matrix4<T> perspectiveMatrix = Matrix4<T>::PerspectiveMatrix(static_cast<T>(45), static_cast<T>(1.0), static_cast<T>(0.1), static_cast<T>(100));
			// Expected values can be pre-calculated based on the implementation details
			Matrix4<T> expectedPerspectiveMatrix(
				static_cast<T>(2.41421356237), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(2.41421356237), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-1.002002002), static_cast<T>(-1),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-0.2002002002), static_cast<T>(0)
			);
			PrintAndCompareMatrix(perspectiveMatrix, "PerspectiveMatrix", expectedPerspectiveMatrix);

			// Test orthographic matrix
			Matrix4<T> orthographicMatrix = Matrix4<T>::OrthographicMatrix(static_cast<T>(-1), static_cast<T>(1), static_cast<T>(-1), static_cast<T>(1), static_cast<T>(0.1), static_cast<T>(100));
			Matrix4<T> expectedOrthographicMatrix(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-0.0202020202), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-1.002002002), static_cast<T>(1)
			);
			PrintAndCompareMatrix(orthographicMatrix, "OrthographicMatrix", expectedOrthographicMatrix);

			// Test unorthographic matrix
			Matrix4<T> unorthoMatrix = Matrix4<T>::UnorthoMatrix(orthographicMatrix, static_cast<T>(-1), static_cast<T>(1), static_cast<T>(-1), static_cast<T>(1), static_cast<T>(0.1), static_cast<T>(100));
			Matrix4<T> expectedUnorthoMatrix(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-49.5), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-50.5), static_cast<T>(1)
			);
			PrintAndCompareMatrix(unorthoMatrix, "UnorthoMatrix", expectedUnorthoMatrix);

			// Test look at matrix
			Vector3<T> eye(static_cast<T>(0), static_cast<T>(0), static_cast<T>(5));
			Vector3<T> at(static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
			Vector3<T> up(static_cast<T>(0), static_cast<T>(1), static_cast<T>(0));
			Matrix4<T> lookAtMatrix = Matrix4<T>::LookAt(eye, at, up);
			Matrix4<T> expectedLookAtMatrix(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(1), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(-5), static_cast<T>(1)
			);
			PrintAndCompareMatrix(lookAtMatrix, "LookAt", expectedLookAtMatrix);
		}

		inline void TestMatrix4()
		{
			std::cout << "Testing Matrix4 with float type:\n";
			TestMatrix4Type<float>();

			std::cout << "\nTesting Matrix4 with double type:\n";
			TestMatrix4Type<double>();

			std::cout << "\nTesting Matrix4 with int type:\n";
			TestMatrix4Type<int>();
		}

#endif // _DEBUG

	}
}

#endif // MJOLNIR_MATRIX4_H