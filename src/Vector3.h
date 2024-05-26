// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_VECTOR3_H
#define MJOLNIR_VECTOR3_H

#include "MjolnirConstants.h"

#if _DEBUG
#include <stdexcept>
#include <iostream>
#endif // _DEBUG 

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Vector3
		{
			T x, y, z;

			~Vector3() = default;

			explicit Vector3(T s = 0):
				x(s),
				y(s),
				z(s)
			{}

			Vector3(T _x, T _y, T _z):
				x(_x),
				y(_y),
				z(_z)
			{}

			Vector3(const Vector3& v):
				x(v.x),
				y(v.y),
				z(v.z)
			{}

			Vector3(Vector3&& v) noexcept:
				x(v.x),
				y(v.y),
				z(v.z)
			{
				v.x = v.y = v.z = T();
			}

			Vector3& operator=(const Vector3& v)
			{
				x = v.x;
				y = v.y;
				z = v.z;
				return *this;
			}

			Vector3& operator=(Vector3&& v) noexcept
			{
				if(this != &v)
				{
					x = v.x;
					y = v.y;
					z = v.z;
					v.x = v.y = v.z = T();
				}
				return *this;
			}

			bool operator==(const Vector3<T>& v) const
			{
				return (x == v.x) && (y == v.y) && (z == v.z);
			}

			bool operator!=(const Vector3<T>& v) const
			{
				return !(*this == v);
			}

			T operator[](const unsigned int index) const
			{
				return *(&x + index);
			}

			T& operator[](const unsigned int index)
			{
				return *(&x + index);
			}

			Vector3 operator-() const
			{
				return Vector3(-x, -y, -z);
			}

			Vector3 operator+(const Vector3& v) const
			{
				return Vector3(
					x + v.x,
					y + v.y,
					z + v.z
				);
			}

			Vector3 operator-(const Vector3& v) const
			{
				return Vector3(
					x - v.x,
					y - v.y,
					z - v.z
				);
			}

			Vector3 operator+=(const Vector3& v)
			{
				x += v.x;
				y += v.y;
				z += v.z;
				return *this;
			}

			Vector3 operator-=(const Vector3& v)
			{
				x -= v.x;
				y -= v.y;
				z -= v.z;
				return *this;
			}

			Vector3 operator*(const T& s) const
			{
				return Vector3(
					x * s,
					y * s,
					z * s
				);
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * vector = product;
			static friend Vector3 operator*(const T& s, const Vector3& v)
			{
				return v * s;
			}

			Vector3 operator*=(const T& s)
			{
				x *= s;
				y *= s;
				z *= s;
				return *this;
			}

			Vector3 operator/(const T& s) const
			{
#ifdef _DEBUG
				if(std::fabs(s) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				return Vector3(
					x / s,
					y / s,
					z / s
				);
			}

			Vector3 operator/=(const T& s)
			{
#ifdef _DEBUG
				if(std::fabs(s) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				x /= s;
				y /= s;
				z /= s;
				return *this;
			}

			T Magnitude() const;
			Vector3<T> Normal() const;
			void Normalize();
			std::string ToString() const;

			static Vector3<T> Normalize(const Vector3<T>& v);
			static T Distance(const Vector3& v1, const Vector3& v2);
			static T DotProduct(const Vector3& v1, const Vector3& v2);
			static T Angle(const Vector3& v1, const Vector3& v2);
			static Vector3 CrossProduct(const Vector3& v1, const Vector3& v2);
			static Vector3 Lerp(const Vector3& v1, const Vector3& v2, const T t);
			static Vector3 ProjectVectorOnToVector(const Vector3& project, const Vector3& target);

			static const Vector3 ZeroVector;
			static const Vector3 UpVector;
			static const Vector3 RightVector;
			static const Vector3 DownVector;
			static const Vector3 LeftVector;
			static const Vector3 ForwardVector;
			static const Vector3 BackVector;
		};

		template <typename T>
		T Vector3<T>::Magnitude() const
		{
			return static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)));
		}

		template <typename T>
		Vector3<T> Vector3<T>::Normal() const
		{
			T a = static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)));
#ifdef _DEBUG
			if(fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			return Vector3(x / a, y / a, z / a);
		}

		template <typename T>
		void Vector3<T>::Normalize()
		{
			T a = static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)));
#ifdef _DEBUG
			if(fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			x /= a;
			y /= a;
			z /= a;
		}

		template<typename T>
		std::string Vector3<T>::ToString() const
		{
			std::ostringstream oss;
			oss << "(" << x << ", " << y << ", " << z << ")";
			return oss.str();
		}

		template <typename T>
		Vector3<T> Vector3<T>::Normalize(const Vector3<T>& v)
		{
			return v.Normal();
		}

		template <typename T>
		T Vector3<T>::Distance(const Vector3& v1, const Vector3& v2)
		{
			const Vector3 diff(v1 - v2);
			return diff.Magnitude();
		}

		template <typename T>
		T Vector3<T>::DotProduct(const Vector3& v1, const Vector3& v2)
		{
			return static_cast<T>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
		}

		template <typename T>
		T Vector3<T>::Angle(const Vector3& v1, const Vector3& v2)
		{
			T a = DotProduct(v1, v2);
			T b = static_cast<T>(sqrt(v1.Magnitude() * v2.Magnitude()));	// Multiplying the magnitude of each vector
			return  static_cast<T>(acos(a / b));	// Taking the cos inverse of the dot product / magnitude result to get the angle in radians
		}

		template <typename T>
		Vector3<T> Vector3<T>::CrossProduct(const Vector3& v1, const Vector3& v2)
		{
			return Vector3(
				v1.y * v2.z - v1.z * v2.y,
				-(v1.x * v2.z - v1.z * v2.x),
				v1.x * v2.y - v1.y * v2.x
			);
		}

		template <typename T>
		Vector3<T> Vector3<T>::Lerp(const Vector3& v1, const Vector3& v2, const T t)
		{
			return Vector3(v1 + t * (v2 - v1));
		}

		template <typename T>
		Vector3<T> Vector3<T>::ProjectVectorOnToVector(const Vector3& project, const Vector3& target)
		{
			T mag = target.Magnitude();
#ifdef _DEBUG
			if(fabs(mag) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero in vector projection!"); // Checking the magnitude of the target vector
			}
#endif
			Vector3 normal(target.x / mag, target.y / mag, target.z / mag);     // Normal of the target vector
			T scalar = DotProduct(project, normal);     // Dot product of the vector to project and the normalized target vector
			return normal * scalar;     // Returning the normal scaled to the length of the projected vector's shadow
		}

		template<typename T>
		const Vector3<T> Vector3<T>::ZeroVector = Vector3<T>(0, 0, 0);

		template<typename T>
		const Vector3<T> Vector3<T>::UpVector = Vector3<T>(0, 1, 0);

		template<typename T>
		const Vector3<T> Vector3<T>::RightVector = Vector3<T>(1, 0, 0);

		template<typename T>
		const Vector3<T> Vector3<T>::DownVector = Vector3<T>(0, -1, 0);

		template<typename T>
		const Vector3<T> Vector3<T>::LeftVector = Vector3<T>(-1, 0, 0);

		template<typename T>
		const Vector3<T> Vector3<T>::ForwardVector = Vector3<T>(0, 0, 1);

		template<typename T>
		const Vector3<T> Vector3<T>::BackVector = Vector3<T>(0, 0, -1);



		// ~~~~~~~~~~~~~~~~~ Unit Tests ~~~~~~~~~~~~~~~~~~~~~~

#if _DEBUG
		template<typename T>
		void PrintVector(const Math::Vector3<T>& v, const std::string& name, const Math::Vector3<T>& expected)
		{
			std::cout << name << ":" << v.ToString() << ""
				<< " | Expected: " << expected.ToString() << "\n";
		}

		template<typename T>
		void PrintScalar(T value, const std::string& name, T expected)
		{
			std::cout << name << ": " << value << " | Expected: " << expected << "\n";
		}

		template<typename T>
		void TestVector3Type(const std::string& typeName)
		{
			std::cout << "Running Vector3 Tests for " << typeName << "\n";

			Math::Vector3<T> v1(static_cast<T>(1), static_cast<T>(2), static_cast<T>(3));
			Math::Vector3<T> v2(static_cast<T>(4), static_cast<T>(5), static_cast<T>(6));

			PrintVector(v1, "v1", Math::Vector3<T>(static_cast<T>(1), static_cast<T>(2), static_cast<T>(3)));
			PrintVector(v2, "v2", Math::Vector3<T>(static_cast<T>(4), static_cast<T>(5), static_cast<T>(6)));

			Math::Vector3<T> sum = v1 + v2;
			PrintVector(sum, "v1 + v2", Math::Vector3<T>(static_cast<T>(5), static_cast<T>(7), static_cast<T>(9)));

			Math::Vector3<T> diff = v1 - v2;
			PrintVector(diff, "v1 - v2", Math::Vector3<T>(static_cast<T>(-3), static_cast<T>(-3), static_cast<T>(-3)));

			Math::Vector3<T> neg = -v1;
			PrintVector(neg, "-v1", Math::Vector3<T>(static_cast<T>(-1), static_cast<T>(-2), static_cast<T>(-3)));

			Math::Vector3<T> scaled = v1 * static_cast<T>(2);
			PrintVector(scaled, "v1 * 2", Math::Vector3<T>(static_cast<T>(2), static_cast<T>(4), static_cast<T>(6)));

			Math::Vector3<T> div = v1 / static_cast<T>(2);
			PrintVector(div, "v1 / 2", Math::Vector3<T>(static_cast<T>(0.5), static_cast<T>(1), static_cast<T>(1.5)));

			T dot = Math::Vector3<T>::DotProduct(v1, v2);
			PrintScalar(dot, "Dot Product (v1 . v2)", static_cast<T>(32));

			Math::Vector3<T> cross = Math::Vector3<T>::CrossProduct(v1, v2);
			PrintVector(cross, "Cross Product (v1 x v2)", Math::Vector3<T>(static_cast<T>(-3), static_cast<T>(6), static_cast<T>(-3)));

			T magnitude = v1.Magnitude();
			PrintScalar(magnitude, "Magnitude of v1", static_cast<T>(sqrt(static_cast<T>(14)))); // sqrt(1^2 + 2^2 + 3^2) = sqrt(14)

			Math::Vector3<T> normal = v1.Normal();
			T invMag = static_cast<T>(1) / magnitude;
			PrintVector(normal, "Normal of v1", Math::Vector3<T>(static_cast<T>(1 * invMag), static_cast<T>(2 * invMag), static_cast<T>(3 * invMag)));

			Math::Vector3<T> lerpResult = Math::Vector3<T>::Lerp(v1, v2, static_cast<T>(0.5));
			PrintVector(lerpResult, "Lerp between v1 and v2 at t = 0.5", Math::Vector3<T>(static_cast<T>(2.5), static_cast<T>(3.5), static_cast<T>(4.5)));

			std::cout << "Vector3 Tests for " << typeName << " Completed\n";
		}

		template<typename T>
		void TestProjection(const std::string& typeName)
		{
			std::cout << "Running Projection Tests for " << typeName << "\n";

			Math::Vector3<T> vectorToProject(static_cast<T>(3), static_cast<T>(4), static_cast<T>(0));
			Math::Vector3<T> targetVector(static_cast<T>(1), static_cast<T>(2), static_cast<T>(2));

			Math::Vector3<T> projectionResult = Math::Vector3<T>::ProjectVectorOnToVector(vectorToProject, targetVector);

			// Calculating expected projection manually
			T mag = static_cast<T>(sqrt(pow(targetVector.x, 2.0) + pow(targetVector.y, 2.0) + pow(targetVector.z, 2.0)));
			Math::Vector3<T> normal(targetVector.x / mag, targetVector.y / mag, targetVector.z / mag);
			T scalar = Math::Vector3<T>::DotProduct(vectorToProject, normal);
			Math::Vector3<T> expectedResult = normal * scalar;

			PrintVector(projectionResult, "Projection of vectorToProject onto targetVector", expectedResult);

			std::cout << "Projection Tests for " << typeName << " Completed\n";
		}

		inline void TestVector3()
		{
			TestVector3Type<float>("float");
			TestVector3Type<double>("double");
			TestVector3Type<int>("int");

			TestProjection<float>("float");
			TestProjection<double>("double");
			TestProjection<int>("int");
		}
#endif // _DEBUG

	}
}

typedef Mjolnir::Math::Vector3<int> Vector3i;
typedef Mjolnir::Math::Vector3<float> Vector3f;
typedef Mjolnir::Math::Vector3<double> Vector3d;

#endif // MJOLNIR_VECTOR3_H