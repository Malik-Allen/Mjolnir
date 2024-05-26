// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_VECTOR4_H
#define MJOLNIR_VECTOR4_H

#include "MjolnirConstants.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Vector4
		{
			T x, y, z, w;

			~Vector4() = default;

			explicit Vector4(T s = 0):
				x(s),
				y(s),
				z(s),
				w(s)
			{}

			Vector4(T _x, T _y, T _z, T _w):
				x(_x),
				y(_y),
				z(_z),
				w(_w)
			{}

			Vector4(const Vector4& v):
				x(v.x),
				y(v.y),
				z(v.z),
				w(v.w)
			{}

			Vector4(Vector4&& v) noexcept:
				x(v.x),
				y(v.y),
				z(v.z),
				w(v.w)
			{
				v.x = v.y = v.z = v.w = T();
			}

			Vector4& operator=(const Vector4& v)
			{
				x = v.x;
				y = v.y;
				z = v.z;
				w = v.w;
				return *this;
			}

			Vector4& operator=(Vector4&& v) noexcept
			{
				if(this != &v)
				{
					x = v.x;
					y = v.y;
					z = v.z;
					w = v.w;
					v.x = v.y = v.z = v.w = T();
				}
				return *this;
			}

			bool operator==(const Vector4<T>& v) const
			{
				return (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w);
			}

			bool operator!=(const Vector4<T>& v) const
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

			Vector4 operator-() const
			{
				return Vector4(-x, -y, -z, -w);
			}

			Vector4 operator+(const Vector4& v) const
			{
				return Vector4(
					x + v.x,
					y + v.y,
					z + v.z,
					w + v.w
				);
			}

			Vector4 operator-(const Vector4& v) const
			{
				return Vector4(
					x - v.x,
					y - v.y,
					z - v.z,
					w - v.w
				);
			}

			Vector4 operator+=(const Vector4& v)
			{
				x += v.x;
				y += v.y;
				z += v.z;
				w += v.w;
				return *this;
			}

			Vector4 operator-=(const Vector4& v)
			{
				x -= v.x;
				y -= v.y;
				z -= v.z;
				w -= v.w;
				return *this;
			}

			Vector4 operator*(const T& s) const
			{
				return Vector4(
					x * s,
					y * s,
					z * s,
					w * s
				);
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * vector = product;
			static friend Vector4 operator*(const T& s, const Vector4& v)
			{
				return v * s;
			}

			Vector4 operator*=(const T& s)
			{
				x *= s;
				y *= s;
				z *= s;
				w *= s;
				return *this;
			}

			Vector4 operator/(const T& s) const
			{
#ifdef _DEBUG
				if(std::fabs(s) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				return Vector4(
					x / s,
					y / s,
					z / s,
					w / s
				);
			}

			Vector4 operator/=(const T& s)
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
				w /= s;
				return *this;
			}

			T Magnitude() const;
			Vector4<T> Normal() const;
			void Normalize();
			std::string ToString() const;

			static Vector4<T> Normalize(const Vector4<T>& v);
			static T Distance(const Vector4& v1, const Vector4& v2);
			static T DotProduct(const Vector4& v1, const Vector4& v2);
			static T Angle(const Vector4& v1, const Vector4& v2);
			static Vector4 Lerp(const Vector4& v1, const Vector4& v2, const T t);
			static Vector4 ProjectVectorOnToVector(const Vector4& project, const Vector4& target);

			static const Vector4 ZeroVector;
			static const Vector4 UpVector;
			static const Vector4 RightVector;
			static const Vector4 DownVector;
			static const Vector4 LeftVector;
			static const Vector4 ForwardVector;
			static const Vector4 BackVector;
			static const Vector4 WVector;
		};

		template <typename T>
		T Vector4<T>::Magnitude() const
		{
			return static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) + pow(w, 2.0)));
		}

		template <typename T>
		Vector4<T> Vector4<T>::Normal() const
		{
			T a = static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) + pow(w, 2.0)));
#ifdef _DEBUG
			if(fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			return Vector4(x / a, y / a, z / a, w / a);
		}

		template <typename T>
		void Vector4<T>::Normalize()
		{
			T a = static_cast<T>(sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) + pow(w, 2.0)));
#ifdef _DEBUG
			if(std::fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			x /= a;
			y /= a;
			z /= a;
			w /= a;
		}

		template<typename T>
		std::string Vector4<T>::ToString() const
		{
			std::ostringstream oss;
			oss << "(" << x << ", " << y << ", " << z << ", " << w << ")";
			return oss.str();
		}

		template <typename T>
		Vector4<T> Vector4<T>::Normalize(const Vector4<T>& v)
		{
			return v.Normal();
		}

		template <typename T>
		T Vector4<T>::Distance(const Vector4& v1, const Vector4& v2)
		{
			const Vector4 diff(v1 - v2);
			return diff.Magnitude();
		}

		template <typename T>
		T Vector4<T>::DotProduct(const Vector4& v1, const Vector4& v2)
		{
			return static_cast<T>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w);
		}

		template <typename T>
		T Vector4<T>::Angle(const Vector4& v1, const Vector4& v2)
		{
			T a = DotProduct(v1, v2);
			T b = static_cast<T>(sqrt(v1.Magnitude() * v2.Magnitude()));	// Multiplying the magnitude of each vector
			return  static_cast<T>(acos(a / b));	// Taking the cos inverse of the dot product / magnitude result to get the angle in radians
		}

		template <typename T>
		Vector4<T> Vector4<T>::Lerp(const Vector4& v1, const Vector4& v2, const T t)
		{
			return Vector4(v1 + t * (v2 - v1));
		}

		template <typename T>
		Vector4<T> Vector4<T>::ProjectVectorOnToVector(const Vector4& project, const Vector4& target)
		{
			T mag = target.Magnitude();
#ifdef _DEBUG
			if(fabs(mag) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero in vector projection!"); // Checking the magnitude of the target vector
			}
#endif
			Vector4 normal(target.x / mag, target.y / mag, target.z / mag, target.w / mag);     // Normal of the target vector
			T scalar = DotProduct(project, normal);     // Dot product of the vector to project and the normalized target vector
			return normal * scalar;     // Returning the normal scaled to the length of the projected vector's shadow
		}

		template<typename T>
		const Vector4<T> Vector4<T>::ZeroVector = Vector4<T>(0, 0, 0, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::UpVector = Vector4<T>(0, 1, 0, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::RightVector = Vector4<T>(1, 0, 0, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::DownVector = Vector4<T>(0, -1, 0, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::LeftVector = Vector4<T>(-1, 0, 0, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::ForwardVector = Vector4<T>(0, 0, 1, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::BackVector = Vector4<T>(0, 0, -1, 0);

		template<typename T>
		const Vector4<T> Vector4<T>::WVector = Vector4<T>(0, 0, 0, 1);

#if _DEBUG
		template<typename T>
		void PrintVector4(const Math::Vector4<T>& v, const std::string& name, const Math::Vector4<T>& expected)
		{
			std::cout << name << ":" << v.ToString() << ""
				<< " | Expected: " << expected.ToString() << "\n";
		}

		template<typename T>
		void PrintScalar4(T value, const std::string& name, T expected)
		{
			std::cout << name << ": " << value << " | Expected: " << expected << "\n";
		}

		template<typename T>
		void TestVector4Type(const std::string& typeName)
		{
			std::cout << "Running Vector4 Tests for " << typeName << "\n";

			Math::Vector4<T> v1(static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(4));
			Math::Vector4<T> v2(static_cast<T>(5), static_cast<T>(6), static_cast<T>(7), static_cast<T>(8));

			PrintVector4(v1, "v1", Math::Vector4<T>(static_cast<T>(1), static_cast<T>(2), static_cast<T>(3), static_cast<T>(4)));
			PrintVector4(v2, "v2", Math::Vector4<T>(static_cast<T>(5), static_cast<T>(6), static_cast<T>(7), static_cast<T>(8)));

			Math::Vector4<T> sum = v1 + v2;
			PrintVector4(sum, "v1 + v2", Math::Vector4<T>(static_cast<T>(6), static_cast<T>(8), static_cast<T>(10), static_cast<T>(12)));

			Math::Vector4<T> diff = v1 - v2;
			PrintVector4(diff, "v1 - v2", Math::Vector4<T>(static_cast<T>(-4), static_cast<T>(-4), static_cast<T>(-4), static_cast<T>(-4)));

			Math::Vector4<T> neg = -v1;
			PrintVector4(neg, "-v1", Math::Vector4<T>(static_cast<T>(-1), static_cast<T>(-2), static_cast<T>(-3), static_cast<T>(-4)));

			Math::Vector4<T> scaled = v1 * static_cast<T>(2);
			PrintVector4(scaled, "v1 * 2", Math::Vector4<T>(static_cast<T>(2), static_cast<T>(4), static_cast<T>(6), static_cast<T>(8)));

			Math::Vector4<T> div = v1 / static_cast<T>(2);
			PrintVector4(div, "v1 / 2", Math::Vector4<T>(static_cast<T>(0.5), static_cast<T>(1), static_cast<T>(1.5), static_cast<T>(2)));

			T dot = Math::Vector4<T>::DotProduct(v1, v2);
			PrintScalar4(dot, "Dot Product (v1 . v2)", static_cast<T>(70));

			T magnitude = v1.Magnitude();
			PrintScalar4(magnitude, "Magnitude of v1", static_cast<T>(sqrt(static_cast<T>(30)))); // sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(30)

			Math::Vector4<T> normal = v1.Normal();
			T invMag = static_cast<T>(1) / magnitude;
			PrintVector4(normal, "Normal of v1", Math::Vector4<T>(static_cast<T>(1 * invMag), static_cast<T>(2 * invMag), static_cast<T>(3 * invMag), static_cast<T>(4 * invMag)));

			Math::Vector4<T> lerpResult = Math::Vector4<T>::Lerp(v1, v2, static_cast<T>(0.5));
			PrintVector4(lerpResult, "Lerp between v1 and v2 at t = 0.5", Math::Vector4<T>(static_cast<T>(3), static_cast<T>(4), static_cast<T>(5), static_cast<T>(6)));

			Math::Vector4<T> projectionResult = Math::Vector4<T>::ProjectVectorOnToVector(v1, v2);
			T projMag = static_cast<T>(sqrt(pow(v2.x, 2.0) + pow(v2.y, 2.0) + pow(v2.z, 2.0) + pow(v2.w, 2.0)));
			Math::Vector4<T> projNormal(v2.x / projMag, v2.y / projMag, v2.z / projMag, v2.w / projMag);
			T projScalar = Math::Vector4<T>::DotProduct(v1, projNormal);
			Math::Vector4<T> expectedResult = projNormal * projScalar;
			PrintVector4(projectionResult, "Projection of v1 onto v2", expectedResult);

			std::cout << "Vector4 Tests for " << typeName << " Completed\n";
		}

		inline void TestVector4()
		{
			TestVector4Type<float>("float");
			TestVector4Type<double>("double");
			TestVector4Type<int>("int");
		}
#endif // _DEBUG

	}
}

typedef Mjolnir::Math::Vector4<int> Vector4i;
typedef Mjolnir::Math::Vector4<float> Vector4f;
typedef Mjolnir::Math::Vector4<double> Vector4d;

#endif // MJOLNIR_VECTOR4_H
