// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_VECTOR2_H
#define MJOLNIR_VECTOR2_H

#include "MjolnirConstants.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Vector2
		{
			T x, y;

			~Vector2() = default;

			explicit Vector2(T s = 0):
				x(s),
				y(s)
			{}

			Vector2(T _x, T _y):
				x(_x),
				y(_y)
			{}

			Vector2(const Vector2<T>& v):
				x(v.x),
				y(v.y)
			{}

			Vector2(Vector2<T>&& v) noexcept:
				x(v.x),
				y(v.y)
			{
				v.x = v.y = T();
			}

			Vector2& operator=(const Vector2<T>& v)
			{
				x = v.x;
				y = v.y;
				return *this;
			}

			Vector2& operator=(Vector2<T>&& v) noexcept
			{
				if(this != &v)
				{
					x = v.x;
					y = v.y;
					v.x = v.y = T();
				}
				return *this;
			}

			bool operator==(const Vector2<T>& v) const
			{
				return (x == v.x) && (y == v.y);
			}

			bool operator!=(const Vector2<T>& v) const
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

			Vector2 operator-() const
			{
				return Vector2(-x, -y);
			}

			Vector2 operator+(const Vector2<T>& v) const
			{
				return Vector2<T>(
					x + v.x,
					y + v.y
				);
			}

			Vector2 operator-(const Vector2<T>& v) const
			{
				return Vector2<T>(
					x - v.x,
					y - v.y
				);
			}

			Vector2 operator+=(const Vector2<T>& v)
			{
				x += v.x;
				y += v.y;
				return *this;
			}

			Vector2 operator-=(const Vector2<T>& v)
			{
				x -= v.x;
				y -= v.y;
				return *this;
			}

			Vector2 operator*(const T& s) const
			{
				return Vector2(
					x * s,
					y * s
				);
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * vector = product;
			static friend Vector2 operator*(const T& s, const Vector2& v)
			{
				return v * s;
			}

			Vector2 operator*=(const T& s)
			{
				x *= s;
				y *= s;
				return *this;
			}

			Vector2 operator/(const T& s) const
			{
#ifdef _DEBUG
				if(std::fabs(s) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				return Vector2(
					x / s,
					y / s
				);
			}

			Vector2 operator/=(const T& s)
			{
#ifdef _DEBUG
				if(std::fabs(s) < NEARLY_ZERO)
				{
					throw std::runtime_error("Error! Dividing by nearly zero!");
				}
#endif
				x /= s;
				y /= s;
				return *this;
			}

			T Magnitude() const;
			T MagnitudeSquared() const;
			Vector2<T> Normal() const;
			void Normalize();
			std::string ToString() const;

			static T Distance(const Vector2& v1, const Vector2& v2);
			static T DotProduct(const Vector2& v1, const Vector2& v2);
			static T Angle(const Vector2& v1, const Vector2& v2);
			static Vector2 Lerp(const Vector2& v1, const Vector2& v2, const T t);
			static Vector2 ProjectVectorOnToVector(const Vector2& project, const Vector2& target);

			static const Vector2 ZeroVector;
			static const Vector2 UpVector;
			static const Vector2 RightVector;
			static const Vector2 DownVector;
			static const Vector2 LeftVector;
		};

		template <typename T>
		T Vector2<T>::Magnitude() const
		{
			return static_cast<T>(sqrt((x * x) + (y * y)));
		}

		template <typename T>
		T Vector2<T>::MagnitudeSquared() const
		{
			return (x * x) + (y * y);
		}

		template <typename T>
		Vector2<T> Vector2<T>::Normal() const
		{
			T a = static_cast<T>(sqrt((x * x) + (y * y)));
#ifdef _DEBUG
			if(fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			return Vector2(x / a, y / a);
		}

		template <typename T>
		void Vector2<T>::Normalize()
		{
			T a = static_cast<T>(sqrt((x * x) + (y * y)));
#ifdef _DEBUG
			if(fabs(a) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			x /= a;
			y /= a;
		}

		template<typename T>
		std::string Vector2<T>::ToString() const
		{
			std::ostringstream oss;
			oss << "(" << x << ", " << y << ")";
			return oss.str();
		}

		template <typename T>
		T Vector2<T>::Distance(const Vector2& v1, const Vector2& v2)
		{
			const Vector2 diff(v1 - v2);
			return diff.Magnitude();
		}

		template <typename T>
		T Vector2<T>::DotProduct(const Vector2& v1, const Vector2& v2)
		{
			return static_cast<T>(v1.x * v2.x + v1.y * v2.y);
		}

		template <typename T>
		T Vector2<T>::Angle(const Vector2& v1, const Vector2& v2)
		{
			T a = DotProduct(v1, v2);
			T b = static_cast<T>(sqrt(v1.Magnitude() * v2.Magnitude()));	// Multiplying the magnitude of each vector
			return  static_cast<T>(acos(a / b));	// Taking the cos inverse of the dot product / magnitude result to get the angle in radians
		}

		template <typename T>
		Vector2<T> Vector2<T>::Lerp(const Vector2& v1, const Vector2& v2, const T t)
		{
			return Vector2(v1 + t * (v2 - v1));
		}

		template <typename T>
		Vector2<T> Vector2<T>::ProjectVectorOnToVector(const Vector2& project, const Vector2& target)
		{
			T mag = target.Magnitude();
#ifdef _DEBUG
			if(fabs(mag) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero in vector projection!"); // Checking the magnitude of the target vector
			}
#endif
			Vector2 normal(target.x / mag, target.y / mag);     // Normal of the target vector
			T scalar = DotProduct(project, normal);     // Dot product of the vector to project and the normalized target vector
			return normal * scalar;     // Returning the normal scaled to the length of the projected vector's shadow
		}

		template<typename T>
		const Vector2<T> Vector2<T>::ZeroVector = Vector2<T>(0, 0);

		template<typename T>
		const Vector2<T> Vector2<T>::UpVector = Vector2<T>(0, 1);

		template<typename T>
		const Vector2<T> Vector2<T>::RightVector = Vector2<T>(1, 0);

		template<typename T>
		const Vector2<T> Vector2<T>::DownVector = Vector2<T>(0, -1);

		template<typename T>
		const Vector2<T> Vector2<T>::LeftVector = Vector2<T>(-1, 0);

#if _DEBUG
		template<typename T>
		void PrintVector2(const Math::Vector2<T>& v, const std::string& name, const Math::Vector2<T>& expected)
		{
			std::cout << name << ":" << v.ToString() << ""
				<< " | Expected: " << expected.ToString() << "\n";
		}

		template<typename T>
		void PrintScalar2(T value, const std::string& name, T expected)
		{
			std::cout << name << ": " << value << " | Expected: " << expected << "\n";
		}

		template<typename T>
		void TestVector2Type(const std::string& typeName)
		{
			std::cout << "Running Vector2 Tests for " << typeName << "\n";

			Math::Vector2<T> v1(static_cast<T>(1), static_cast<T>(2));
			Math::Vector2<T> v2(static_cast<T>(3), static_cast<T>(4));

			PrintVector2(v1, "v1", Math::Vector2<T>(static_cast<T>(1), static_cast<T>(2)));
			PrintVector2(v2, "v2", Math::Vector2<T>(static_cast<T>(3), static_cast<T>(4)));

			Math::Vector2<T> sum = v1 + v2;
			PrintVector2(sum, "v1 + v2", Math::Vector2<T>(static_cast<T>(4), static_cast<T>(6)));

			Math::Vector2<T> diff = v1 - v2;
			PrintVector2(diff, "v1 - v2", Math::Vector2<T>(static_cast<T>(-2), static_cast<T>(-2)));

			Math::Vector2<T> neg = -v1;
			PrintVector2(neg, "-v1", Math::Vector2<T>(static_cast<T>(-1), static_cast<T>(-2)));

			Math::Vector2<T> scaled = v1 * static_cast<T>(2);
			PrintVector2(scaled, "v1 * 2", Math::Vector2<T>(static_cast<T>(2), static_cast<T>(4)));

			Math::Vector2<T> div = v1 / static_cast<T>(2);
			PrintVector2(div, "v1 / 2", Math::Vector2<T>(static_cast<T>(0.5), static_cast<T>(1)));

			T dot = Math::Vector2<T>::DotProduct(v1, v2);
			PrintScalar2(dot, "Dot Product (v1 . v2)", static_cast<T>(11));

			T magnitude = v1.Magnitude();
			PrintScalar2(magnitude, "Magnitude of v1", static_cast<T>(sqrt(static_cast<T>(5)))); // sqrt(1^2 + 2^2) = sqrt(5)

			Math::Vector2<T> normal = v1.Normal();
			T invMag = static_cast<T>(1) / magnitude;
			PrintVector2(normal, "Normal of v1", Math::Vector2<T>(static_cast<T>(1 * invMag), static_cast<T>(2 * invMag)));

			Math::Vector2<T> lerpResult = Math::Vector2<T>::Lerp(v1, v2, static_cast<T>(0.5));
			PrintVector2(lerpResult, "Lerp between v1 and v2 at t = 0.5", Math::Vector2<T>(static_cast<T>(2), static_cast<T>(3)));

			Math::Vector2<T> projectionResult = Math::Vector2<T>::ProjectVectorOnToVector(v1, v2);
			T projMag = static_cast<T>(sqrt(pow(v2.x, 2.0f) + pow(v2.y, 2.0f)));
			Math::Vector2<T> projNormal(v2.x / projMag, v2.y / projMag);
			T projScalar = Math::Vector2<T>::DotProduct(v1, projNormal);
			Math::Vector2<T> expectedResult = projNormal * projScalar;
			PrintVector2(projectionResult, "Projection of v1 onto v2", expectedResult);

			std::cout << "Vector2 Tests for " << typeName << " Completed\n";
		}

		inline void TestVector2()
		{
			TestVector2Type<float>("float");
			TestVector2Type<double>("double");
			TestVector2Type<int>("int");
		}
#endif // _DEBUG
	}
}

typedef Mjolnir::Math::Vector2<int> Vector2i;
typedef Mjolnir::Math::Vector2<float> Vector2f;
typedef Mjolnir::Math::Vector2<double> Vector2d;

#endif // MJOLNIR_VECTOR2_H
