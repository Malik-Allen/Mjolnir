// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_RAY_3D_H
#define MJOLNIR_RAY_3D_H

#include "Vector3.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Ray3D
		{
			Vector3<T> origin;
			Vector3<T> direction;

			virtual ~Ray3D() = default;

			Ray3D():
				origin(Vector3<T>::ZeroVector),
				direction(Vector3<T>::UpVector)
			{}

			Ray3D(const Ray3D& ray):
				origin(ray.origin),
				direction(ray.direction)
			{}

			Ray3D(Ray3D&& ray) noexcept
			{
				origin = ray.origin;
				direction = ray.direction;

				ray.origin = Vector3<T>();
				ray.direction = Vector3<T>();
			}

			Ray3D(const Vector3<T>& _origin, const Vector3<T>& _direction):
				origin(_origin),
				direction(_direction.Normal())
			{}

			Ray3D& operator=(const Ray3D& ray)
			{
				origin = ray.origin;
				direction = ray.direction;
				return *this;
			}

			Ray3D& operator=(Ray3D&& ray) noexcept
			{
				origin = ray.origin;
				direction = ray.direction;

				ray.origin = Vector3<T>();
				ray.direction = Vector3<T>();

				return *this;
			}

			bool operator==(const Ray3D& ray) const
			{
				return (origin == ray.origin) && (direction == ray.direction);
			}

			bool operator!=(const Ray3D& ray) const
			{
				return !(*this == ray);
			}

			Vector3<T> GetPoint(const T distance)
			{
				return origin + (distance * direction);
			}

			std::string ToString() const
			{
				std::ostringstream oss;
				oss << "(Origin:" << origin.ToString() << ", Direction:" << direction.ToString() << ")";
				return oss.str();
			}
		};


#if _DEBUG
		void TestRay3D()
		{
			std::cout << "Testing Ray3D\n";
			// Default Constructor
			Ray3D<float> ray1;
			assert(ray1.origin == Vector3<float>::ZeroVector);
			assert(ray1.direction == Vector3<float>::UpVector);
			std::cout << "Default constructor test passed\n";

			// Parameterized Constructor
			Vector3<float> origin(1.0f, 2.0f, 3.0f);
			Vector3<float> direction(4.0f, 5.0f, 6.0f);
			Ray3D<float> ray2(origin, direction);
			assert(ray2.origin == origin);
			assert(ray2.direction == direction.Normal());
			std::cout << "Parameterized constructor test passed\n";

			// Copy Constructor
			Ray3D<float> ray3(ray2);
			assert(ray3.origin == ray2.origin);
			assert(ray3.direction == ray2.direction);
			std::cout << "Copy constructor test passed\n";

			// Move Constructor
			Ray3D<float> ray4(std::move(ray2));
			assert(ray4.origin == origin);
			assert(ray4.direction == direction.Normal());
			// assert(ray2.origin == Vector3<float>());  // Ensure moved-from object is reset
			// assert(ray2.direction == Vector3<float>());
			std::cout << "Move constructor test passed\n";

			// Copy Assignment Operator
			Ray3D<float> ray5;
			ray5 = ray3;
			assert(ray5.origin == ray3.origin);
			assert(ray5.direction == ray3.direction);
			std::cout << "Copy assignment operator test passed\n";

			// Move Assignment Operator
			Ray3D<float> ray6;
			ray6 = std::move(ray3);
			assert(ray6.origin == origin);
			assert(ray6.direction == direction.Normal());
			// assert(ray3.origin == Vector3<float>());  // Ensure moved-from object is reset
			// assert(ray3.direction == Vector3<float>());
			std::cout << "Move assignment operator test passed\n";

			// Comparison Operators
			Ray3D<float> ray7(origin, direction);
			assert(ray7 == ray6);
			assert(!(ray7 != ray6));
			Ray3D<float> ray8(Vector3<float>(7.0f, 8.0f, 9.0f), direction);
			assert(ray8 != ray7);
			std::cout << "Comparison operators test passed\n";

			// GetPoint Method
			Vector3<float> point = ray4.GetPoint(2.0f);
			Vector3<float> expectedPoint = origin + (2.0f * direction.Normal());
			assert(point == expectedPoint);
			std::cout << "GetPoint method test passed\n";

			// ToString Method
			std::string str = ray4.ToString();
			std::string expectedStr = "(Origin:(1, 2, 3), Direction:(0.455842, 0.569803, 0.683763))"; // Assuming direction is normalized to these values
			assert(str == expectedStr);
			std::cout << "ToString method test passed\n\n";
		}
#endif
	}
}

typedef Mjolnir::Math::Ray3D<float> Ray3Df;
typedef Mjolnir::Math::Ray3D<double> Ray3Dd;

#endif // MJOLNIR_RAY_3D_H