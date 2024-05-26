// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_RAY_2D_H
#define MJOLNIR_RAY_2D_H

#include "Vector2.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Ray2D
		{
			Vector2<T> origin;
			Vector2<T> direction;
			
			virtual ~Ray2D() = default;

			Ray2D():
				origin(Vector2<T>::ZeroVector),
				direction(Vector2<T>::UpVector)
			{}

			Ray2D(const Ray2D& ray):
				origin(ray.origin),
				direction(ray.direction)
			{}

			Ray2D(Ray2D&& ray) noexcept
			{
				origin = ray.origin;
				direction = ray.direction;

				ray.origin = Vector2<T>();
				ray.direction = Vector2<T>();
			}

			Ray2D(const Vector2<T>& _origin, const Vector2<T>& _direction):
				origin(_origin),
				direction(_direction.Normal())
			{}

			Ray2D& operator=(const Ray2D& ray)
			{
				origin = ray.origin;
				direction = ray.direction;
				return *this;
			}

			Ray2D& operator=(Ray2D&& ray) noexcept
			{
				origin = ray.origin;
				direction = ray.direction;

				ray.origin = Vector2<T>();
				ray.direction = Vector2<T>();

				return *this;
			}

			bool operator==(const Ray2D& ray) const
			{
				return (origin == ray.origin) && (direction == ray.direction);
			}

			bool operator!=(const Ray2D& ray) const
			{
				return !(*this == _m);
			}

			Vector2<T> GetPoint(const T distance)
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
		void TestRay2D()
		{
			// Default Constructor
			Ray2D<float> ray1;
			assert(ray1.origin == Vector2<float>::ZeroVector);
			assert(ray1.direction == Vector2<float>::UpVector);
			std::cout << "Default constructor test passed\n";

			// Parameterized Constructor
			Vector2<float> origin(1.0f, 2.0f);
			Vector2<float> direction(3.0f, 4.0f);
			Ray2D<float> ray2(origin, direction);
			assert(ray2.origin == origin);
			assert(ray2.direction == direction.Normal());
			std::cout << "Parameterized constructor test passed\n";

			// Copy Constructor
			Ray2D<float> ray3(ray2);
			assert(ray3.origin == ray2.origin);
			assert(ray3.direction == ray2.direction);
			std::cout << "Copy constructor test passed\n";

			// Move Constructor
			Ray2D<float> ray4(std::move(ray2));
			assert(ray4.origin == origin);
			assert(ray4.direction == direction.Normal());
			// assert(ray2.origin == Vector2<float>());  // Ensure moved-from object is reset
			// assert(ray2.direction == Vector2<float>());
			std::cout << "Move constructor test passed\n";

			// Copy Assignment Operator
			Ray2D<float> ray5;
			ray5 = ray3;
			assert(ray5.origin == ray3.origin);
			assert(ray5.direction == ray3.direction);
			std::cout << "Copy assignment operator test passed\n";

			// Move Assignment Operator
			Ray2D<float> ray6;
			ray6 = std::move(ray3);
			assert(ray6.origin == origin);
			assert(ray6.direction == direction.Normal());
			// assert(ray3.origin == Vector2<float>());  // Ensure moved-from object is reset
			// assert(ray3.direction == Vector2<float>());
			std::cout << "Move assignment operator test passed\n";

			// GetPoint Method
			Vector2<float> point = ray4.GetPoint(2.0f);
			Vector2<float> expectedPoint = origin + (2.0f * direction.Normal());
			assert(point == expectedPoint);
			std::cout << "GetPoint method test passed\n";

			// ToString Method
			std::string str = ray4.ToString();
			std::string expectedStr = "(Origin:(1, 2), Direction:(0.6, 0.8))"; // Assuming direction is normalized to (0.6, 0.8)
			assert(str == expectedStr);
			std::cout << "ToString method test passed\n";
		}
#endif
	}
}

typedef Mjolnir::Math::Ray2D<float> Ray2Df;
typedef Mjolnir::Math::Ray2D<double> Ray2Dd;

#endif // MJOLNIR_RAY_2D_H