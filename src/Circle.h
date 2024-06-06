// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_CIRCLE_H
#define MJOLNIR_CIRCLE_H

#include "Vector2.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Circle
		{
			Vector2<T> center;
			T radius;

			virtual ~Circle() = default;
			
			explicit Circle(T _radius = 1):
				center(Vector2<T>::ZeroVector),
				radius(_radius)
			{}

			Circle(Vector2<T> _center, T _radius):
				center(_center),
				radius(_radius)
			{}

			Circle(const Circle& circle):
				center(circle.center),
				radius(circle.radius)
			{}

			Circle(Circle&& circle) noexcept
			{
				center = circle.center;
				radius = circle.radius;

				circle.center = Vector2<T>();
				circle.radius = T();
			}

			Circle& operator=(const Circle& circle)
			{
				center = circle.center;
				radius = circle.radius;
				return *this;
			}

			Circle& operator=(Circle&& circle) noexcept
			{
				center = circle.center;
				radius = circle.radius;

				circle.center = Vector2<T>();
				circle.radius = T();
				return *this;
			}

			bool operator==(const Circle& circle) const
			{
				return (center == circle.center) && (radius == circle.radius);
			}

			bool operator!=(const Circle& circle) const
			{
				return !(*this == circle);
			}

			std::string ToString() const
			{
				std::ostringstream oss;
				oss << "(Center:" << center.ToString() << ", Radius:" << radius << ")";
				return oss.str();
			}

			enum ProjectionMode
			{
				KeepInside = 0,
				ProjectOntoSurface
			};

			static Vector2<T> ClosestPoint(const Vector2<T>& point, const Circle<T>& circle, const ProjectionMode& projectionMode);
		};

		template<typename T>
		inline Vector2<T> Circle<T>::ClosestPoint(const Vector2<T>& point, const Circle<T>& circle, const ProjectionMode& projectionMode)
		{
			Vector2<T> centerToPoint = point - circle.center;

			const T radius = circle.radius;
			const T distanceSquared = centerToPoint.MagnitudeSquared();

			if(distanceSquared < radius * radius)
			{
				if(projectionMode == ProjectionMode::KeepInside)
				{
					return point;
				}
				else if(projectionMode == ProjectionMode::ProjectOntoSurface)
				{
					if(std::fabs(distanceSquared) < NEARLY_ZERO)
					{
						centerToPoint = Vector2<T>::UpVector;	// Some arbitrary direction, in the case when we are passed a point that lies on the center of the circle.
					}
				}
			}
			
			return circle.center + (centerToPoint.Normal() * radius);
		}

#if _DEBUG
		void TestCircle() {
			std::cout << "\nPerforming Circle Unit Tests\n";

			// Default Constructor
			Circle<float> circle1(5.0f);
			assert(circle1.center == Vector2<float>::ZeroVector);
			assert(circle1.radius == 5.0f);
			std::cout << "Default constructor test passed\n";

			// Parameterized Constructor
			Vector2<float> center(1.0f, 2.0f);
			Circle<float> circle2(center, 10.0f);
			assert(circle2.center == center);
			assert(circle2.radius == 10.0f);
			std::cout << "Parameterized constructor test passed\n";

			// Copy Constructor
			Circle<float> circle3(circle2);
			assert(circle3.center == circle2.center);
			assert(circle3.radius == circle2.radius);
			std::cout << "Copy constructor test passed\n";

			// Move Constructor
			Circle<float> circle4(std::move(circle2));
			assert(circle4.center == center);
			assert(circle4.radius == 10.0f);
			assert(circle2.center == Vector2<float>());  // Ensure moved-from object is reset
			assert(circle2.radius == 0.0f);
			std::cout << "Move constructor test passed\n";

			// Copy Assignment Operator
			Circle<float> circle5;
			circle5 = circle3;
			assert(circle5.center == circle3.center);
			assert(circle5.radius == circle3.radius);
			std::cout << "Copy assignment operator test passed\n";

			// Move Assignment Operator
			Circle<float> circle6;
			circle6 = std::move(circle3);
			assert(circle6.center == center);
			assert(circle6.radius == 10.0f);
			assert(circle3.center == Vector2<float>());  // Ensure moved-from object is reset
			assert(circle3.radius == 0.0f);
			std::cout << "Move assignment operator test passed\n";

			// Comparison Operators
			Circle<float> circle7(center, 10.0f);
			assert(circle7 == circle6);
			assert(!(circle7 != circle6));
			Circle<float> circle8(Vector2<float>(7.0f, 8.0f), 15.0f);
			assert(circle8 != circle7);
			std::cout << "Comparison operators test passed\n";

			// ToString Method
			std::string str = circle6.ToString();
			std::string expectedStr = "(Center:(1, 2), Radius:10)";
			assert(str == expectedStr);
			std::cout << "ToString method test passed\n";

			// ClosestPoint Method
			Vector2<float> point(12.0f, 2.0f);
			Vector2<float> closestPoint = Circle<float>::ClosestPoint(point, circle6, Circle<float>::ProjectionMode::ProjectOntoSurface);
			Vector2<float> expectedClosestPoint(11.0f, 2.0f);  // Calculate manually based on the circle properties
			assert(closestPoint == expectedClosestPoint);
			std::cout << "ClosestPoint method test passed\n";
		}
#endif // _DEBUG
	}
}

typedef Mjolnir::Math::Circle<int> CircleI;
typedef Mjolnir::Math::Circle<float> CircleF;
typedef Mjolnir::Math::Circle<double> CircleD;

#endif // MJOLNIR_CIRCLE_H