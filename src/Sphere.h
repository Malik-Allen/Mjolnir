// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_SPHERE_H
#define MJOLNIR_SPHERE_H

#include "Vector3.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Sphere
		{
			Vector3<T> center;
			T radius;

			virtual ~Sphere() = default;

			explicit Sphere(T _radius = 1):
				center(Vector3<T>::ZeroVector),
				radius(_radius)
			{}

			Sphere(Vector3<T> _center, T _radius):
				center(_center),
				radius(_radius)
			{}

			Sphere(const Sphere& sphere):
				center(sphere.center),
				radius(sphere.radius)
			{}

			Sphere(Sphere&& sphere) noexcept
			{
				center = sphere.center;
				radius = sphere.radius;

				sphere.center = Vector3<T>();
				sphere.radius = T();
			}

			Sphere& operator=(const Sphere& sphere)
			{
				center = sphere.center;
				radius = sphere.radius;
				return *this;
			}

			Sphere& operator=(Sphere&& sphere) noexcept
			{
				center = sphere.center;
				radius = sphere.radius;

				sphere.center = Vector3<T>();
				sphere.radius = T();
				return *this;
			}

			bool operator==(const Sphere& sphere) const
			{
				return (center == sphere.center) && (radius == sphere.radius);
			}

			bool operator!=(const Sphere& sphere) const
			{
				return !(*this == sphere);
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

			static Vector3<T> ClosestPoint(const Vector3<T>& point, const Sphere<T>& sphere, const ProjectionMode& projectionMode);
		};

		template<typename T>
		inline Vector3<T> Sphere<T>::ClosestPoint(const Vector3<T>& point, const Sphere<T>& sphere, const ProjectionMode& projectionMode)
		{
			Vector3<T> centerToPoint = point - sphere.center;

			const T radius = sphere.radius;
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
						centerToPoint = Vector3<T>::UpVector;	// Some arbitrary direction, in the case when we are passed a point that lies on the center of the sphere.
					}
				}
			}

			return sphere.center + (centerToPoint.Normal() * radius);
		}

#if _DEBUG
		inline void TestSphere() {
			std::cout << "\nPerforming Sphere Unit Tests\n";
			// Default Constructor
			Sphere<float> sphere1(5.0f);
			assert(sphere1.center == Vector3<float>::ZeroVector);
			assert(sphere1.radius == 5.0f);
			std::cout << "Default constructor test passed\n";

			// Parameterized Constructor
			Vector3<float> center(1.0f, 2.0f, 3.0f);
			Sphere<float> sphere2(center, 10.0f);
			assert(sphere2.center == center);
			assert(sphere2.radius == 10.0f);
			std::cout << "Parameterized constructor test passed\n";

			// Copy Constructor
			Sphere<float> sphere3(sphere2);
			assert(sphere3.center == sphere2.center);
			assert(sphere3.radius == sphere2.radius);
			std::cout << "Copy constructor test passed\n";

			// Move Constructor
			Sphere<float> sphere4(std::move(sphere2));
			assert(sphere4.center == center);
			assert(sphere4.radius == 10.0f);
			assert(sphere2.center == Vector3<float>());  // Ensure moved-from object is reset
			assert(sphere2.radius == 0.0f);
			std::cout << "Move constructor test passed\n";

			// Copy Assignment Operator
			Sphere<float> sphere5;
			sphere5 = sphere3;
			assert(sphere5.center == sphere3.center);
			assert(sphere5.radius == sphere3.radius);
			std::cout << "Copy assignment operator test passed\n";

			// Move Assignment Operator
			Sphere<float> sphere6;
			sphere6 = std::move(sphere3);
			assert(sphere6.center == center);
			assert(sphere6.radius == 10.0f);
			assert(sphere3.center == Vector3<float>());  // Ensure moved-from object is reset
			assert(sphere3.radius == 0.0f);
			std::cout << "Move assignment operator test passed\n";

			// Comparison Operators
			Sphere<float> sphere7(center, 10.0f);
			assert(sphere7 == sphere6);
			assert(!(sphere7 != sphere6));
			Sphere<float> sphere8(Vector3<float>(7.0f, 8.0f, 9.0f), 15.0f);
			assert(sphere8 != sphere7);
			std::cout << "Comparison operators test passed\n";

			// ToString Method
			std::string str = sphere6.ToString();
			std::string expectedStr = "(Center:(1, 2, 3), Radius:10)";
			assert(str == expectedStr);
			std::cout << "ToString method test passed\n";

			// ClosestPoint Method
			Vector3<float> point(12.0f, 2.0f, 3.0f);
			Vector3<float> closestPoint = Sphere<float>::ClosestPoint(point, sphere6, Sphere<float>::ProjectionMode::ProjectOntoSurface);
			Vector3<float> expectedClosestPoint(11.0f, 2.0f, 3.0f);  // Calculate manually based on the sphere properties
			assert(closestPoint == expectedClosestPoint);
			std::cout << "ClosestPoint method test passed\n";
		}
#endif // _DEBUG
	}
}

typedef Mjolnir::Math::Sphere<float> SphereF;
typedef Mjolnir::Math::Sphere<double> SphereD;

#endif // MJOLNIR_SPHERE_H