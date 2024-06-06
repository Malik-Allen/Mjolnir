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

			explicit Sphere(T _radius):
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

			Sphere& operator=(Sphere& sphere) noexcept
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
			Vector3<T> centerToPoint = sphere.center - point;

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
	}
}

typedef Mjolnir::Math::Sphere<float> SphereF;
typedef Mjolnir::Math::Sphere<double> SphereD;

#endif // MJOLNIR_SPHERE_H