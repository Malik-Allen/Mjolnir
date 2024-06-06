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
			
			explicit Circle(T _radius):
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

			Circle& operator=(Circle& circle) noexcept
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
			Vector2<T> centerToPoint = circle.center - point;

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
	}
}

typedef Mjolnir::Math::Circle<int> CircleI;
typedef Mjolnir::Math::Circle<float> CircleF;
typedef Mjolnir::Math::Circle<double> CircleD;

#endif // MJOLNIR_CIRCLE_H