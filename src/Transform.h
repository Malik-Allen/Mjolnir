// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_TRANSFORM_H
#define MJOLNIR_TRANSFORM_H

#include "Quaternion.h"

namespace Mjolnir
{
	namespace Math
	{
		struct Transform
		{
			virtual ~Transform() = default;

			Transform():
				translation(Vector3f(0.0f)),
				orientation(Quatf::Identity),
				scale(Vector3f(1.0f)),
				modelMatrix(Matrix4f::Identity)
			{}

			Transform(const Vector3f& position, const Quatf& rotation, const Vector3f& _scale):
				translation(position),
				orientation(rotation),
				scale(_scale),
				modelMatrix(Matrix4f::Identity)
			{}

			Transform(const Transform& transform)
			{
				translation = transform.translation;
				orientation = transform.orientation;
				scale = transform.scale;
				modelMatrix = transform.modelMatrix;
			}

			Transform(Transform&& transform) noexcept
			{
				translation = transform.translation;
				orientation = transform.orientation;
				scale = transform.scale;
				modelMatrix = transform.modelMatrix;

				transform.Reset();
			}
			
			Transform& operator=(const Transform& transform) 
			{
				translation = transform.translation;
				orientation = transform.orientation;
				scale = transform.scale;
				modelMatrix = transform.modelMatrix;
				return *this;
			}

			Transform& operator=(Transform&& transform) noexcept
			{
				translation = transform.translation;
				orientation = transform.orientation;
				scale = transform.scale;
				modelMatrix = transform.modelMatrix;

				transform.Reset();
				return *this;
			}

			void SetPosition(const Vector3f& position);
			void Rotate(const Vector3f& axis, const float angle);
			void SetRotation(const Quatf& rotation);
			void SetScale(const Vector3f& scale);

			Vector3f GetPosition() const;
			Vector3f GetScale() const;
			Quatf GetRotation() const;
			Matrix4f GetModelMatrix() const;

			std::string ToString() const;

		protected:
			Vector3f translation;
			Quatf orientation;
			Vector3f scale;
			Matrix4f modelMatrix;

		private:
			void Update();
			void Reset();
		};

#if _DEBUG
		void TestTransform();
#endif // _DEBUG
	}
}

typedef Mjolnir::Math::Transform Transform;

#endif // MJOLNIR_TRANSFORM_H