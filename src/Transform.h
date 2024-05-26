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
				m_translation(Vector3f(0.0f)),
				m_orientation(Quatf::Identity),
				m_scale(Vector3f(1.0f)),
				m_modelMatrix(Matrix4f::Identity)
			{}

			Transform(const Vector3f& position, const Quatf& rotation, const Vector3f& scale):
				m_translation(position),
				m_orientation(rotation),
				m_scale(scale),
				m_modelMatrix(Matrix4f::Identity)
			{}

			Transform(const Transform& transform) = default;
			Transform(Transform&& transform) = default;
			Transform& operator=(const Transform& transform) = default;
			Transform& operator=(Transform&& transform) = default;

			void SetPosition(const Vector3f& position);
			void SetRotation(const Vector3f& axis, const float angle);
			void SetRotation(const Quatf& rotation);
			void SetScale(const Vector3f& scale);

			Vector3f GetPosition() const;
			Vector3f GetScale() const;
			Quatf GetRotation() const;

			std::string ToString() const;

		protected:
			Vector3f m_translation;
			Quatf m_orientation;
			Vector3f m_scale;
			Matrix4f m_modelMatrix;

		private:
			void Update();
		};

#if _DEBUG
		void TestTransform();
#endif // _DEBUG
	}
}

#endif // MJOLNIR_TRANSFORM_H