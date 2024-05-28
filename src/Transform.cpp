// MIT License, Copyright (c) 2024 Malik Allen

#include "Transform.h"

namespace Mjolnir
{
	namespace Math
	{
		void Transform::SetPosition(const Vector3f& position)
		{
			translation = position;
			Update();
		}

		void Transform::SetRotation(const Vector3f & axis, const float angle)
		{
			orientation.RotateAxis(axis, angle);
			Update();
		}

		void Transform::SetRotation(const Quatf& rotation)
		{
			orientation = rotation;
			Update();
		}

		void Transform::SetScale(const Vector3f & _scale)
		{
			scale = _scale;
			Update();
		}

		Vector3f Transform::GetPosition() const
		{
			return translation;
		}

		Vector3f Transform::GetScale() const
		{
			return scale;
		}

		Quatf Transform::GetRotation() const
		{
			return orientation;
		}

		std::string Transform::ToString() const
		{
			std::ostringstream oss;
			oss << "[p:" << translation.ToString() << ", r:" << orientation.ToString() << ", s:" << scale.ToString() << "]" << "\n" << modelMatrix.ToString();
			return oss.str();
		}

		void Transform::Update()
		{
			// In order to obtain a uniform model matrix, we multiply our matrices in this order: Scale -> Rotate -> Translate
			// In matrix multiplication the order is read from right to left L <-- R
			modelMatrix = Matrix4f::Translate(translation) * orientation.RotationMatrix4() * Matrix4f::Scale(scale);
		}

#if _DEBUG
		void TestTransform()
		{
			Math::Transform transform;
			std::cout <<  "init:" << transform.ToString() << "\n";

			transform.SetPosition(Vector3f(-51.0f, -1124.0f, 86887.0f));
			std::cout <<  "set position:" << transform.ToString() << "\n";

			transform.SetScale(Vector3f(1.0f, 20.0f, 1001.0f));
			std::cout <<  "set scale:" << transform.ToString() << "\n";

			transform.SetRotation(Vector3f(0.0f, 1.0f, 0.0f), 45.0f);
			std::cout <<  "set rotation:" << transform.ToString() << "\n";

			Quatf orientation(transform.GetRotation());

			Quatf rotation(Vector3f(0.0f, 1.0f, 0.0f), -90.0f);
			transform.SetRotation((orientation * rotation));

			std::cout <<  "acummulate rotation:" << transform.ToString() << "\n";
		}
#endif // _DEBUG
	}
}