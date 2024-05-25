// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_QUATERNION_H
#define MJOLNIR_QUATERNION_H

#include "Matrix3.h"
#include "Matrix4.h"

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Quaternion
		{
			// Quaternion
			// a + b*i + c*j + d*k
			// In this representation:
			// a: Scalar part (real(r) component)
			// b, c, d: Vector parts (imaginary components)
			// i, j, k: Unit vectors in quaternion space

			T i, j, k, r;

			~Quaternion() = default;

			explicit Quaternion(T s = 0.0)
			{
				if(s == 0.0)
				{
					i = s;
					j = s;
					k = s;
					r = 1.0;
				}
				else
				{
					i = s;
					j = s;
					k = s;
					r = s;
				}
			}

			Quaternion(T _x, T _y, T _z, T _w):
				i(_x),
				j(_y),
				k(_z),
				r(_w)
			{}

			Quaternion(const Quaternion& q):
				i(q.i),
				j(q.j),
				k(q.k),
				r(q.r)
			{}

			Quaternion(Quaternion&& q) noexcept:
				i(q.i),
				j(q.j),
				k(q.k),
				r(q.r)
			{
				q.i = q.j = q.k = q.r = T();
			}

			Quaternion(const Vector3<T>& axis, const T angle)
			{
				i = 0.0f;
				j = 0.0f;
				k = 0.0f;
				r = 1.0f;

				RotateAxis(axis, angle);
			}

			Quaternion& operator=(const Quaternion& q)
			{
				i = q.i;
				j = q.j;
				k = q.k;
				r = q.r;
				return *this;
			}

			Quaternion& operator=(Quaternion&& q) noexcept
			{
				if(this != &q)
				{
					i = q.i;
					j = q.j;
					k = q.k;
					r = q.r;
					q.i = q.j = q.k = q.r = T();
				}
				return *this;
			}

			T operator[](const unsigned int index) const
			{
				return *(&i + index);
			}

			T& operator[](const unsigned int index)
			{
				return *(&i + index);
			}

			Quaternion operator-() const
			{
				return Quaternion(-i, -j, -k, -r);
			}

			Quaternion operator+(const Quaternion& q) const
			{
				return Quaternion(
					i + q.i,
					j + q.j,
					k + q.k,
					r + q.r
				);
			}

			Quaternion operator-(const Quaternion& q) const
			{
				return Quaternion(
					i - q.i,
					j - q.j,
					k - q.k,
					r - q.r
				);
			}

			Quaternion operator+=(const Quaternion& q)
			{
				i += q.i;
				j += q.j;
				k += q.k;
				r += q.r;
				return *this;
			}

			Quaternion operator-=(const Quaternion& q)
			{
				i -= q.i;
				j -= q.j;
				k -= q.k;
				r -= q.r;
				return *this;
			}

			Quaternion operator*(const Quaternion& q) const
			{
				// Derived from: https://en.wikipedia.org/wiki/Quaternion#Algebraic_properties
				// *this Quaternion represents p in the following equation

				const T pScalarPart = this->r;
				const T qScalarPart = q.r;

				Vector3<T> pVectorPart(this->i, this->j, this->k);
				Vector3<T> qVectorPart(q.i, q.j, q.k);

				const T resultScalarPart = (pScalarPart * qScalarPart) - Vector3<T>::DotProduct(pVectorPart, qVectorPart);

				const Vector3<T> resultVectorPart = (pScalarPart * qVectorPart) + (qScalarPart * pVectorPart) + (Vector3<T>::CrossProduct(pVectorPart, qVectorPart));

				return Quaternion(
					resultVectorPart.x,
					resultVectorPart.y,
					resultVectorPart.z,
					resultScalarPart
				);
			}

			Quaternion operator*(const T& s) const
			{
				return Quaternion(
					i * s,
					j * s,
					k * s,
					r * s
				);
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * quaternion = product;
			static friend Quaternion operator*(const T& s, const Quaternion& q)
			{
				return q * s;
			}

			Quaternion operator*=(const T& s)
			{
				i *= s;
				j *= s;
				k *= s;
				r *= s;
				return *this;
			}

			Quaternion<T> Inverse() const;
			Quaternion<T> Conjugate() const;
			Quaternion<T> Normal() const;
			void Normalize();
			std::string ToString() const;
			Vector3<T> RotateVector(const Vector3<T>& v) const;
			Matrix3<T> RotationMatrix3() const;
			Matrix4<T> RotationMatrix4() const;
			void RotateAxis(const Vector3<T>& axis, const T& angle);

			static Quaternion<T> Inverse(const Quaternion<T>& q);
			static Quaternion<T> Conjugate(const Quaternion<T>& q);
			static Quaternion<T> Normalize(const Quaternion<T>& q);
			static T AngularDistance(const Quaternion& q1, const Quaternion& q2);
			static Quaternion Lerp(const Quaternion& q1, const Quaternion& q2, const T t);
			static Quaternion SLerp(const Quaternion& q1, const Quaternion& q2, const T t);
			

			static const Quaternion<T> Identity;
		};

		template <typename T>
		Quaternion<T> Quaternion<T>::Inverse() const
		{
			// Derived-from: https://en.wikipedia.org/wiki/Quaternion#Algebraic_properties
			Quaternion conj(-i, -j, -k, r);
			const T magSqrd = pow(conj.i, 2.0f) + pow(conj.j, 2.0f) + pow(conj.k, 2.0f) + pow(conj.r, 2.0f);
#ifdef _DEBUG
			if(std::fabs(magSqrd) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			return Quaternion(conj.i / magSqrd, conj.j / magSqrd, conj.k / magSqrd, conj.r / magSqrd);
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Conjugate() const
		{
			return Quaternion<T>(-i, -j, -k, r);
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Normal() const
		{
			const T mag = sqrt(pow(i, 2.0f) + pow(j, 2.0f) + pow(k, 2.0f) + pow(r, 2.0f));
#ifdef _DEBUG
			if(std::fabs(mag) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			return Quaternion(i / mag, j / mag, k / mag, r / mag);
		}

		template <typename T>
		void Quaternion<T>::Normalize()
		{
			const T mag = sqrt(pow(i, 2.0f) + pow(j, 2.0f) + pow(k, 2.0f) + pow(r, 2.0f));
#ifdef _DEBUG
			if(std::fabs(mag) < NEARLY_ZERO)
			{
				throw std::runtime_error("Error! Dividing by nearly zero!");
			}
#endif
			i /= mag;
			j /= mag;
			k /= mag;
			r /= mag;
		}

		template<typename T>
		std::string Quaternion<T>::ToString() const
		{
			std::ostringstream oss;
			oss << "(" << i << ", " << j << ", " << k << ", " << r << ")";
			return oss.str();
		}

		template <typename T>
		Vector3<T> Quaternion<T>::RotateVector(const Vector3<T>& v) const
		{
			Quaternion<T> q1(*this);
			Quaternion<T> q2(v.x, v.y, v.z, 0);
			Quaternion<T> result = q1 * q2 * Quaternion::Inverse(q1);
			return Vector3<T>(result.i, result.j, result.k);
		}

		template<typename T>
		inline Matrix3<T> Quaternion<T>::RotationMatrix3() const
		{
			return Matrix3<T>(
				1.0f - 2.0f * (j * j + k * k), 2.0f * (i * j - k * r), 2.0f * (i * k + j * r),
				2.0f * (i * j + k * r), 1.0f - 2.0f * (i * i + k * k), 2.0f * (j * k - i * r),
				2.0f * (i * k - j * r), 2.0f * (j * k + i * r), 1.0f - 2.0f * (i * i + j * j)
			);
		}

		template<typename T>
		inline Matrix4<T> Quaternion<T>::RotationMatrix4() const
		{
			return Matrix4<T>(
				1.0f - 2.0f * (j * j + k * k), 2.0f * (i * j - k * r), 2.0f * (i * k + j * r), 0.0f,
				2.0f * (i * j + k * r), 1.0f - 2.0f * (i * i + k * k), 2.0f * (j * k - i * r), 0.0f,
				2.0f * (i * k - j * r), 2.0f * (j * k + i * r), 1.0f - 2.0f * (i * i + j * j), 0.0f,
				0.0f, 0.0f, 0.0f, 1.0f
			);
		}

		template<typename T>
		inline void Quaternion<T>::RotateAxis(const Vector3<T>& axis, const T& angle)
		{
			// Derived-from: https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
			Vector3<T> norm = axis.Normal();
			const T radians = angle * DEGREES_TO_RADIANS;
			const T sinHalfAngle = sin(radians/2.0f);
			const T cosHalfAngle = cos(radians/2.0f);
			i = norm.x * sinHalfAngle;
			j = norm.y * sinHalfAngle;
			k = norm.z * sinHalfAngle;
			r = cosHalfAngle;
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Inverse(const Quaternion<T>& q)
		{
			return q.Inverse();
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Conjugate(const Quaternion<T>& q)
		{
			return q.Conjugate();
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Normalize(const Quaternion<T>& q)
		{
			return q.Normal();
		}

		template <typename T>
		T Quaternion<T>::AngularDistance(const Quaternion& q1, const Quaternion& q2)
		{
			// Derived-from: https://math.stackexchange.com/questions/90081/quaternion-distance

			Quaternion norm1 = q1.Normal();
			Quaternion norm2 = q2.Normal();

			const T innerProduct = (norm1.i * norm2.i + norm1.j * norm2.j + norm1.k * norm2.k + norm1.r * norm2.r);
			return acos((2 * innerProduct * innerProduct) - 1.0f);
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::Lerp(const Quaternion& q1, const Quaternion& q2, const T t)
		{
			// Derived-from: https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=4736a2ea55426919d7c47ed78c33e5a153c40f6e

			const T lerp = fmax(0.0f, fmin(t, 1.0f));
			const T oneMinusLerp = 1.0f - lerp;
			return (oneMinusLerp * q1) + (lerp * q2);
		}

		template <typename T>
		Quaternion<T> Quaternion<T>::SLerp(const Quaternion& q1, const Quaternion& q2, const T t)
		{
			// Derived-from: https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=4736a2ea55426919d7c47ed78c33e5a153c40f6e

			T radians = Quaternion<T>::AngularDistance(q1, q2);

			Quaternion norm1 = q1.Normal();
			Quaternion norm2 = q2.Normal();

			const T slerp = fmax(0.0f, fmin(t, 1.0f));
			const T oneMinusSlerp = 1.0f - slerp;

			const T smallAngle = (45.0f * DEGREES_TO_RADIANS);

			if(fabs(radians) <= smallAngle)
			{
				// If angle is small enough just perform Linear Interpolation
				return (oneMinusSlerp * q1) + (slerp * q2);
			}

			const T sinAngle = sin(radians);
			const T sinSlerp = sin(slerp);
			const T sinOneMinusSlerp = sin(oneMinusSlerp);

			return (((sinOneMinusSlerp * radians) / sinAngle) * norm1) + (((sinSlerp * radians) / sinAngle) * norm2);
		}

		template<typename T>
		const Quaternion<T> Quaternion<T>::Identity = Quaternion<T>(0, 0, 0, 1);

		// ~~~~~~~~~~~~~~~~~~~~~~~ Unit Tests ~~~~~~~~~~~~~~~~~~~~~~~~~
#if _DEBUG
		template<typename T>
		void PrintQuaternion(const Math::Quaternion<T>& q, const std::string& name)
		{
			std::cout << name << ": " << q.ToString() << "\n";
		}

		template<typename T>
		void TestQuaternionOperations()
		{
			using namespace Math;

			// Test default constructor
			Quaternion<T> q1;
			Quaternion<T> expectedQ1(0, 0, 0, 0);
			PrintQuaternion(q1, "q1 (Default Constructor)");
			PrintQuaternion(expectedQ1, "Expected q1");

			// Test parameterized constructor
			Quaternion<T> q2(1, 2, 3, 4);
			Quaternion<T> expectedQ2(1, 2, 3, 4);
			PrintQuaternion(q2, "q2 (Parameterized Constructor)");
			PrintQuaternion(expectedQ2, "Expected q2");

			// Test copy constructor
			Quaternion<T> q3(q2);
			Quaternion<T> expectedQ3(1, 2, 3, 4);
			PrintQuaternion(q3, "q3 (Copy Constructor)");
			PrintQuaternion(expectedQ3, "Expected q3");

			// Test move constructor
			Quaternion<T> q4(std::move(q3));
			Quaternion<T> expectedQ4(1, 2, 3, 4);
			PrintQuaternion(q4, "q4 (Move Constructor)");
			PrintQuaternion(expectedQ4, "Expected q4");

			// Test assignment operator
			Quaternion<T> q5;
			q5 = q2;
			Quaternion<T> expectedQ5(1, 2, 3, 4);
			PrintQuaternion(q5, "q5 (Assignment Operator)");
			PrintQuaternion(expectedQ5, "Expected q5");

			// Test negation operator
			Quaternion<T> q6 = -q2;
			Quaternion<T> expectedQ6(-1, -2, -3, -4);
			PrintQuaternion(q6, "q6 (Negation Operator)");
			PrintQuaternion(expectedQ6, "Expected q6");

			// Test addition operator
			Quaternion<T> q7 = q2 + q2;
			Quaternion<T> expectedQ7(2, 4, 6, 8);
			PrintQuaternion(q7, "q7 (Addition Operator)");
			PrintQuaternion(expectedQ7, "Expected q7");

			// Test subtraction operator
			Quaternion<T> q8 = q2 - q2;
			Quaternion<T> expectedQ8(0, 0, 0, 0);
			PrintQuaternion(q8, "q8 (Subtraction Operator)");
			PrintQuaternion(expectedQ8, "Expected q8");

			// Test multiplication by scalar
			Quaternion<T> q9 = q2 * static_cast<T>(2);
			Quaternion<T> expectedQ9(2, 4, 6, 8);
			PrintQuaternion(q9, "q9 (Multiplication by Scalar)");
			PrintQuaternion(expectedQ9, "Expected q9");

			// Test multiplication by quaternion
			Quaternion<T> q10(0, 1, 0, 0);
			Quaternion<T> q11(0, 0, 1, 0);
			Quaternion<T> q12 = q10 * q11; // Should result in a quaternion representing a 90 degree rotation around z-axis
			Quaternion<T> expectedQ12(1, 0, 0, 0);
			PrintQuaternion(q12, "q12 (Multiplication by Quaternion)");
			PrintQuaternion(expectedQ12, "Expected q12");

			// Test inverse
			Quaternion<T> q13(-3, 0, 2, 1);
			Quaternion<T> invQ13 = q13.Inverse();
			Quaternion<T> expectedInvQ13(0.2142f, 0.0f, -0.1428f, 0.07142f); // Inverse of unit quaternion is itself
			PrintQuaternion(invQ13, "invQ13 (Inverse)");
			PrintQuaternion(expectedInvQ13, "Expected invQ13");

			// Test conjugate
			Quaternion<T> q14(1, 2, 3, 4);
			Quaternion<T> conjQ14 = q14.Conjugate();
			Quaternion<T> expectedConjQ14(-1, -2, -3, 4);
			PrintQuaternion(conjQ14, "conjQ14 (Conjugate)");
			PrintQuaternion(expectedConjQ14, "Expected conjQ14");

			// Test normalization
			Quaternion<T> q15(1, 2, 3, 4);
			q15.Normalize();
			T magQ15 = static_cast<T>(sqrt(1 * 1 + 2 * 2 + 3 * 3 + 4 * 4));
			Quaternion<T> expectedQ15(1 / magQ15, 2 / magQ15, 3 / magQ15, 4 / magQ15);
			PrintQuaternion(q15, "q15 (Normalized)");
			PrintQuaternion(expectedQ15, "Expected q15");

			// Test rotate vector
			Vector3<T> v1(1, 0, 0);
			Quaternion<T> q16(0, 0, 1, 0); // 180 degree rotation around y-axis
			Vector3<T> rotatedV1 = q16.RotateVector(v1);
			Vector3<T> expectedRotatedV1(-1, 0, 0);
			std::cout << "rotatedV1: (" << rotatedV1.x << ", " << rotatedV1.y << ", " << rotatedV1.z << ")\n";
			std::cout << "Expected rotatedV1: (" << expectedRotatedV1.x << ", " << expectedRotatedV1.y << ", " << expectedRotatedV1.z << ")\n";

			// Test angular distance
			Quaternion<T> q17(15, 0, 1, 0);
			Quaternion<T> q18(0, 1, -18, 0);
			T angularDistance = Quaternion<T>::AngularDistance(q17, q18);
			T expectedAngularDistance = static_cast<T>(PI);
			std::cout << "Angular Distance: " << angularDistance << " | Expected: " << expectedAngularDistance << "\n";

			// Test LERP
			Quaternion<T> q19(1, 1, 0, 1);
			Quaternion<T> q20(1, 1, 1, 1);
			Quaternion<T> lerpQ = Quaternion<T>::Lerp(q19, q20, static_cast<T>(0.5));
			Quaternion<T> expectedLerpQ(0.0f, 0.7071f, 0, 0.7071f);
			PrintQuaternion(lerpQ, "LERP Quaternion");
			PrintQuaternion(expectedLerpQ, "Expected LERP Quaternion");

			// Test SLERP
			Quaternion<T> slerpQ = Quaternion<T>::SLerp(q19, q20, static_cast<T>(0.5));
			Quaternion<T> expectedSlerpQ = Quaternion<T>::Normalize(Quaternion<T>(0, 1, 0, 1));
			PrintQuaternion(slerpQ, "SLERP Quaternion");
			PrintQuaternion(expectedSlerpQ, "Expected SLERP Quaternion");

			// Test Rotate Axis
			Vector3<T> axis(0, 1, 0);
			T angle = 180;
			std::cout << "rotationAxis: (" << axis.x << ", " << axis.y << ", " << axis.z << ")\n";
			std::cout << "angle: (" << angle<< ")\n";
			Quaternion<T> q21;
			q21.RotateAxis(axis, angle);
			Quaternion<T> expectedRotationQuat = Quaternion<T>(0, 1, 0, 0);
			PrintQuaternion(q21, "Rotate Axis");
			PrintQuaternion(expectedRotationQuat, "Expected Rotate Axis Quaternion");
		}
		inline void TestQuaternion()
		{
			std::cout << "Testing Quaternions with float type:\n";
			TestQuaternionOperations<float>();

			std::cout << "\nTesting Quaternions with double type:\n";
			TestQuaternionOperations<double>();
		}
#endif // _DEBUG

	}
}

#endif // MJOLNIR_QUATERNION_H