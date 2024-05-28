// MIT License, Copyright (c) 2024 Malik Allen

#ifndef NON_LINEAR_TRANSFORMS_H
#define NON_LINEAR_TRANSFORMS_H

#include <cmath>

namespace Mjolnir
{
	namespace Math
	{
		// Inspired-by: Fast & Funky 1D Non-linear Transformations - GDC 2015 Talk https://youtu.be/mr5xkf6zSzk?si=Ln33jYuoIl_dQdVj

		// Library of efficient mathematical curves
		// All passed values are expected to be in the inclusive range of [0,...1]
		template<typename T>
		struct NonLinearTransforms1D
		{
			template<typename T>
			inline static T SmoothStart2(T x)
			{
				return x * x;
			}

			template<typename T>
			inline static T SmoothStart3(T x)
			{
				return x * x * x;
			}

			template<typename T>
			inline static T SmoothStart4(T x)
			{
				return x * x * x * x;
			}

			template<typename T>
			inline static T SmoothStart5(T x)
			{
				return x * x * x * x * x;
			}

			template<typename T>
			inline static T Flip(T x)
			{
				return 1.0f - x;
			}

			template<typename T>
			inline static T SmoothStop2(T x)
			{
				return Flip(SmoothStart2(Flip(x)));
			}

			template<typename T>
			inline static T SmoothStop3(T x)
			{
				return Flip(SmoothStart3(Flip(x)));
			}

			template<typename T>
			inline static T SmoothStop4(T x)
			{
				return Flip(SmoothStart4(Flip(x)));
			}

			template<typename T>
			inline static T SmoothStop5(T x)
			{
				return Flip(SmoothStart5(Flip(x)));
			}

			template<typename T>
			inline static T Mix(T a, T b, T x)
			{
				return (a + x * (b - a));
			}

			template<typename T>
			inline static T BlendStep2(T x)
			{
				return Mix(SmoothStart2(x), SmoothStop2(x), x);
			}

			template<typename T>
			inline static T BlendStep3(T x)
			{
				return Mix(SmoothStart3(x), SmoothStop3(x), x);
			}

			template<typename T>
			inline static T BlendStep4(T x)
			{
				return Mix(SmoothStart4(x), SmoothStop4(x), x);
			}

			template<typename T>
			inline static T BlendStep5(T x)
			{
				return Mix(SmoothStart5(x), SmoothStop5(x), x);
			}

			template<typename T>
			inline static T Arch(T x)
			{
				return x * (Flip(x));
			}

			template<typename T>
			inline static T SmoothStartArch(T x)
			{
				return x * Arch(x);
			}

			template<typename T>
			inline static T SmoothStopArch(T x)
			{
				return Flip(x) * Arch(x);
			}

			template<typename T>
			inline static T BlendStepArch(T x)
			{
				return Mix(SmoothStartArch(x), SmoothStopArch(x), x);
			}

			template<typename T>
			inline static T BellCurve(T x)
			{
				return SmoothStopArch(x) * SmoothStartArch(x);
			}

			template<typename T>
			inline static T BounceClampBottom(T x)
			{
				return std::fabs(x);
			}

			template<typename T>
			inline static T BounceClampTop(T x)
			{
				return Flip(std::fabs(Flip(x)));
			}

			template<typename T>
			inline static T BounceClampBottomTop(T x)
			{
				return BounceClampTop(BounceClampBottom(x));
			}

			// Cubic (3rd) Bezier curve through A,B,C,D where A (start) and D (end) are assumed to be 0 and 1
			template<typename T>
			inline static T NormalizedBezier3(T b, T c, T x)
			{
				// Derived-from: https://en.wikipedia.org/wiki/B%C3%A9zier_curve
				T s = 1.0f - x;
				T x2 = x * x;
				T s2 = s * s;
				T x3 = x2 * x;
				return (3.0f * b * s2 * x) + (3.0f * c * s * x2) + x3;
			}

			// Cubic (3rd) Hermite curve through A,B,C,D where A (start) and D (end) are assumed to be 0 and 1
			inline static T NormalizedHermite3(T b, T c, T x)
			{
				// Derived-from: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
				T x2 = x * x;
				T x3 = x2 * x;

				T m = b * (x3 - (2.0f * x2) + x);
				T p = c * (x3 - x2);
				return m + p;
			}

		private:
			NonLinearTransforms1D() = delete;
			~NonLinearTransforms1D() = delete;
			NonLinearTransforms1D(const NonLinearTransforms1D& NonLinearTransforms1D) = default;
			NonLinearTransforms1D(NonLinearTransforms1D&& NonLinearTransforms1D) = default;
			NonLinearTransforms1D& operator=(const NonLinearTransforms1D& NonLinearTransforms1D) = default;
			NonLinearTransforms1D& operator=(NonLinearTransforms1D&& NonLinearTransforms1D) = default;
		};
	}
}

typedef Mjolnir::Math::NonLinearTransforms1D<float> NonLinearTransforms1Df;
typedef Mjolnir::Math::NonLinearTransforms1D<double> NonLinearTransforms1Dd;

#endif // NON_LINEAR_TRANSFORMS_H