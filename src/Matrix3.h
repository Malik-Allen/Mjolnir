// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_MATRIX3_H
#define MJOLNIR_MATRIX3_H

#include "Vector3.h"

	// Shoutout and thank you to Scott
	/// These are the default vectors of the eye (camera) according to OpenGL and the right hand rule
	///								+Y   -Z
	///	                             |  /
	///   origin(0.0f,0.0f,0.0f);    | /
	///   up(0.0f,1.0f,0.0f);        |/_____+X
	///   forward(0.0f,0.0f,-1.0f);  /
	///                             /
	///                           +Z
	///

	/// Let's just make sure that all is clear about how this matrix is layed out. 

	/// How a matrix is really layed out is pretty much abitrary but we need to agree
	/// and the world has agreed (except for Microsoft) on the right-hand rule. 
	/// First, the 4x4 matrix is really just an array of 16 numbers.  
	/// We need to think of the array as a matrix in the following way
	/// 4x4 matrix - COLUMN MAJOR - The OpenGL, science, physics, mathematics and engineering way. 
	///	0  4  8  12        [0][0]  [1][0]  [2][0]  [3][0]   
	///	1  5  9  13  (or)  [0][1]  [1][1]  [2][1]  [3][1]
	///	2  6  10 14        [0][2]  [1][2]  [2][2]  [3][2]
	///	3  7  11 15        [0][3]  [1][3]  [2][3]  [3][3]

namespace Mjolnir
{
	namespace Math
	{
		template<typename T>
		struct Matrix3
		{
			~Matrix3() = default;

			explicit Matrix3(T s = 0)
			{
				if(s == 0)
				{
					SetIdentity();
				}
				else
				{
					m[0] = s;	m[1] = s;	m[2] = s;
					m[3] = s;	m[4] = s;	m[5] = s;
					m[6] = s;	m[7] = s;	m[8] = s;
				}
			}

			Matrix3(const Matrix3<T>& _m)
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];
				m[3] = _m[3];   m[4] = _m[4];   m[5] = _m[5];
				m[6] = _m[6];   m[7] = _m[7];   m[8] = _m[8];
			}

			Matrix3(Matrix3<T>&& _m) noexcept
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];
				m[3] = _m[3];   m[4] = _m[4];   m[5] = _m[5];
				m[6] = _m[6];   m[7] = _m[7];   m[8] = _m[8];

				_m[0] = T();    _m[1] = T();    _m[2] = T();
				_m[3] = T();    _m[4] = T();    _m[5] = T();
				_m[6] = T();    _m[7] = T();    _m[8] = T();
			}

			Matrix3(
				T xx, T xy, T xz,
				T yx, T yy, T yz,
				T zx, T zy, T zz)
			{
				m[0] = xx;  m[1] = xy;  m[2] = xz;
				m[3] = yx;  m[4] = yy;  m[5] = yz;
				m[6] = zx;  m[7] = zy;  m[8] = zz;
			}

			Matrix3(const Vector3<T>& v)
			{
				m[0] = v.x; m[1] = 0.0; m[2] = 0.0;
				m[3] = 0.0; m[4] = v.y; m[5] = 0.0;
				m[6] = 0.0; m[7] = 0.0; m[8] = v.z;
			}

			Matrix3& operator=(const Matrix3<T>& _m)
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];
				m[3] = _m[3];   m[4] = _m[4];   m[5] = _m[5];
				m[6] = _m[6];   m[7] = _m[7];   m[8] = _m[8];
				return *this;
			}

			Matrix3& operator=(Matrix3<T>&& _m) noexcept
			{
				m[0] = _m[0];   m[1] = _m[1];   m[2] = _m[2];
				m[3] = _m[3];   m[4] = _m[4];   m[5] = _m[5];
				m[6] = _m[6];   m[7] = _m[7];   m[8] = _m[8];

				_m[0] = T();    _m[1] = T();    _m[2] = T();
				_m[3] = T();    _m[4] = T();    _m[5] = T();
				_m[6] = T();    _m[7] = T();    _m[8] = T();
				return *this;
			}

			bool operator==(const Matrix3<T>& _m) const
			{
				return (m == _m.m);
			}

			bool operator!=(const Matrix3<T>& _m) const
			{
				return !(*this == _m);
			}

			T operator[](const unsigned int index) const
			{
				return *(m + index);
			}

			T& operator[](const unsigned int index)
			{
				return *(m + index);
			}

			Matrix3 operator*(const T& s) const
			{
				Matrix3 result = *this;
				result.m[0] *= s;   result.m[1] *= s;   result.m[2] *= s;
				result.m[3] *= s;   result.m[4] *= s;   result.m[5] *= s;
				result.m[6] *= s;   result.m[7] *= s;   result.m[8] *= s;
				return result;
			}

			Matrix3 operator*(const Matrix3& _m) const
			{
				Matrix3 result;
				for(int i = 0; i < 3; ++i)
				{
					for(int k = 0; k < 3; ++k)
					{
						result[i * 3 + k] = (m[0 * 3 + k] * _m[i * 3 + 0]) + (m[1 * 3 + k] * _m[i * 3 + 1]) + (m[2 * 3 + k] * _m[i * 3 + 2]);
					}
				}
				return result;
			}

			// Allows for multiplication with scalar as first variable i.e. scalar * vector = product;
			static friend Matrix3 operator*(const T& s, const Matrix3& _m)
			{
				return _m * s;
			}

			Matrix3& operator*=(const T& s)
			{
				*this = *this * s;
				return *this;
			}

			Matrix3& operator*=(const Matrix3& _m)
			{
				*this = *this * _m;
				return *this;
			}

			Vector3<T> operator*(const Vector3<T>& v) const
			{
				const T x = v.x * m[0] + v.y * m[3] + v.z * m[6];
				const T y = v.x * m[1] + v.y * m[4] + v.z * m[7];
				const T z = v.x * m[2] + v.y * m[5] + v.z * m[8];
				return Vector3<T>(x, y, z);
			}

			void SetIdentity()
			{
				m[0] = 1;	m[1] = 0;	m[2] = 0;
				m[3] = 0;	m[4] = 1;	m[5] = 0;
				m[6] = 0;	m[7] = 0;	m[8] = 1;
			}

			void SetIdentity(T s)
			{
				m[0] = s;	m[1] = 0;	m[2] = 0;
				m[3] = 0;	m[4] = s;	m[5] = 0;
				m[6] = 0;	m[7] = 0;	m[8] = s;
			}

			std::string ToString() const;

			static const Matrix3 Identity;

		private:
			T m[9];
		};

		template<typename T>
		std::string Matrix3<T>::ToString() const
		{
			std::ostringstream oss;
			for (int i = 0; i < 3; ++i) {
				oss << "(";
				for (int j = 0; j < 3; ++j) {
					oss << m[i * 3 + j];
					if (j < 2) oss << ", ";
				}
				oss << ")\n";
			}
			return oss.str();
		}

		template<typename T>
		const Matrix3<T> Matrix3<T>::Identity = Matrix3<T>(
			1, 0, 0, 
			0, 1, 0, 
			0, 0, 1);

		// ~~~~~~~~~~~~~~~~~ Unit Tests ~~~~~~~~~~~~~~~~~~~~~~~~

#if _DEBUG
	// Helper function to print matrix and compare with expected values
		template<typename T>
		void PrintMatrix(const Matrix3<T>& matrix, const std::string& name, const Matrix3<T>& expected)
		{
			std::cout << name << ":\n";
			bool pass = true;
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
				{
					if(matrix[i * 3 + j] != expected[i * 3 + j])
					{
						pass = false;
					}
				}
			}
			std::cout << matrix.ToString();
			std::cout << "Expected Matrix:\n";
			std::cout << expected.ToString();
			std::cout << (pass ? "Pass" : "Fail") << "\n\n";
		}

		template<typename T>
		void TestMatrix3Type()
		{
			// Test default constructor
			Matrix3<T> mat1;
			Matrix3<T> expected1(static_cast<T>(0));
			PrintMatrix(mat1, "Default Constructor", expected1);

			// Test constructor with scalar
			Matrix3<T> mat2(static_cast<T>(5));
			Matrix3<T> expected2(
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5),
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5),
				static_cast<T>(5), static_cast<T>(5), static_cast<T>(5)
			);
			PrintMatrix(mat2, "Constructor with Scalar", expected2);

			// Test copy constructor
			Matrix3<T> mat3 = mat2;
			PrintMatrix(mat3, "Copy Constructor", expected2);

			// Test move constructor
			Matrix3<T> mat4 = std::move(mat3);
			PrintMatrix(mat4, "Move Constructor", expected2);

			// Test parameterized constructor
			Matrix3<T> mat5(
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3),
				static_cast<T>(4), static_cast<T>(5), static_cast<T>(6),
				static_cast<T>(7), static_cast<T>(8), static_cast<T>(9)
			);
			Matrix3<T> expected5(
				static_cast<T>(1), static_cast<T>(2), static_cast<T>(3),
				static_cast<T>(4), static_cast<T>(5), static_cast<T>(6),
				static_cast<T>(7), static_cast<T>(8), static_cast<T>(9)
			);
			PrintMatrix(mat5, "Parameterized Constructor", expected5);

			// Test assignment operator
			Matrix3<T> mat6;
			mat6 = mat5;
			PrintMatrix(mat6, "Assignment Operator", expected5);

			// Test move assignment operator
			Matrix3<T> mat7;
			mat7 = std::move(mat6);
			PrintMatrix(mat7, "Move Assignment Operator", expected5);

			// Test identity loading
			Matrix3<T> mat8;
			mat8.SetIdentity();
			Matrix3<T> expected8(
				static_cast<T>(1), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(1), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(1)
			);
			PrintMatrix(mat8, "Set Identity", expected8);

			// Test identity loading with scalar
			Matrix3<T> mat9;
			mat9.SetIdentity(static_cast<T>(3));
			Matrix3<T> expected9(
				static_cast<T>(3), static_cast<T>(0), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(3), static_cast<T>(0),
				static_cast<T>(0), static_cast<T>(0), static_cast<T>(3)
			);
			PrintMatrix(mat9, "Set Identity with Scalar", expected9);

			// Test matrix multiplication
			Matrix3<T> mat10 = mat5 * mat8;
			PrintMatrix(mat10, "Matrix Multiplication", expected5); // mat5 * Identity should be mat5

			// Test matrix multiplication with scalar
			Matrix3<T> mat11 = mat5 * static_cast<T>(2);
			Matrix3<T> expected11(
				static_cast<T>(2), static_cast<T>(4), static_cast<T>(6),
				static_cast<T>(8), static_cast<T>(10), static_cast<T>(12),
				static_cast<T>(14), static_cast<T>(16), static_cast<T>(18)
			);
			PrintMatrix(mat11, "Matrix Multiplication with Scalar", expected11);

			// Test compound multiplication assignment
			Matrix3<T> mat12 = mat5;
			mat12 *= mat8;
			PrintMatrix(mat12, "Compound Multiplication Assignment", expected5); // mat5 * Identity should be mat5
		}

		inline void TestMatrix3()
		{
			std::cout << "Testing Matrix3 with float type:\n";
			TestMatrix3Type<float>();

			std::cout << "\nTesting Matrix3 with double type:\n";
			TestMatrix3Type<double>();

			std::cout << "\nTesting Matrix3 with int type:\n";
			TestMatrix3Type<int>();
		}
#endif // _DEBUG

		

}
}

#endif // MJOLNIR_MATRIX3_H