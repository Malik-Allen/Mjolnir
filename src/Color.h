// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_COLOR_H
#define MJOLNIR_COLOR_H

namespace Mjolnir
{
	namespace Math
	{
		struct Color
		{ 
			unsigned int r : 8;
			unsigned int g : 8;
			unsigned int b : 8;
			unsigned int a : 8;

			explicit Color(unsigned int _r, unsigned int _g, unsigned int _b, unsigned int _a = 255):
				r(_r),
				g(_g),
				b(_b),
				a(_a)
			{}

			Color(const Color& other):
				r(other.r),
				g(other.g),
				b(other.b),
				a(other.a)
			{}

			Color(Color&& other) noexcept:
				r(other.r),
				g(other.g),
				b(other.b),
				a(other.a)
			{
				other.r = other.g = other.b = other.a = 0;
			}

			~Color() = default;

			Color& operator=(const Color& other)
			{
				r = other.r;
				g = other.g;
				b = other.b;
				a = other.a;
				return *this;
			}

			Color& operator=(Color&& other) noexcept
			{
				if(this != &other)
				{
					r = other.r;
					g = other.g;
					b = other.b;
					a = other.a;
					other.r = other.g = other.b = other.a = 0;
				}
				return *this;
			}

			static const Color Black;
			static const Color Red;
			static const Color Blue;
			static const Color Green;
			static const Color White;
			static const Color Cyan;
			static const Color Magenta;
			static const Color Yellow;
		};
	}
}

typedef Mjolnir::Math::Color Color;

#endif // !MJOLNIR_COLOR_H