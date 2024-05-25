// MIT License, Copyright (c) 2024 Malik Allen

#ifndef MJOLNIR_MATHCONSTANTS_H
#define MJOLNIR_MATHCONSTANTS_H

#include <cmath>
#include <string>
#include <sstream>

#ifndef NEARLY_ZERO
#define NEARLY_ZERO 1.0e-7f
#endif

#ifndef PI
#define PI 3.14159265358979323846f
#endif

#ifndef TWO_PI
#define TWO_PI (2.0f * PI)
#endif

#ifndef HALF_PI
#define HALF_PI (PI / 2.0f)
#endif

#ifndef DEGREES_TO_RADIANS
#define DEGREES_TO_RADIANS (PI / 180.0f)
#endif

#ifndef RADIANS_TO_DEGREES
#define RADIANS_TO_DEGREES (180.0f / PI)
#endif

#ifndef GOLDEN_RATIO
#define GOLDEN_RATIO ((1.0f + SQRT_5) / 2.0f)
#endif

#endif // MJOLNIR_MATHCONSTANTS_H