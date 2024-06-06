// MIT License, Copyright (c) 2024 Malik Allen

#include "Mjolnir.h"

using namespace Mjolnir;

void RunTests()
{
#if _DEBUG
    Math::TestVector2();
    Math::TestVector3();
    Math::TestVector4();
    Math::TestMatrix3();
    Math::TestMatrix4();
    Math::TestQuaternion();
    Math::TestTransform();
    Math::TestRay2D();
    Math::TestRay3D();
    Math::TestCircle();
    Math::TestSphere();
#endif    
}

int main()
{
    RunTests();
    return 0;
}