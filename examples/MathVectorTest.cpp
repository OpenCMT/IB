
#include "testing-prototype.h"

#include <iostream>
#include <math.h>


#include <Math/Dense.h>
#include <Math/Geometry.h>

using namespace uLib;

int Vector4f0(Vector4f c)
{
  c(3) = 0;
  if ( fabs(c(0)) < 0.001 && fabs(c(1)) < 0.001 && fabs(c(2)) < 0.001 )
    return 0;
  else
    return 1;
}


int main()
{
  BEGIN_TESTING(IBMath);

  //////////////////////////////////////////////////////////////////////////////
  ///////////////// VECTOR TESTING /////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  std::cout << "Testing IBVectors\n";

  HPoint3f p1(1,1,1);             // point with const float ctr
  HVector3f v1;  v1 << 1,1,1,1;   // vector with void ctr and comma assignment
  //  std::cout << "[test Hpoint - Hvector] : " << p1.transpose() << " - " <<
  //               v1.transpose() << " = " << (p1-v1).transpose() << "\n";
  TEST0( Vector4f0(p1 - v1));


  //////////////////////////////////////////////////////////////////////////////
  ///////////////// GEOMETRY TESTING ///////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  std::cout << "Testing IBGeometry\n";

  Geometry Geo;

  Geo.SetPosition(1,1,1);
  Geo.AddWorldEulerYZYRotation(0,0,0);
  HPoint3f wp = Geo.GetWorldPoint(HPoint3f(1,1,1));
  TEST0( Vector4f0(wp - HPoint3f(2,2,2)) );

  Geo.AddWorldEulerYZYRotation(M_PI_2,0,0);
  wp = Geo.GetWorldPoint(HPoint3f(1,1,1));
  TEST0( Vector4f0(wp - HPoint3f(2,2,0)) );

  Geo.AddWorldEulerYZYRotation(0,M_PI_2,0);
  wp = Geo.GetWorldPoint(HPoint3f(1,1,1));
  TEST0( Vector4f0(wp - HPoint3f(2,2,2)) );

  Geo.AddWorldEulerYZYRotation(0,0,M_PI_2);
  wp = Geo.GetWorldPoint(HPoint3f(1,1,1));
  //  std::cout << "Geometry matrix\n" << Geo.GetTransform() << "\n";
  //  std::cout << "World 1,1,1 coords\n" << wp << "\n";
  TEST0( Vector4f0(wp - HPoint3f(0,2,2)) );

  // TESTING FLIP AXES //

  Geo.SetPosition(0,0,0);
  Geo.AddWorldEulerYZYRotation(-M_PI_2,-M_PI_2,-M_PI_2); // reset previous
  Geo.AddWorldEulerYZYRotation(M_PI_2,0,0);              // PI_2 along X
  Geo.FlipAxes(0,2); // flip  X-Z
  HPoint3f p = Geo.GetWorldPoint(1,0,0);
  TEST0( Vector4f0(p - HVector3f(1,0,0)) ); // after flip and rotation X->X






  END_TESTING;
}


