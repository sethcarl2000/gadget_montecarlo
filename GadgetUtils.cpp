
#include <GadgetUtils.hpp> 
#include <cmath>
#include <utility> 
#include <limits> 

using GadgetUtils::Vec3; 
using namespace std; 


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________


//___________________________________________________________________________________________________________
double GadgetUtils::DistanceToSphere(const Vec3& X, const Vec3& S, double R2)
{ 
    double s_x = (S * X)/S.mag(); 

    double x_x = X * X;

    //argument under the sqrt(...). if its < 0, then this ray does not intersect the circle (for t>0). 
    double sqrt_arg = s_x*s_x + (R2 - x_x); 

    if (sqrt_arg < 0.) return std::numeric_limits<double>::quiet_NaN(); 

    //the branch of the quadratic eqn. we want changes depending on whether or not we're inside the sphere
    return x_x > R2 ? - s_x - sqrt(sqrt_arg) : - s_x + sqrt(sqrt_arg);  
}   
//___________________________________________________________________________________________________________
