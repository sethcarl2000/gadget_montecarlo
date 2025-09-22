#ifndef GadgetUtils_H
#define GadgetUtils_H

#include <array> 
#include <cmath>
#include <utility> 
#include <limits> 

//some utilities 
namespace GadgetUtils
{
    struct Vec3 { 
        std::array<double,3> data;
        
        //direct access operator with no bounds checking
        double& operator[](unsigned int i) { return data[i]; }

        //modifiable element access
        double& x() noexcept { return data[0]; }
        double& y() noexcept { return data[1]; }
        double& z() noexcept { return data[2]; }

        //const element access
        double x() const noexcept { return data[0]; }
        double y() const noexcept { return data[1]; }
        double z() const noexcept { return data[2]; }

        //addition operator
        Vec3 operator+(const Vec3& rhs) const noexcept { 
            
            std::array<double,3>&& dat = { 
                data[0]+rhs.data[0],
                data[1]+rhs.data[1],
                data[2]+rhs.data[2] 
            };
            
            return Vec3{std::move(dat)}; 
        }
        
        //dot product operator
        double operator*(const Vec3& rhs) const noexcept {
            return 
                data[0] * rhs.data[0] +
                data[1] * rhs.data[1] + 
                data[2] * rhs.data[2];  
        }

        //magnitude of this vector 
        inline double mag() const noexcept { 
            return std::sqrt((*this) * (*this)); 
        }

        //returns a normalized version of this Vec3 with mag()=1. 
        Vec3 unit() const noexcept {

            double&& length = mag(); 
            
            std::array<double,3>&& dat = { 
                data[0]/ length,
                data[1]/ length,
                data[2]/ length
            };
            return Vec3{std::move(dat)}; 
        }
    };

    //compute the distance to a sphere centered at the origin, with radius sqrt(R2). if the ray defined as
    //  r(t)_i = x_i + s_i*t 
    //does NOT intersect with the sphere, then it returns nan. s
    double DistanceToSphere(const Vec3& X, const Vec3& S, double R2=1.) { 

        double s_x = (S * X)/S.mag(); 

        double x_x = X * X;

        //argument under the sqrt(...). if its < 0, then this ray does not intersect the circle (for t>0). 
        double sqrt_arg = s_x*s_x + (R2 - x_x); 

        if (sqrt_arg < 0.) return std::numeric_limits<double>::quiet_NaN(); 

        //the branch of the quadratic eqn. we want changes depending on whether or not we're inside the sphere
        return x_x > R2 ? - s_x - sqrt(sqrt_arg) : - s_x + sqrt(sqrt_arg);  
    }   
};


#endif 