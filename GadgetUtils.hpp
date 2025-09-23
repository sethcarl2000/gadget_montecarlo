#ifndef GadgetUtils_H
#define GadgetUtils_H

#include <array> 
#include <cmath> 

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
        inline Vec3 operator+(const Vec3& rhs) const noexcept{
            return Vec3{{
                data[0] + rhs.data[0],
                data[1] + rhs.data[1],
                data[2] + rhs.data[2]
            }};
        } 
        
        //dot product operator
        double operator*(const Vec3& rhs) const noexcept {
            return 
                data[0] * rhs.data[0] +
                data[1] * rhs.data[1] + 
                data[2] * rhs.data[2]; 
        };
        
        //scalar multiplication operator
        inline Vec3 operator*(double mult) const noexcept {
            return Vec3{{
                data[0] * mult, 
                data[1] * mult,
                data[2] * mult
            }};
        }

        //magnitude of this vector 
        inline double mag() const noexcept {return std::sqrt((*this) * (*this));}; 

        //returns a normalized version of this Vec3 with mag()=1. 
        inline Vec3 unit() const noexcept {
            double mult = 1./mag(); 
            return (*this) * mult; 
        }; 
    };

    enum EventType : char {
        kNone           = 0, 
        kExit           = 1 << 0,
        kElastic_235    = 1 << 1,
        kElastic_238    = 1 << 2,
        kAbsorb_235     = 1 << 3
    };

    struct Neutron {
        Vec3 pos; 
        Vec3 momentum; 
        double energy;
        EventType event{kNone}; 
    };

    //compute the distance to a sphere centered at the origin, with radius sqrt(R2). if the ray defined as
    //  r(t)_i = x_i + s_i*t 
    //does NOT intersect with the sphere, then it returns nan. s
    double DistanceToSphere(const Vec3& X, const Vec3& S, double R2=1.); 
};


#endif 