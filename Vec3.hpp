#ifndef Vec3_H
#define Vec3_H

//_________________________________________________________________
//
// class: Vec3
// Very basic class to quickly store / compute 3-vector information
//
//_________________________________________________________________

#include <array> 
#include <cmath> 

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


#endif 