#ifndef GadgetUtils_H
#define GadgetUtils_H

#include <array> 
#include <cmath> 
#include <Vec3.hpp> 
#include <vector> 
#include <TRandom3.h> 
#include "EventType.hpp" 
#include "Neutron.hpp"

//some utilities 
namespace GadgetUtils
{
    //compute the distance to a sphere centered at the origin, with radius sqrt(R2). if the ray defined as
    //  r(t)_i = x_i + s_i*t 
    //does NOT intersect with the sphere, then it returns nan. s
    double DistanceToSphere(const Vec3& X, const Vec3& S, double R2=1.); 

    void SimualteGenerations(
        const std::vector<EventType> event_types, 
        const int n_generations, 
        const int n_simulations,
        double sphere_rad=10., //units in cm 
        double number_density=4.82e22 //units in cm^-3. here is the value for natural, solid uranium
    );
    
    constexpr double barns_to_cm = 1.e-24; 
};


#endif 