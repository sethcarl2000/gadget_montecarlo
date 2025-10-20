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

    //returns the computed k-value 
    double SimualteGenerations(
        const std::vector<EventType> event_types, 
        const int n_generations, 
        const int n_simulations,
        double sphere_rad=10., //units in cm 
        double number_density=4.82e22, //units in cm^-3. here is the value for natural, solid uranium, 
        const bool draw_generations_gif=false
    );

    //throws a random int according to a poisson dist, given lambda = <n>., and rndm is a random dist in [0, 1).  
    int ThrowPoisson(double rndm, double lambda); 

    //once an elastic scattering event has been determined, determine a new direction, energy, etc.
    //'M' is the mass of the nucleus doing the scattering. 
    Neutron ElasticScatter(Neutron n, double M, TRandom3* rand); 

    //simulate inelastic scattering. this is done by assuming that: 
    // 1. inelastic scattering is always isotropic, and 
    // 2. the amount of energy subtracted from the neutron is uniform and isotropic. 
    Neutron InelasticScatter(Neutron n, TRandom3* rand); 

    void GenerateNeutrons(std::vector<Neutron>& neutron_buffer, int N, Vec3 pos, TRandom3* rand); 

    void Fission(std::vector<Neutron>& buffer, Neutron n, TRandom3* rand); 

    constexpr double barns_to_cm = 1.e-24; 

};


#endif 