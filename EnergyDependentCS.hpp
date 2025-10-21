#ifndef EnergyDependentCS_H
#define EnergyDependentCS_H

#include <cmath> 

#include "EventType.hpp" 
#include "GadgetUtils.hpp"

//_______________________________________________________________________________________
std::vector<EventType> EnergyDependentCS(double enrichment_frac, double num_density=4.82e22)
{
    //create all possible types of events (except leaving the sphere)

    std::vector<EventType> types; 

    // 235u elastic 
    types.push_back({
        EventType::kElastic, 
        [](double MeV){  return 2.; 
            double Elog = std::log( MeV ) / 2.30258509; 
            return 5.19295 * std::exp( -0.564431 * Elog );
        }, 
        enrichment_frac*num_density
    }); 

    // 235u fission 
    types.push_back({
        EventType::kFission, 
        [](double MeV){ 
            double Elog = std::log( MeV ) / 2.30258509; 
            return (0.6 * Elog * Elog) + 1.; 
        }, 
        enrichment_frac*num_density
    }); 

    // 235u inelastic 
    /*types.push_back({
        EventType::kInelastic, 
        [](double MeV){ 
            double Elog = std::log( MeV ) / 2.30258509; 
            return std::exp( (-0.7*Elog*Elog) + (-0.6*Elog) + 0.6 );  
        }, 
        enrichment_frac*num_density
    });*/  

    // 238u elastic 
    types.push_back({
        EventType::kElastic, 
        [](double MeV){ return 2.; 
            double Elog = std::log( MeV ) / 2.30258509; 
            return 1.15 * 5.19295 * std::exp( -0.564431 * Elog );
        }, 
        (1. - enrichment_frac)*num_density
    });

    // 238u fission
    types.push_back({
        EventType::kFission, 
        [](double MeV){ 
            double Elog = std::log( MeV ) / 2.30258509; 
            return std::exp( (-1.*Elog*Elog) + (-3.3*Elog) + -2.2 );  
        }, 
        (1. - enrichment_frac)*num_density
    });

    // 238u inelastic 
    types.push_back({
        EventType::kInelastic, 
        [](double MeV){ 
            double Elog = std::log( MeV ) / 2.30258509; 
            return std::exp( (-1.3*Elog*Elog) + (-0.8*Elog) + 1. );  
        }, 
        (1. - enrichment_frac)*num_density
    });   

    return types; 
}
//_______________________________________________________________________________________

#endif 