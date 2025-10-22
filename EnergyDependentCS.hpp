#ifndef EnergyDependentCS_H
#define EnergyDependentCS_H

#include <cmath> 

#include "EventType.hpp" 
#include "GadgetUtils.hpp"
#include <map> 
#include <string> 
#include <memory> 
#include <cmath> 

class EnergyDependentCS {
private:
    double MeV_min, MeV_max;
    double MeV_spacing_log; 
    double MeV_span_log; 

    //these points are spaced GEOMETRICALLY, so the (geometric) spacing between points is = pow( MeV_max/MeV_min, 1/points.size() ). 
    std::vector<double> points; 
public: 
    EnergyDependentCS(double _min, double _max, const std::vector<double>& _points); 

    //comptue the energy-dependent cross section (in barns)
    double Compute_CS(double MeV) const; 
}; 


class EnergyDependentCSManager {
private: 

    void Fill_EDCS_Uranium(); 

public: 
    //list of ptrs to energydependent CS objects
    std::map<std::string, EnergyDependentCS*> fEDCS; 

    EnergyDependentCSManager(); 
    ~EnergyDependentCSManager(); 

    std::vector<EventType> MakeEvents(const double enrichment_frac, const double number_density=4.82e22);   

}; 

#endif 