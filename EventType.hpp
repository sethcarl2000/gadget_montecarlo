#ifndef EventType_H
#define EventType_H

//_________________________________________________________________________________________
// 
//  class: EventType
//  Given a label below, it can calculate a cross section (or mean free path) for a given process
//
//_________________________________________________________________________________________

#include <vector> 
#include <utility> 
#include <functional> 

class EventType {
public: 

    enum Flag : int {
        kNone = 0,
        kExit,
        kElastic,
        kInelastic, 
        kFission
    }; 

private: 
    
    //type of event this process is 
    Flag fEventFlag{kNone};
    
    //for now, the cross section will be a dummy fuction
    std::function<double(double)> fCrossSection;
    
    //number density of the relevant type of nuclei, in 1/cm^3 
    double fNumberDensity{-1.}; 

    //nuclear mass of the particle of this event. this is needed for elastic-scattering energy loss computations.
    //units must be in MeV/c^2 
    double fMass; 

public: 

    EventType(Flag event_type, std::function<double(double)> cross_section, double number_density=-1., double mass=939.)
     :  fEventFlag{event_type}, 
        fCrossSection{cross_section}, 
        fNumberDensity{number_density}, 
        fMass{mass} {}; 
    
    ~EventType() {}; 
    
    //initialize all event types
    static std::vector<EventType> Init(double enrichment_frac = 0.800, double num_density=4.82e22); 

    Flag GetType() const { return fEventFlag; };

    //nuclear mass of the species of nucleus for this event. Units in MeV/c^2
    double GetMass() const { return fMass; } 

    //get the cross section
    double CrossSection(double MeV) const { return fCrossSection(MeV); }

    double MeanFreePath(double MeV) const { return 1./(fCrossSection(MeV) * fNumberDensity); }
};

#endif 