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

    enum Flag {
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

public: 

    EventType(Flag event_type, std::function<double(double)> cross_section, double number_density=-1.)
     :  fEventFlag{event_type}, 
        fCrossSection{cross_section}, 
        fNumberDensity{number_density} {}; 
    
    ~EventType() {}; 
    
    //initialize all event types
    static std::vector<EventType> Init(double enrichment_frac = 0.800, double num_density=4.82e22); 

    Flag GetType() const { return fEventFlag; };

    //get the cross section
    double CrossSection(double MeV) const { return fCrossSection(MeV); }

    double MeanFreePath(double MeV) const { return 1./(fCrossSection(MeV) * fNumberDensity); }
};

#endif 