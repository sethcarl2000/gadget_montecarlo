#include "EventType.hpp" 
#include "GadgetUtils.hpp"

using namespace std; 

//_______________________________________________________________________________________
vector<EventType> EventType::Init(double enrichment_frac, double num_density)
{
    //create all possible types of events (except leaving the sphere)

    vector<EventType> types; 

    // 235u elastic 
    types.push_back({
        kElastic, 
        [](double MeV){ return 8.4 * GadgetUtils::barns_to_cm; }, 
        enrichment_frac*num_density
    }); 

    // 235u fission 
    types.push_back({
        kFission, 
        [](double MeV){ return 1.6 * GadgetUtils::barns_to_cm; }, 
        enrichment_frac*num_density
    }); 

    // 238u elastic 
    types.push_back({
        kElastic, 
        [](double MeV){ return 9.4 * GadgetUtils::barns_to_cm; }, 
        (1. - enrichment_frac)*num_density
    }); 

    // 238u fission
    types.push_back({
        kFission, 
        [](double MeV){ return 0.6 * GadgetUtils::barns_to_cm; }, 
        (1. - enrichment_frac)*num_density
    }); 

    return types; 
}
//_______________________________________________________________________________________