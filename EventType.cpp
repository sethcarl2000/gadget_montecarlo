#include "EventType.hpp" 
#include "GadgetUtils.hpp"

using namespace std; 

//_______________________________________________________________________________________
vector<EventType> EventType::Init(double enrichment_frac, double num_density)
{
    //create all possible types of events (except leaving the sphere)

    vector<EventType> types; 

    types.push_back({
        kElastic_235, 
        [](double MeV){ return 8.4 * GadgetUtils::barns_to_cm; }, 
        enrichment_frac*num_density
    }); 

    types.push_back({
        kFission_235, 
        [](double MeV){ return 1.6 * GadgetUtils::barns_to_cm; }, 
        enrichment_frac*num_density
    }); 


    types.push_back({
        kElastic_238, 
        [](double MeV){ return 9.4 * GadgetUtils::barns_to_cm; }, 
        (1. - enrichment_frac)*num_density
    }); 

    types.push_back({
        kFission_238, 
        [](double MeV){ return 0.6 * GadgetUtils::barns_to_cm; }, 
        (1. - enrichment_frac)*num_density
    }); 

    return types; 
}
//_______________________________________________________________________________________