#include "EventType.hpp" 
#include "GadgetUtils.hpp"

using namespace std; 

//_______________________________________________________________________________________
vector<EventType> EventType::Init(double enrichment_frac, double num_density)
{
    //create all possible types of events (except leaving the sphere)

    vector<EventType> types; 

    //for now, i'm just going to use the approximate values for 2 MeV
    auto CS_235u_elastic = [](double MeV) {
        return 6.5 * GadgetUtils::barns_to_cm; 
    };  
    types.push_back({kElastic_235, CS_235u_elastic, enrichment_frac*num_density}); 

    auto CS_235u_fission = [](double MeV) {
        return 1.5 * GadgetUtils::barns_to_cm; 
    };
    types.push_back({kAbsorb_235, CS_235u_fission, enrichment_frac*num_density}); 

    auto CS_238u_elastic = [](double MeV) {
        return 8.0 * GadgetUtils::barns_to_cm; 
    };
    types.push_back({kElastic_238, CS_238u_elastic, (1. - enrichment_frac)*num_density}); 

    return types; 
}
//_______________________________________________________________________________________