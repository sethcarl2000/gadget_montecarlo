#include <iostream> 
#include <cmath> 
#include <GadgetUtils.hpp> 
#include <LatticeFreePath.hpp> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TRandom3.h> 
#include <TGraph.h> 
#include <cstdio> 
#include <vector> 
#include <Neutron.hpp>
#include <EventType.hpp> 
#include <Vec3.hpp>
#include <iostream> 

using namespace std;  
using GadgetUtils::DistanceToSphere; 

int main(int argc, char* argv[]) 
{
    
    //
    vector<EventType> event_types = EventType::Init(0.8);

    GadgetUtils::SimualteGenerations(event_types, 7, 1e5, 12.5); 

    return 0;
}