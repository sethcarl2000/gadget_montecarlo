#include <iostream> 
#include <cmath> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TPad.h> 
#include <TF1.h> 
#include <TH1F.h> 
#include <TRandom3.h> 
#include <TGraph.h> 
#include <TLegend.h> 
#include <cstdio> 
#include <vector> 
#include <string> 
#include <iostream> 

#include "GadgetUtils.hpp" 
#include "Neutron.hpp"
#include "Vec3.hpp"
#include "EventType.hpp" 

using namespace std;  
using GadgetUtils::DistanceToSphere; 

int test_draw() 
{

    auto event_types = EventType::Init(1.00); 
    int n_generations = 2;  

    GadgetUtils::SimualteGenerations(event_types, n_generations, 1, 15., 4.82e22, false);
    
    return 0;
}