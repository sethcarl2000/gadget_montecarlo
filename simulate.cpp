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
    vector<EventType> event_types = EventType::Init(0.500);

    vector<double> pts_radius, pts_k; 
    double radius = 1.; //units in cm 
    const double dr = 0.75; 
    const int npts = 30; 
    
    for (int i=0; i<npts; i++) {

        pts_radius.push_back( radius ); 
        pts_k     .push_back( GadgetUtils::SimualteGenerations(event_types, 4, 1e5, radius) ); 
        
        printf("radius/k-value: %3.2f/%1.6f", radius, pts_k.back()); cout << endl; 
        radius += dr; 
    }

    auto canv = new TCanvas("c", "canvas"); 

    auto graph = new TGraph(pts_radius.size(), pts_radius.data(), pts_k.data()); 
    graph->Draw(); 

    graph->SetTitle("k-value vs radius;radius;k-value"); 
    canv->SaveAs("k_vs_radius_enrich-0.500.png"); 

    return 0;
}