#include <iostream> 
#include <cmath> 
#include <GadgetUtils.hpp> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TRandom3.h> 
#include <cstdio> 
#include <vector> 

using namespace std; 
using GadgetUtils::Vec3; 
using GadgetUtils::Neutron; 
using GadgetUtils::DistanceToSphere; 

int main(int argc, char* argv[]) 
{
    //
    const int n_simulate = 1e7; 

    printf("simulating %i events...", n_simulate); cout << flush; 
    
    TRandom3 rand; 

    TH1D* hist = new TH1D("h", "distance to sphere edge", 200, 0., 2.); 

    vector<Neutron> event_buffer{}; 

    for (int i=0; i<n_simulate; i++) {

        Vec3 pos{{2., 0., 0.}}; 

        //generate a random position inside the sphere
        while (pos.mag() > 1.) pos.data = {1. - 2.*rand.Rndm(), 1. - 2.*rand.Rndm(), 1. - 2.*rand.Rndm()}; 

        //Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 
        Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 

        hist->Fill( DistanceToSphere(pos, dir) ); 
    }

    cout << "done." << endl; 

    auto canv = new TCanvas("c", "Test hist"); 

    hist->Draw("E"); 

    canv->SaveAs("test.png"); 

    return 0; 
}