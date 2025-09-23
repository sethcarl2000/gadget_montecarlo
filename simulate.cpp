#include <iostream> 
#include <cmath> 
#include <GadgetUtils.hpp> 
#include <LatticeFreePath.hpp> 
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
    const int n_simulate = 1e6; 

    printf("simulating %i events...", n_simulate); cout << flush; 
    
    TRandom3 rand; 

    TH1D* hist = new TH1D("h", "distance to sphere edge", 200, 0., 30.); 

    vector<Neutron> event_buffer{}; 

    double sphere_rad = 0.5; 
    double lattice_spacing = 1.; 

    for (int i=0; i<n_simulate; i++) {

        Vec3 pos{{0., 0., 0.}}; 

        //generate a random position inside the sphere
        while (pos.mag() < sphere_rad) pos.data = {
            lattice_spacing*(0.5 - rand.Rndm()), 
            lattice_spacing*(0.5 - rand.Rndm()), 
            lattice_spacing*(0.5 - rand.Rndm())
        }; 

        //Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 
        Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 
        dir = dir.unit(); 

        double length = LatticePathLength(pos, dir, sphere_rad, lattice_spacing, 1e3); 

        /*if (length != length) { cout << "length: nan\n"; }
        else                  { 
            cout << "length: " << length << "\n"; 
        }*/ 

        if (length == length) hist->Fill( length ); 
    }

    cout << "done." << endl; 

    auto canv = new TCanvas("c", "Test hist"); 

    hist->Draw("E");
    hist->DrawCopy("SAME L");  

    canv->SaveAs("test.png"); 

    return 0; 
}