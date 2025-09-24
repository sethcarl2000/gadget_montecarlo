#include <iostream> 
#include <cmath> 
#include <GadgetUtils.hpp> 
#include <LatticeFreePath.hpp> 
#include <TH1D.h> 
#include <TCanvas.h>
#include <TPad.h>  
#include <TGraph.h> 
#include <TF1.h> 
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
    
    TRandom3 rand; 

    double sphere_rad = 0.5; 
    double lattice_spacing = 1.; 

    const size_t n_samples  = 20; 
    const int n_simulate    = 1e3;
    const double step       = 0.95; 

    vector<double> vec_sphere_size{}; 
    vector<double> vec_avg_pathlen{}; 

    for (size_t sample=0; sample<n_samples; sample++) {

        vec_sphere_size.push_back(sphere_rad); 

        double avg = 0.; 
        int n_nan=0; 

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

            double length = LatticePathLength(pos, dir, sphere_rad, lattice_spacing, 1e6); 

            /*if (length != length) { cout << "length: nan\n"; }
            else                  { 
                cout << "length: " << length << "\n"; 
            }*/ 

            if (length == length) { avg += length; } else { n_nan++; } 
        }

        vec_avg_pathlen.push_back( avg / ((double)(n_simulate - n_nan)) ); 

        printf("point %2i/%i, sphere size: %.4f, fraction nan: %.4f, avg pathlen: %4.4f",
            sample, n_samples,  
            sphere_rad, 
            ((double)n_nan)/((double)n_simulate),
            vec_avg_pathlen.back()
        ); 

        //decrease the sphere size
        sphere_rad *= step; 

        cout << endl; 
    }
    
    cout << "done." << endl; 

    auto canv = new TCanvas("c", "Test hist"); 

    auto graph = new TGraph(vec_avg_pathlen.size(), vec_sphere_size.data(), vec_avg_pathlen.data()); 

    gPad->SetLogx(1); 
    gPad->SetLogy(1); 

    graph->SetMarkerStyle(kOpenCircle); 
    graph->Draw(); 

    auto mfp = new TF1("mfp", "1./(TMath::Pi()*x*x)", vec_sphere_size.back(), vec_sphere_size.front()); 

    mfp->Draw("SAME"); 

    canv->SaveAs("sphere_vs_pathlen.png"); 

    return 0; 
}