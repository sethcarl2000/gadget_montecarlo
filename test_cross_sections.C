#include "GadgetUtils.hpp"
#include "EnergyDependentCS.hpp"
#include <TGraph.h> 
#include <TCanvas.h> 
#include <vector> 
#include <TPad.h> 

using namespace std; 

int test_cross_sections()
{
    auto edcs_mgr = new EnergyDependentCSManager(); 

    auto edcs_test = edcs_mgr->fEDCS["238u_Elastic"]; 
    
    double x_min = 1e-4; 
    double x_max = 20.; 

    double y_min = 1e-8; 
    double y_max = 20; 

    const int npts = 100; 

    double pts_x[npts], pts_y[npts]; 

    const double spacing = pow( x_max/x_min, 1./((double)npts-1) ); 

    double x = x_min; 
    for (int i=0; i<npts; i++) {

        pts_x[i] = x; 
        pts_y[i] = edcs_test->Compute_CS(x); 

        x *= spacing; 
    }

    auto c = new TCanvas;
    gPad->SetLogx(1); 

    auto graph = new TGraph(npts, pts_x, pts_y); 
    graph->SetTitle("Elastic cross section #sigma;#sigma (b);Incident Neutron Energy (MeV)"); 
    graph->Draw(); 

    //c->SaveAs("FissionCS_test.png");

    return 0; 
}

