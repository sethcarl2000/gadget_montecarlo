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
#include "EnergyDependentCS.hpp"

using namespace std;  
using GadgetUtils::DistanceToSphere; 

int test_enrichment_level(const char* path_png) 
{
    //create the Energy-dependent CS manager 
    auto edcs_mgr = new EnergyDependentCSManager(); 

    struct DrawnGraphData_t {
        string legend_name; 
        double enrichment; 
        
        unsigned int color{1}, marker_style{0}; 
        TGraph *graph=nullptr; 
    }; 
    
    //
    const vector<DrawnGraphData_t> enrichment_levels{
        {"0.72 %", 0.00720, kBlack, kOpenCircle}, 
        {"20. %",    0.200, kBlue,  kOpenTriangleDown},  
        {"50. %",    0.500, kRed,   kOpenTriangleUp}, 
        {"85. %",    0.850, kBlack, kOpenSquare}, 
        {"100.%",    1.000, kBlue,  kFullCircle}
    };

    //units in cm 
    const double radius0 = 1.; 
    const double radius1 = 20.; 
    const int npts = 20; 

    const double uranium_density = 19.1e-3; //kg/cm^3 
    auto ball_mass = [uranium_density](double R) { return uranium_density * (TMath::Pi() * 4. / 3.)*pow( R, 3 ); }; 

    const int n_generations = 3; 
    const int n_simulations = 1e5; 
    
    auto canv = new TCanvas("c", "canvas", 1500, 650);
    canv->Divide(2,1);
     
    canv->cd(1);  
    auto hist_frame_rad = gPad->DrawFrame(0., 0., radius1*1.05, 3.); 
    hist_frame_rad->SetTitle(Form("K-value vs radius (%.2e trials of %i gens.);radius (cm);k-value", 
        (double)n_simulations, n_generations));

    canv->cd(2);  
    auto hist_frame_kg = gPad->DrawFrame(0., 0., ball_mass(radius1)*1.05, 3.); 
    hist_frame_kg->SetTitle(Form("K-value vs mass (%.2e trials of %i gens.);mass (kg);k-value", 
        (double)n_simulations, n_generations));
    
    hist_frame_kg->GetXaxis()->SetNdivisions(8);

    auto legend1 = new TLegend(0., 0.7, 0.5, 0.9);
    legend1->SetHeader("enrichment fraction");

    auto legend2 = new TLegend(0., 0.7, 0.5, 0.9);
    legend2->SetHeader("enrichment fraction");


    for (auto& g_data : enrichment_levels) {

        const double enrichment = g_data.enrichment;
        cout << "enrichment: " << enrichment << "~~~~~~~~~~~~~~~~~" << endl;
            
        vector<EventType> event_types = edcs_mgr->MakeEvents(enrichment); 
        //vector<EventType> event_types = EnergyDependentCS(enrichment); 

        vector<double> pts_radius, pts_mass, pts_k; 
        
        const double dr = (radius1 - radius0)/((double)npts-1); 
        double radius = radius0; 

        for (int i=0; i<npts; i++) {

            pts_radius.push_back( radius ); 
            pts_mass  .push_back( ball_mass(radius) ); 
            pts_k     .push_back( GadgetUtils::SimualteGenerations(event_types, n_generations, 1e5, radius) ); 
            
            printf("radius/k-value: %3.2f/%1.6f", radius, pts_k.back()); cout << endl; 
            radius += dr; 
        }
        auto graph = g_data.graph; 
        
        graph = new TGraph(pts_radius.size(), pts_radius.data(), pts_k.data()); 
        
        graph->SetTitle("k-value vs radius;radius;k-value"); 
        graph->SetMarkerColor(g_data.color); 
        graph->SetLineColor(g_data.color); 
        graph->SetMarkerStyle(g_data.marker_style); 
        graph->SetMarkerSize(0.5);

        canv->cd(1); 
        graph->Draw("SAME PL"); 
        
        legend1->AddEntry(graph, g_data.legend_name.c_str()); 
        
        //now make the graph for the mass plot
        graph = new TGraph(pts_mass.size(), pts_mass.data(), pts_k.data()); 
        
        graph->SetMarkerColor(g_data.color); 
        graph->SetLineColor(g_data.color); 
        graph->SetMarkerStyle(g_data.marker_style); 
        graph->SetMarkerSize(0.5);

        canv->cd(2); 
        graph->Draw("SAME PL"); 

        legend2->AddEntry(graph, g_data.legend_name.c_str()); 
        
    }

    canv->cd(1); legend1->Draw(); 
    canv->cd(2); legend2->Draw(); 

    canv->SaveAs(path_png); 

    return 0;
}