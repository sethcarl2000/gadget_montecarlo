#include <iostream> 
#include <cmath> 
#include <GadgetUtils.hpp> 
#include <LatticeFreePath.hpp> 
#include <TH1D.h> 
#include <TCanvas.h> 
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
    const int n_simulate = 1e6; 
    const int max_it     = 1e6; 
    
    TRandom3 rand; 

    pair<double,double> draw_range{ 0., 15. }; 

    TH1D* hist_template = new TH1D("h", "path_length / lattice spacing", 125, draw_range.first, draw_range.second); 

    struct SphereSpacingPt { double rad; TH1D* hist=nullptr; };

    vector<SphereSpacingPt> points{
        {.rad = 0.500},
        {.rad = 0.250},
        {.rad = 0.125}
    };

    const double hist_scale_factor = (((double)n_simulate) * (draw_range.second-draw_range.first))/((double)hist_template->GetXaxis()->GetNbins()); 

    double sphere_rad = 0.5; 
    double lattice_spacing = 1.; 

    int i=0; 
    for (auto& point : points) {
        
        auto& hist = point.hist;

        int n_nan=0; 

        hist = (TH1D*)hist_template->Clone(Form("h_%i",++i)); 
        
        printf("(%i/%zi) simulating %i events...", i, points.size(), n_simulate); cout << flush; 

        for (int i=0; i<n_simulate; i++) {

            Vec3 pos{{0., 0., 0.}}; 

            //generate a random position inside the sphere
            while (pos.mag() < point.rad) pos.data = {
                lattice_spacing*(0.5 - rand.Rndm()), 
                lattice_spacing*(0.5 - rand.Rndm()), 
                lattice_spacing*(0.5 - rand.Rndm())
            }; 

            //Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 
            Vec3 dir{{rand.Gaus(), rand.Gaus(), rand.Gaus()}}; 
            dir = dir.unit(); 

            double length = LatticePathLength(pos, dir, point.rad, lattice_spacing, max_it); 

            if (length == length) { hist->Fill( length ); } else { n_nan++; } 
        }

        hist->Scale( 1./hist_scale_factor ); 

        printf("done. frac nan: %.4f", ((double)n_nan)/((double)n_simulate)); cout << endl; 
    }
    

    cout << "done." << endl; 

    auto canv = new TCanvas("c", "Test hist"); 

    gPad->SetLogy(1); 

    auto hist0 = points[0].hist; hist0->SetLineColor(kBlack);   hist0->Draw("E L");
    auto hist1 = points[1].hist; hist1->SetLineColor(kRed);     hist1->Draw("SAME E L");
    auto hist2 = points[2].hist; hist2->SetLineColor(kBlue);    hist2->Draw("SAME E L");

    auto mfp = new TF1("mfp", "TMath::Pi()*pow([0],2)*exp(-x * TMath::Pi() * pow([0],2) )", draw_range.first, draw_range.second);

    mfp->SetLineStyle(kDotted); mfp->SetLineWidth(1.); 
    
    mfp->SetParameter(0, points[0].rad); mfp->SetLineColor(kBlack); mfp->DrawCopy("SAME");  
    mfp->SetParameter(0, points[1].rad); mfp->SetLineColor(kRed);   mfp->DrawCopy("SAME");
    mfp->SetParameter(0, points[2].rad); mfp->SetLineColor(kBlue);  mfp->DrawCopy("SAME");
    
    canv->SaveAs("pathlen_vs_spheresize.png"); 

    return 0; 
}