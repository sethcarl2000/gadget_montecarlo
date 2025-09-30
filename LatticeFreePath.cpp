#include <LatticeFreePath.hpp>
#include <vector>
#include <utility> 
#include <limits> 
#include <cmath> 
#include <cstdio> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TRandom3.h> 
#include <TGraph.h>

using namespace std; 
using GadgetUtils::DistanceToSphere; 


//find the intersect of this ray from within the cube of side length 1.  
pair<Vec3,Vec3> FindCubeIntersect(const Vec3& pos, const Vec3& dir, double side_length)
{   
    Vec3 best_norm; 
    double t = 1e30; 

    side_length *= 0.50; 

    for (const auto& norm : gCube_normals) {

        double t_side = (side_length - (norm * pos)) / (norm * dir); 

        //if this is the smallest 't' we have found. 
        if (t_side > 0. && t_side < t) { best_norm = norm; t = t_side; } 
    }

    return {pos + (dir * t), best_norm};
}

//______________________________________________________________________________________________________________
double LatticePathLength(Vec3 pos, Vec3 dir, double sphere_R, double lattice_spacing, long int max_iterations)
{
    long int iteration=0; 

    double R2 = pow( sphere_R, 2 ); 

    //if the position is inside the sphere, then quit. 
    if (pos.mag() < sphere_R) return std::numeric_limits<double>::quiet_NaN(); 

    double length = 0.; 

    const double dt = 1e-6; 

    while (iteration++ < max_iterations) {

        /*printf("it %4li pos {%+.3f %+.3f %+.3f} dir {%+.3f %+.3f %+.3f}\n", 
            iteration, 
            pos.x(), pos.y(), pos.z(),
            dir.x(), dir.y(), dir.z()
        );*/ 

        double length_to_sphere = DistanceToSphere( pos, dir, R2 ); 

        //we have intersected with the sphere
        if (length_to_sphere == length_to_sphere && length_to_sphere > 0.) { 
            
            return length + length_to_sphere; 

        } else { //we have NOT intersected the sphere. 

            auto pos_and_norm = FindCubeIntersect(pos, dir, lattice_spacing); 
            
            length += (pos_and_norm.first + (pos*-1.)).mag(); 

            //loop this track to the opposite side of the crystal. 
            pos = pos_and_norm.first + (pos_and_norm.second * -lattice_spacing); 

            /*if (fabs(pos.x()) > lattice_spacing/2. || 
                fabs(pos.y()) > lattice_spacing/2. || 
                fabs(pos.z()) > lattice_spacing/2.)   {
                    
                printf("error - exit of lattice. pos{%+.4f, %+.4f, %.4f}\n", pos.x(), pos.y(), pos.z()); 
            } */ 

            length += dt; 

            pos = pos + (dir * dt); 
        } 
    }

    //if we got here, we didn't find an intercept within 'max_iterations' 
    return std::numeric_limits<double>::quiet_NaN(); 
}
//______________________________________________________________________________________________________________
int SimulateLatticePath()
{
    TRandom3 rand; 

    double sphere_rad = 0.5; 
    double lattice_spacing = 1.; 

    const size_t n_samples  = 100; 
    const int n_simulate    = 1e5;
    const double step       = 0.9616; 

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

            double length = LatticePathLength(pos, dir, sphere_rad, lattice_spacing, 1e5); 

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
