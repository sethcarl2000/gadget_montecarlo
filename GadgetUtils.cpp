#include "GadgetUtils.hpp"
#include <cmath>
#include <utility> 
#include <limits> 
#include <TRandom3.h> 
#include <string> 
#include <iostream> 
#include <TLine.h> 
#include <TEllipse.h> 
#include <TPad.h> 
#include <TObject.h> 
#include <TClass.h> 
#include <TText.h> 
#include <TH1F.h> 
#include <TCanvas.h> 
#include <thread> 
#include <mutex> 
#include <chrono> 
#include <random> 

using namespace std; 

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
double GadgetUtils::DistanceToSphere(const Vec3& X, const Vec3& S, double R2)
{ 
    double s_x = (S * X)/S.mag(); 

    double x_x = X * X;

    //argument under the sqrt(...). if its < 0, then this ray does not intersect the circle (for t>0). 
    double sqrt_arg = s_x*s_x + (R2 - x_x); 

    if (sqrt_arg < 0.) return std::numeric_limits<double>::quiet_NaN(); 

    //the branch of the quadratic eqn. we want changes depending on whether or not we're inside the sphere
    return x_x > R2 ? - s_x - sqrt(sqrt_arg) : - s_x + sqrt(sqrt_arg);  
}   

TCanvas *gDraw_canv=nullptr; 
TH1F *gDraw_pad=nullptr; 

//___________________________________________________________________________________________________________
double GadgetUtils::SimualteGenerations(    const std::vector<EventType> event_types, 
                                            const int n_generations, 
                                            const int n_simulations,
                                            double sphere_rad,          //units in cm 
                                            double number_density,       //units in cm^-3. 
                                            const bool draw_generations_gif
                                            )
{
    //take the neutron and pick a random direction 
    const double neutron_energy = 2.0; // MeV
    const int neutrons_per_fission = 3; 

    //this might change if I add a moderator
    const double M = 238. * Neutron::mass; 

    const double R2 = pow( sphere_rad, 2 ); 

    //const vector<EventType> event_types = EventType::Init(0.8, number_density); 

    //average number of neutrons at the end of 'n_generations' generations
    long int N_total =0.; 

    //call this wil a seed of '0', which tells TRandom3 to pick a random seed 
    TRandom3 rand; 

    for (int trial=0; trial<n_simulations; trial++) {

        int n_neutrons_generate=1; 
        vector<Neutron> nextgen_buffer{};

        //pick a random starting place within the sphere
        Vec3 start_pos{sphere_rad*2., 0., 0.}; 
        while (start_pos.mag() > sphere_rad) {
            start_pos[0] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
            start_pos[1] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
            start_pos[2] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
        }

        //add this first neutron
        GadgetUtils::GenerateNeutrons(nextgen_buffer, 1, start_pos, &rand); 

        double N=0.; 

        if (draw_generations_gif) {
            cout << "creating canvas..." << endl; 
            gDraw_canv = new TCanvas("c", "drawing", 600,600); 
            gDraw_pad = gDraw_canv->DrawFrame(-sphere_rad,-sphere_rad, sphere_rad,sphere_rad); 
        }

        for (int g=0; g<n_generations; g++) {

            //printf("<~~~  generation %4i...\n", g); 

            //if all the neutrons from the last generation left the sphere, then exit. 
            if (nextgen_buffer.empty()) break; 

            //make a copy of this vector, and clear it out. 
            vector<Neutron> neutrons{nextgen_buffer}; 
            nextgen_buffer.clear(); 

            int i_n=0; 
            //transport all the neutrons
            for (auto& neutron : neutrons) {
                neutron = GadgetUtils::TransportNeutron(&rand, neutron, nextgen_buffer, event_types, sphere_rad); 
            } 

            if (draw_generations_gif) {
                
                //write the generation name, and continue
                auto text = new TText(0.25, 0.10, Form("gen: %i", g)); 
                text->Draw("SAME"); 

                auto circ = new TEllipse(0.,0.,  sphere_rad,sphere_rad); 
                circ->SetLineColor(kBlue); 
                circ->SetLineWidth(2); 
                circ->Draw("SAME"); 

                //gDraw_canv->Modified(); 
                //gDraw_canv->Update();
                //gDraw_canv->SaveAs("output.gif+100"); 

                //gDraw_canv->DrawFrame(-sphere_rad,-sphere_rad, sphere_rad,sphere_rad); 
                
                /*/delete all the lines from past generations, and continue. 
                for (auto tobj : *(gDraw_canv->GetListOfPrimitives())) {
                    if (tobj->IsA() == TClass::GetClass<TLine>() || 
                        tobj->IsA() == TClass::GetClass<TEllipse>() || 
                        tobj->IsA() == TClass::GetClass<TText>()
                    ) delete tobj; 
                }*/ 
                //this_thread::sleep_for(1500ms); 
                cout << "done with gen " << g << " neutron buffer size: " << nextgen_buffer.size() << endl; 
            }

        }
        //check to see how many neutrons there are.
        N_total += (long int)nextgen_buffer.size(); 

    }//for (int trial=0; ...)
    
    double N_avg = ((double)N_total)/((double)n_simulations); 

    double k_effective = pow( N_avg, 1./((double)n_generations) );

    printf("after %i, simulations of %i generations each, average is %.5f, so k=%.6f\n", 
        n_simulations,
        n_generations, 
        N_avg, 
        k_effective 
    ); 

    //now, find the new kind of event. 
    return k_effective; 
}   
//___________________________________________________________________________________________________________
int GadgetUtils::ThrowPoisson(double rndm, double lambda)
{
    double n=0.; 
    double p = exp(-lambda);
    double cdf = p;   
    while( cdf < rndm ) {
        n += 1.; 
        p *= lambda / n;
        cdf += p;  
    }
    return (int)n; 
}
//___________________________________________________________________________________________________________
Neutron GadgetUtils::ElasticScatter(Neutron n, double M, TRandom3* rand)
{
    // First, generate a random, unit vector. 
    Vec3 u_rand = Vec3{
        rand->Gaus(), 
        rand->Gaus(), 
        rand->Gaus()
    }.unit(); 

    // Now, find the component of 'u_rand' the neutron's momentum perpendicular to the new unit vector
    Vec3 p_parallel = n.momentum.unit(); 

    auto p_perp = u_rand + ( u_rand * (u_rand * p_parallel) * -1. ); 
    p_perp = p_perp.unit(); 

    // now we have a unit vector in our original direction (p_parallel) and another one with a random phi, 
    // which is perpendicular to the original direction (p_perp).

    // Now, we must select a random cos(X). 
    // To simulate how elastic scattering becomes more forward-biased at higher energies, 
    // we will use an expoenential distribution. 
    double MeV = n.GetKineticEnergy(); 

    double angle_bias_param = 3.22291 * MeV; 

    // this generates cosX on the interval: cosX \in [-1, 1). 
    double cosX; 
    if (MeV > 0.1) {
        
        do { 
            cosX = -(1./angle_bias_param) * log( 1. - rand->Rndm() ); 
        } while (cosX > 2.); 

        cosX = 1. - cosX; 

    } else {
        
        //at these energies (MeV < 0.1), scattering is basically isotropic. 
        cosX = 1. - 2.*rand->Rndm(); 
    }

    // so, now we've chosen cosX (with 'X' being defiend as the angle between our new dir. and our old dir.).
    // Now, we must choose the new energy. when M >> m_n (scattering center mass bigger than neutron mass), 
    // the following is a good approx.: 
    double p_new_mag = (( M + (Neutron::mass*cosX) )/( M + Neutron::mass )) * n.momentum.mag(); 

    // so, now let's pick the new momentum: 
    Vec3 p_new = (p_parallel * cosX) + (p_perp * sqrt(1. - (cosX*cosX))); 

    n.momentum = p_new * p_new_mag; 

    return n; 
}
//___________________________________________________________________________________________________________
Neutron GadgetUtils::InelasticScatter(Neutron n, TRandom3* rand)
{
    // this assumes that inelastic scattering is isotropic, AND that the energy a neutron looses 
    // is randomly, uniformly distributed on the interval: [0, neutron_energy]. 

    // First, generate a random, unit vector. 
    Vec3 u_rand = Vec3{
        rand->Gaus(), 
        rand->Gaus(), 
        rand->Gaus()
    }.unit(); 

    double energy = n.GetKineticEnergy(); 

    double new_energy = energy - (min<double>( 3., energy ) * rand->Rndm()); 

    n.momentum = u_rand * sqrt( 2. * Neutron::mass * new_energy ); 

    return n; 
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GadgetUtils::GenerateNeutrons(vector<Neutron>& neutron_buffer, int N, Vec3 pos, TRandom3* rand)
{
    //when we use a maxwellian dist for the 'temp', we need this as a way to give all the neutrons proper momenta 
    const double PFNS_momentum_sigma = 35.22; // MeV/c 

    for (int i=0; i<N; i++) {

        Vec3 momentum{
            rand->Gaus() * PFNS_momentum_sigma, 
            rand->Gaus() * PFNS_momentum_sigma, 
            rand->Gaus() * PFNS_momentum_sigma
        }; 

        neutron_buffer.push_back({ pos, momentum }); 
    }
}
//___________________________________________________________________________________________________________
void GadgetUtils::Fission(vector<Neutron>& buffer, Neutron n, TRandom3* rand)
{
    //use this empirical formula to compute the number of fissile neutrons
    // using a poisson dist from 
    double n_avg = (0.1543 * n.GetKineticEnergy()) + 2.299; 

    int n_new_neutrons = GadgetUtils::ThrowPoisson(rand->Rndm(), n_avg); 

    GadgetUtils::GenerateNeutrons(buffer, n_new_neutrons, n.pos, rand); 

    return; 
}
//___________________________________________________________________________________________________________
Neutron GadgetUtils::TransportNeutron(
    TRandom3* rand,                         //the random number generator 
    Neutron n, 
    std::vector<Neutron>& neutron_buffer,   //the 'buffer' of the next generation's neutrons
    const std::vector<EventType>& events,   //the list of possible 'events' corresponding to this material 
    const double rad1,                      //the outer-radius of this material 
    const double rad0,                      //the inner-radius of this material. if rad0<0, then an 'inner raidus' is not considered. 
    bool* is_inside_exit                    //the fcn sets this to 'true' if the neutron leaves the hollow shell on the inner surface. 
                                            // otherwise (the neutron exits on the outer surface), this is set to 'false'.  
                                            // this is set to 'nullptr' if the final event-type is NOT an exit. 
)
{
    //check if the TRandom3 ptr passed is valid
    if (rand==nullptr) {
        throw logic_error("in <GadgetUtils::TransportNeutron>: TRandom3 ptr passed is null!");
    }

    //clear this ptr. we do this to make sure that this fcn can manage this ptr. 
    if (is_inside_exit!=nullptr) delete is_inside_exit; 

    //see whether or not this material is a 'hollow shell', with a material of a different type
    const bool simulate_inner_sphere = (rad0 > 0.); 

    double MeV = n.GetKineticEnergy(); 

    //arbitrary big number that would never be a reasonable physical result from any of the path-length computations below
    const double big_num = 1.e100;

    const unsigned int max_n_steps = 1e4; 

    long int step=0; 
    while (step++ < 1e4) {
        
        //compute all possible events.
        
        //first, we compute the length this neutron needs to exit this bulk material entirely 
        bool inside_exit; 
        double x_exit_out = GadgetUtils::DistanceToSphere(n.pos, n.momentum, rad1*rad1); 
        double x_exit_in;; 
        if (simulate_inner_sphere) {

            //compute the distance until this neutron exits this hollow shell (out of the inner surface)
            x_exit_in = GadgetUtils::DistanceToSphere(n.pos, n.momentum, rad0*rad0); 

            //if x_exit_in == NaN, then that means that this neutron misses the inner shell entirely. 
            if (x_exit_in != x_exit_in) x_exit_in = big_num; 
        
        } else { x_exit_in=big_num; }  

        double min_x; 
        if (x_exit_in < x_exit_out) { inside_exit=true;   min_x=x_exit_in; }
        else                        { inside_exit=false;  min_x=x_exit_out; }

        //a 'nullptr' corresponds to a type of 'exit' 
        EventType::Flag min_event_flag = EventType::kExit; 
        
        //nuclear mass corresponding to the species of nucleus for the 'min' even
        double min_event_mass=0.; 

        //now, get the event types which DON'T correspond to an exit 
        for (const auto& ev : events) {

            double x_ev = -log(1. - rand->Rndm()) * ev.MeanFreePath(MeV); 
            
            if (x_ev < min_x ) {
                min_x = x_ev; 
                min_event_flag = ev.GetType();
                min_event_mass = ev.GetMass();  
            }
        }

        //otherwise, update the position of this neutron 
        n.total_path += min_x;                                  //update the neutron's path len
        n.event       = min_event_flag;                         //update the neutron's event type           
        n.pos         = n.pos + (n.momentum.unit() * min_x);    //update the neutron's position
        n.n_collisions++; 

        //figure out what event type the 'min-x' event was. 
        switch (min_event_flag) {

            case EventType::kExit       : {
                //if we're simulating a hollow shell, then report whether the 'exit' occured on the inner or outer surface
                if (simulate_inner_sphere) is_inside_exit = new bool(inside_exit); 
                return n; //quit the transportation, we're done. 
            }
            case EventType::kFission    : {
                GadgetUtils::Fission(neutron_buffer, n, rand);
                return n; //quit the transportation, we're done. 
            }
            case EventType::kElastic    : { 
                n = GadgetUtils::ElasticScatter(n, min_event_mass, rand); 
                break; 
            }
            case EventType::kInelastic  : { 
                n = GadgetUtils::InelasticScatter(n, rand); 
                break; 
            }
            default                     : {
                fprintf(stderr, "Warning in <GadgetUtils::TransportNeutron>: unknown event flag encountered: %u\n", min_event_flag); 
                return n; 
            }
        }//switch (event_type)

    }//while (step++ < 1e4)

    //we should not get here! either the 'max_n_steps' is too small, or something freaky is going on with our cross-section model. 
    fprintf(stderr, "Warning in <GadgetUtils::TransportNeutron>: "
                    "maximum number of steps reached witout valid exit condition. max_n_steps=%u\n", max_n_steps); 
    return n; 
}
//___________________________________________________________________________________________________________
