
#include "GadgetUtils.hpp"
#include <cmath>
#include <utility> 
#include <limits> 
#include <TRandom3.h> 
#include <string> 
#include <iostream> 

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

//___________________________________________________________________________________________________________
void GadgetUtils::SimualteGenerations(  const std::vector<EventType> event_types, 
                                        const int n_generations, 
                                        const int n_simulations,
                                        double sphere_rad,          //units in cm 
                                        double number_density       //units in cm^-3. 
                                        )
{
    //take the neutron and pick a random direction 
    const double neutron_energy = 2.0; //MeV
    const int neutrons_per_fission = 3; 

    const double R2 = pow( sphere_rad, 2 ); 

    TRandom3 rand; 

    //Get the flag of the event type 
    auto EventFlagName = [](EventType::Flag flag){
        string str; 
        switch (flag) {
            case EventType::kNone : str="none"; break; 
            case EventType::kElastic_235 : str="235u(elastic)"; break; 
            case EventType::kElastic_238 : str="238u(elastic)"; break; 
            case EventType::kAbsorb_235 : str="235u(absorb)"; break; 
            case EventType::kExit : str="exit-sphere"; break; 
            default : str="N/A"; 
        }
        return str; 
    };


    //const vector<EventType> event_types = EventType::Init(0.8, number_density); 

    //average number of neutrons at the end of 'n_generations' generations
    double N_avg =0.; 

    for (int trial=0; trial<n_simulations; trial++) {


        int n_neutrons_generate=1; 
        vector<Neutron> nextgen_buffer{};
        
        //________________________________________________________________________
        //generate 'n' neutrons
        auto GenerateNeutrons = [&](int n, Vec3 pos, double momentum){
            
            //generate a neutron
            for (int i=0; i<n; i++) {
                
                Vec3 P{rand.Gaus(), rand.Gaus(), rand.Gaus()}; 

                //find a random spot from within the sphere 
                nextgen_buffer.push_back({ pos, P.unit() * momentum });
            }
            return; 
        };   
        //________________________________________________________________________

        
        //________________________________________________________________________
        //generate 'n' 
        auto TransportNeutron = [&](Neutron& neutron, long int max_step=1e5){

            if (neutron.pos.mag() > sphere_rad) {
                fprintf(stderr, "Error - neutron ended up outside sphere!\n"); 
                return EventType::kNone; 
            }
        
            double MeV = neutron.momentum.mag(); 

            long int step=0; 
            while (step++ < max_step) {
                
                
                //compute all possible events
                double x_exit = GadgetUtils::DistanceToSphere(neutron.pos, neutron.momentum, R2); 
                //find which event type this is
                double min_x = x_exit; 
                EventType::Flag min_event = EventType::kExit; 

                
                for (const auto& ev : event_types) {

                    double x_ev = -log(1. - rand.Rndm()) * ev.MeanFreePath(MeV); 
                    
                    if (x_ev < min_x ) {
                        min_x = x_ev; 
                        min_event = ev.GetType(); 
                    }
                }
                

                //update the position
                neutron.total_path += min_x; 
                neutron.event       = min_event; 
                neutron.pos         = neutron.pos + (neutron.momentum.unit() * min_x); 
                neutron.n_collisions++; 

                /*printf("   neutron n.coll., pathlen, event-type: %u, %-3.5f, %s\n", 
                    step, 
                    neutron.total_path, 
                    EventFlagName(neutron.event).data()
                );*/  

                if (neutron.event == EventType::kAbsorb_235 || 
                    neutron.event == EventType::kExit) 
                    break; 

                //otherwise, change pick a new, random direction for the neutron 
                neutron.momentum[0] = rand.Gaus();
                neutron.momentum[1] = rand.Gaus();
                neutron.momentum[2] = rand.Gaus();

                neutron.momentum = neutron.momentum.unit() * MeV; 
            }

            return neutron.event; 
        };
        //________________________________________________________________________
        

        //pick a random starting place within the sphere
        Vec3 start_pos{sphere_rad*2., 0., 0.}; 
        while (start_pos.mag() > sphere_rad) {
            start_pos[0] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
            start_pos[1] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
            start_pos[2] = sphere_rad*( 1. - 2.*rand.Rndm() ); 
        }

        //add this first neutron
        GenerateNeutrons(1, start_pos, neutron_energy); 

        double N=0.; 

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
                
                //printf("  ~~~  transporting neutron %4i...\n", i_n++); 

                EventType::Flag end_state = TransportNeutron(neutron); 

                //check if the neutron ended with a fission
                if (end_state == EventType::kAbsorb_235) 
                    GenerateNeutrons(neutrons_per_fission, neutron.pos, neutron.momentum.mag()); 
                
            }
        }

        //check to see how many neutrons there are. 

        N_avg += (double)nextgen_buffer.size(); 
    }
    
    N_avg *= 1./((double)n_simulations); 

    printf("after %i, simulations of %i generations each, average is %.5f, so k=%.6f\n", 
        n_simulations,
        n_generations, 
        N_avg, 
        pow( N_avg, 1./((double)n_generations) ) 
    ); 


    //now, find the new kind of event. 
}   
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
