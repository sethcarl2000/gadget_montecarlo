#ifndef Neutron_H
#define Neutron_H

//___________________________________________________________
//
//  struct: Neutron
//  basic struct with neutron info 
//
//____________________________________________________________
#include <Vec3.hpp>
#include <EventType.hpp> 

struct Neutron {
    Vec3 pos; 
    Vec3 momentum; 
    double total_path=0.; 
    unsigned short int n_collisions=0; 
    EventType::Flag event{EventType::kNone}; 
};

#endif 