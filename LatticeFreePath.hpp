#ifndef LatticeFreePath_H
#define LatticeFreePath_H

#include <GadgetUtils.hpp> 
#include <optional> 
#include <vector>
#include <Vec3.hpp> 

const std::vector<Vec3> gCube_normals{
    Vec3{{+1., 0., 0.}},
    Vec3{{-1., 0., 0.}},
    Vec3{{0., +1., 0.}},
    Vec3{{0., -1., 0.}},
    Vec3{{0., 0., +1.}},
    Vec3{{0., 0., -1.}}
}; 

//if this ray intersects with the square, a Vec3 of the intersect point is returned. if not, it returns nothing (hence std::optional)
std::pair<Vec3,Vec3> FindCubeIntersect(const Vec3& pos, const Vec3& dir, double side_length=1.); 

//for a given sphere size and cubic lattice spacing, and position relative to the center of the sphere, find the path length to the sphere. 
double LatticePathLength(   Vec3 pos, 
                            Vec3 dir, 
                            double sphere_R=0.5, 
                            double lattice_spacing=1., 
                            long int max_iterations=100000 ); 

//simulate many path-length calculations on the lattice. 
int SimulateLatticePath(); 



#endif 