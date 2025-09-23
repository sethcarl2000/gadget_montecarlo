#ifndef LatticeFreePath_H
#define LatticeFreePath_H

#include <GadgetUtils.hpp> 
#include <optional> 
#include <vector>

const std::vector<GadgetUtils::Vec3> gCube_normals{
    GadgetUtils::Vec3{{+1., 0., 0.}},
    GadgetUtils::Vec3{{-1., 0., 0.}},
    GadgetUtils::Vec3{{0., +1., 0.}},
    GadgetUtils::Vec3{{0., -1., 0.}},
    GadgetUtils::Vec3{{0., 0., +1.}},
    GadgetUtils::Vec3{{0., 0., -1.}}
}; 

//if this ray intersects with the square, a Vec3 of the intersect point is returned. if not, it returns nothing (hence std::optional)
std::pair<GadgetUtils::Vec3,GadgetUtils::Vec3> FindCubeIntersect(const GadgetUtils::Vec3& pos, const GadgetUtils::Vec3& dir, double side_length=1.); 

//for a given sphere size and cubic lattice spacing, and position relative to the center of the sphere, find the path length to the sphere. 
double LatticePathLength(   GadgetUtils::Vec3 pos, 
                            GadgetUtils::Vec3 dir, 
                            double sphere_R=0.5, 
                            double lattice_spacing=1., 
                            long int max_iterations=100000 ); 




#endif 