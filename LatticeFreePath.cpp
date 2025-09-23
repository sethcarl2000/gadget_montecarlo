#include <LatticeFreePath.hpp>
#include <vector>
#include <utility> 
#include <limits> 
#include <cmath> 
#include <cstdio> 

using namespace std; 
using GadgetUtils::Vec3;
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

        } 
    }

    //if we got here, we didn't find an intercept within 'max_iterations' 
    return std::numeric_limits<double>::quiet_NaN(); 
}
