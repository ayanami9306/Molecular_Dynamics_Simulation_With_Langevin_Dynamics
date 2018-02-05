#include "Molecular_Model.hpp"

void Model_Segment::Single_Step()
{
    time_Now += deltaT;
    Velocity_Verlet_Step();
    Compute_Forces();
    Velocity_Verlet_After_Step();
}
