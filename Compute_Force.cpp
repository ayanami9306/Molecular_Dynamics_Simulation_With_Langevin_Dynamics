#include "Molecular_Model.hpp"

void Model_Segment::Compute_Forces()
{
    //neighbor-list method
    //atom pair constructed by table?->space inefficient
    //separate list?
    
    //check if reconstruction is necessary
    accu_movement += sqrt(vvMax)*deltaT;
    if(accu_movement > 0.5*radius_NebrShell)
    {
        accu_movement = 0;
        Build_NebrList();
    }
    vvMax = 0;
    
    //calculate force
    //Lennard-Jones Potential and Bond Stretching Potential
    pot_step += Calc_Nonbonded_Potential() + Calc_Bonded_Potential_AND_Apply_Langevin();
    
    // F = ma
    if(segment_mass != 1)
    {
        //#pragma omp parallel for
        for(int i=0; i<nParticle; i++)
        {
            Segment[i].acceleration[0] *= inv_segment_mass;
            Segment[i].acceleration[1] *= inv_segment_mass;
            Segment[i].acceleration[2] *= inv_segment_mass;
        }
    }
}
