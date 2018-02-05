#include "Molecular_Model.hpp"

void Model_Segment::Velocity_Verlet_Step()
{
    for(int i=0; i<nParticle; i++)
    {
        Node *Temp = &Segment[i];
        Temp->coordinate[0] += Temp->velocity[0]*deltaT + Temp->acceleration[0]*deltaT2_half;
         Temp->coordinate[1] += Temp->velocity[1]*deltaT + Temp->acceleration[1]*deltaT2_half;
         Temp->coordinate[2] += Temp->velocity[2]*deltaT + Temp->acceleration[2]*deltaT2_half;
        
        Temp->velocity[0] += Temp->acceleration[0] * deltaT_half;
        Temp->velocity[1] += Temp->acceleration[1] * deltaT_half;
        Temp->velocity[2] += Temp->acceleration[2] * deltaT_half;
        
        Temp->acceleration[0] = -zeta * segment_mass * Temp->velocity[0];
        Temp->acceleration[1] = -zeta * segment_mass * Temp->velocity[1];
        Temp->acceleration[2] = -zeta * segment_mass * Temp->velocity[2];
    }
}

void Model_Segment::Velocity_Verlet_After_Step()
{
    double e_kin = 0;
    double vCM[3] = {0, 0, 0};
    for(int i=0;i<nParticle;i++)
    {
        Segment[i].velocity[0] += Segment[i].acceleration[0]*deltaT_half;
        Segment[i].velocity[1] += Segment[i].acceleration[1]*deltaT_half;
        Segment[i].velocity[2] += Segment[i].acceleration[2]*deltaT_half;
        
        e_kin += Get_V2(i);
        vCM[0] += Segment[i].velocity[0];
        vCM[1] += Segment[i].velocity[1];
        vCM[2] += Segment[i].velocity[2];
    }
    kin_step += 0.5 * e_kin * segment_mass;
    vCM[0] *= inv_nParticle;
    vCM[1] *= inv_nParticle;
    vCM[2] *= inv_nParticle;
    vCM_Total += sqrt(vCM[0]*vCM[0] + vCM[1]*vCM[1] + vCM[2]*vCM[2]);
}
