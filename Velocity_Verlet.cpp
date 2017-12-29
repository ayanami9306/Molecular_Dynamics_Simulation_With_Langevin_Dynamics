//
//  Velocity_Verlet.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Velocity_Verlet_Step()
{
    vvMax = 0;
    double e_kin = 0;
    double vCM[3] = {0, 0, 0};
//#pragma omp parallel for reduction(+:e_kin, vCM[0:3])
    {
        for(int i=0; i<nParticle; i++)
        {
            Node *Temp = &Segment[i];
            
            Temp->velocity[0] += Temp->acceleration[0]*deltaT_half;
            Temp->velocity[1] += Temp->acceleration[1]*deltaT_half;
            Temp->velocity[2] += Temp->acceleration[2]*deltaT_half;
            
            double temp_vv = Get_V2(i);
            e_kin += temp_vv;
            //#pragma omp critical(vvMax)
            if(temp_vv>vvMax) vvMax = temp_vv;
            vCM[0] += Segment[i].velocity[0];
            vCM[1] += Segment[i].velocity[1];
            vCM[2] += Segment[i].velocity[2];
            
            //calc position
            //x(t+deltaT) = x(t) + v(t)*deltaT + 0.5*a(t)*deltaT^2
            Temp->coordinate[0] += Temp->velocity[0]*deltaT + Temp->acceleration[0] * deltaT2_half;
            Temp->coordinate[1] += Temp->velocity[1]*deltaT + Temp->acceleration[1] * deltaT2_half;
            Temp->coordinate[2] += Temp->velocity[2]*deltaT + Temp->acceleration[2] * deltaT2_half;
            
            //apply periodic boundary
            //Periodic_Boundary(Temp->coordinate);

            double temp_acc[3] =
            {
                Temp->acceleration[0], Temp->acceleration[1], Temp->acceleration[2]
            };
            Temp->acceleration[0] = -zeta * segment_mass * Temp->velocity[0];
            Temp->acceleration[1] = -zeta * segment_mass * Temp->velocity[1];
            Temp->acceleration[2] = -zeta * segment_mass * Temp->velocity[2];
            
            //calc half time velocity
            Temp->velocity[0] += temp_acc[0] * deltaT_half;
            Temp->velocity[1] += temp_acc[1] * deltaT_half;
            Temp->velocity[2] += temp_acc[2] * deltaT_half;
        }
    }
    kin_step += 0.5 * e_kin * segment_mass;
    vCM[0] *= inv_nParticle;
    vCM[1] *= inv_nParticle;
    vCM[2] *= inv_nParticle;
    vCM_Total += sqrt(vCM[0]*vCM[0] + vCM[1]*vCM[1] + vCM[2]*vCM[2]);
}
