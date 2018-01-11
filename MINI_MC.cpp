//
//  MINI_MC.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 23..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"
void Model_Segment::MINI_MC(int nCycle)
{
    for(int mccount=0; mccount<nCycle/2; mccount++)
    {
        for(int i=0; i<nParticle; i++)
        {
            int rSegment = (int)(((double)rand()/RAND_MAX)*nParticle);
            double theta = 2*PI*((double)rand()/RAND_MAX), phi = 2*PI*((double)rand()/RAND_MAX), radius = rMax*((double)rand()/RAND_MAX);
            double temp_coordinate[3] =
            {
                Segment[rSegment].coordinate[0],
                Segment[rSegment].coordinate[1],
                Segment[rSegment].coordinate[2]
            };
            Segment[rSegment].coordinate[0] += radius * sin(theta) * cos(phi);
            Segment[rSegment].coordinate[1] += radius * sin(theta) * sin(phi);
            Segment[rSegment].coordinate[2] += radius * cos(theta);
            
            bool check = true;
            for(int j=0; j<Segment[rSegment].linked_segment_num;j++)
            {
                double pair_distnace[3];
                double num2 = Segment[rSegment].linked_segment[j];
                double distance = sqrt(Get_Distance2(rSegment, num2, pair_distnace));
                if(distance > 1.02*bond_length_0 || distance < 0.98*bond_length_0) check = false;
            }
            if(!check)
            {
                Segment[rSegment].coordinate[0] = temp_coordinate[0];
                Segment[rSegment].coordinate[1] = temp_coordinate[1];
                Segment[rSegment].coordinate[2] = temp_coordinate[2];
            }
        }
        //if(!(mccount%step_AVG)) Mol2_File_Write(false);
    }
    for(int mccount=0; mccount<nCycle/2; mccount++)
    {
        for(int i=0; i<nParticle; i++)
        {
            double pot_pre = 0, pot_after = 0;
            int rSegment = (int)(((double)rand()/RAND_MAX)*nParticle);
            for(int j=0;j<nParticle;j++)
            {
                if(j!=rSegment)
                {
                    double pair_distance[3];
                    double distance2 = Get_Distance2(rSegment, j, pair_distance);
                    if(distance2 <= rcut2)
                    {
                        double inv_distance2 = 1.0 / distance2;
                        double inv_distance6 = pow(inv_distance2, 3.0);
                        pot_pre += 4 * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
                    }
                }
            }
            for(int j=0; j<Segment[rSegment].linked_segment_num;j++)
            {
                int num2 = Segment[rSegment].linked_segment[j];
                double bond_distance[3];
                double distance = sqrt(Get_Distance2(rSegment, num2, bond_distance));
                pot_pre += 0.5*k_spring*pow((distance-bond_length_0), 2.0);
            }
            double theta = 2*PI*((double)rand()/RAND_MAX), phi = 2*PI*((double)rand()/RAND_MAX), radius = rMax*((double)rand()/RAND_MAX);
            double temp_coordinate[3] =
            {
                Segment[rSegment].coordinate[0],
                Segment[rSegment].coordinate[1],
                Segment[rSegment].coordinate[2]
            };
            Segment[rSegment].coordinate[0] += radius * sin(theta) * cos(phi);
            Segment[rSegment].coordinate[1] += radius * sin(theta) * sin(phi);
            Segment[rSegment].coordinate[2] += radius * cos(theta);
            
            for(int j=0;j<nParticle;j++)
            {
                if(j!=rSegment)
                {
                    double pair_distance[3];
                    double distance2 = Get_Distance2(rSegment, j, pair_distance);
                    if(distance2 <= rcut2)
                    {
                        double inv_distance2 = 1.0 / distance2;
                        double inv_distance6 = pow(inv_distance2, 3.0);
                        pot_after += 4 *  inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
                    }
                }
            }
            for(int j=0; j<Segment[rSegment].linked_segment_num;j++)
            {
                int num2 = Segment[rSegment].linked_segment[j];
                double bond_distance[3];
                double distance = sqrt(Get_Distance2(rSegment, num2, bond_distance));
                pot_after += 0.5*k_spring*pow((distance-bond_length_0), 2.0);
            }
            
            double delta_pot = pot_after - pot_pre;
            if(exp(-delta_pot) < (double)rand()/RAND_MAX)
            {
                Segment[rSegment].coordinate[0] = temp_coordinate[0];
                Segment[rSegment].coordinate[1] = temp_coordinate[1];
                Segment[rSegment].coordinate[2] = temp_coordinate[2];
            }
        }
        //if(!(mccount%step_AVG)) Mol2_File_Write(false);
    }
}
