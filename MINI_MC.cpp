#include "Molecular_Model.hpp"
int MC_LIST[20000][10000];
int num_MCLIST[20000];
void Model_Segment::MC_NEIGHBORLIST()
{
    double r2_Nebr = pow(rcut + radius_NebrShell, 2.0);
    for(int i=0; i<nParticle;i++) num_MCLIST[i] = 0;
    for(int i=0; i<nParticle-1; i++)
    {
        for(int j=i+1; j<nParticle; j++)
        {
            double temp[3];
            double dist2 = Get_Distance2(i, j, temp);
            if(dist2 < r2_Nebr)
            {
                MC_LIST[i][num_MCLIST[i]++] = j;
                MC_LIST[j][num_MCLIST[j]++] = i;
            }
        }
    }
}
void Model_Segment::MINI_MC(int nCycle)
{
    MC_NEIGHBORLIST();
    int BUILD_COUNT = (int)(radius_NebrShell/(rMax));
    for(int mccount=0; mccount<nCycle; mccount++)
    {
        if(mccount % BUILD_COUNT) MC_NEIGHBORLIST();
        for(int i=0; i<nParticle; i++)
        {
            //calc before potential
            double pot_pre = 0, pot_after = 0;
            int rSegment = (int)(((double)rand()/RAND_MAX)*nParticle);
            for(int j=0;j<num_MCLIST[rSegment];j++)
            {
                int num2 = MC_LIST[rSegment][j];
                double pair_distance[3];
                double distance2 = Get_Distance2(rSegment, num2, pair_distance);
                if(distance2 <= rcut2)
                {
                    double inv_distance2 = 1.0 / distance2;
                    double inv_distance6 = pow(inv_distance2, 3.0);
                    pot_pre += 4 * epsilon * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
                }
            }
            for(int j=0; j<Segment[rSegment].linked_segment_num;j++)
            {
                int num_prev = Segment[rSegment].linked_segment[j];
                double bond_prev[3], bond_next[3];
                //ba vector
                double distance_prev = sqrt(Get_Distance2(num_prev, rSegment, bond_prev));
                if(distance_prev < bond_length_FENE_0) pot_pre += -0.5*k_FENE*pow(bond_length_FENE_0,2.0)*log(1-pow(distance_prev/bond_length_FENE_0, 2.0));
                else pot_pre+=10E99;
#if APPLY_BOND_ANGLE
                if(Segment[rSegment].segment_type)
                {
                    for(int k=j+1; k<Segment[rSegment].linked_segment_num; k++)
                    {
                        int num_next = Segment[rSegment].linked_segment[k];
                        if(Segment[num_next].segment_type || Segment[num_prev].segment_type)
                        {
                            //bc vector
                            double distance_next = sqrt(Get_Distance2(num_next, rSegment, bond_next));
                            
                            double Cos_Value = scalar_multiply(bond_prev, bond_next) / (distance_prev * distance_next);
                            double bond_angle;
                            if(Cos_Value >= 1) bond_angle = 0;
                            else if(Cos_Value <= -1) bond_angle = PI;
                            else bond_angle = acos(Cos_Value);
                            pot_pre += bond_angle_k * pow(bond_angle-bond_angle_0, 2.0);
                        }
                    }
                }
#endif
            }
            
            //calc after potential
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
            
            for(int j=0;j<num_MCLIST[rSegment];j++)
            {
                int num2 = MC_LIST[rSegment][j];
                double pair_distance[3];
                double distance2 = Get_Distance2(rSegment, num2, pair_distance);
                if(distance2 <= rcut2)
                {
                    double inv_distance2 = 1.0 / distance2;
                    double inv_distance6 = pow(inv_distance2, 3.0);
                    pot_after += 4 * epsilon * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
                }
            }
            for(int j=0; j<Segment[rSegment].linked_segment_num;j++)
            {
                int num_prev = Segment[rSegment].linked_segment[j];
                double bond_prev[3], bond_next[3];
                //ba vector
                double distance_prev = sqrt(Get_Distance2(num_prev, rSegment, bond_prev));
                if(distance_prev < bond_length_FENE_0) pot_after += -0.5*k_FENE*pow(bond_length_FENE_0,2.0)*log(1-pow(distance_prev/bond_length_FENE_0, 2.0));
                else pot_after+=10E99;
#if APPLY_BOND_ANGLE
                if(Segment[rSegment].segment_type)
                {
                    for(int k=j+1; k<Segment[rSegment].linked_segment_num; k++)
                    {
                        int num_next = Segment[rSegment].linked_segment[k];
                        if(Segment[num_next].segment_type || Segment[num_prev].segment_type)
                        {
                            //bc vector
                            double distance_next = sqrt(Get_Distance2(num_next, rSegment, bond_next));
                            
                            double Cos_Value = scalar_multiply(bond_prev, bond_next) / (distance_prev * distance_next);
                            double bond_angle;
                            if(Cos_Value >= 1) bond_angle = 0;
                            else if(Cos_Value <= -1) bond_angle = PI;
                            else bond_angle = acos(Cos_Value);
                            pot_after += bond_angle_k * pow(bond_angle-bond_angle_0, 2.0);
                        }
                    }
                }
#endif
            }
            
            double delta_pot = pot_after - pot_pre;
            if(exp(-delta_pot) < (double)rand()/RAND_MAX)
            {
                Segment[rSegment].coordinate[0] = temp_coordinate[0];
                Segment[rSegment].coordinate[1] = temp_coordinate[1];
                Segment[rSegment].coordinate[2] = temp_coordinate[2];
            }
        }
        if(mccount%step_AVG == 0) Mol2_File_Write(false);
    }
}
