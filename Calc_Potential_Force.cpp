#include "Molecular_Model.hpp"

double Model_Segment::Calc_Lennard_Jones_Potential()
{
    double pot_temp = 0;
//#pragma omp parallel for reduction (+:pot_temp)
    for(int i=0; i<num_NebrList; i++)
    {
        int num1 = NebrList_1[i], num2 = NebrList_2[i];
        double pair_distance[3];
        double distance2 = Get_Distance2(num1, num2, pair_distance);
        if(distance2 <= rcut2)
        {
            double inv_distance2 = 1.0 / distance2;
            double inv_distance6 = pow(inv_distance2, 3.0);
            
            //shifted lennard-jones form
            //V_LJ_trunc(r) = V_LJ(r) - V_LJ(r_c) (if r <= r_c),
            //so,f_LJ = 48*epsilon*((sigma/r)^12 - 0.5*(sigma/r)^6)*(vector_r)/r^2
            //f_LJ = -du/(dr)*vector_r/abs(r)
            //else , 0
            double pair_force = 48.0 * epsilon * inv_distance2 * inv_distance6 * (inv_distance6 - 0.5);
            pair_distance[0] *= pair_force;
            pair_distance[1] *= pair_force;
            pair_distance[2] *= pair_force;
            //not use 'for' statement for speed
            
            //#pragma omp atomic
            Segment[num1].acceleration[0] += pair_distance[0];
            //#pragma omp atomic
            Segment[num1].acceleration[1] += pair_distance[1];
            //#pragma omp atomic
            Segment[num1].acceleration[2] += pair_distance[2];
            //#pragma omp atomic
            Segment[num2].acceleration[0] -= pair_distance[0];
            //#pragma omp atomic
            Segment[num2].acceleration[1] -= pair_distance[1];
            //#pragma omp atomic
            Segment[num2].acceleration[2] -= pair_distance[2];
            //calc potential energy : lennard-jones potential
            pot_temp += 4 * epsilon * inv_distance6 * (inv_distance6 - 1.0) - potential_rcut;
        }
    }
    return pot_temp;
}

double Model_Segment::Calc_Bond_Length_Potential_AND_Apply_Langevin()
{
    double pot_temp = 0;
    double random_num[10000];
    for(int i=0;i<nParticle*3;i++) random_num[i] = Rand_Standard_Normal_Dist();
//#pragma omp parallel reduction(+:                                                                                                                    pot_temp)
    {
//#pragma omp for
        for(int i=0; i<nParticle; i++)
        {
            for(int j=0; j<Segment[i].linked_segment_num; j++)
            {
                if(i < Segment[i].linked_segment[j])
                {
                    int num2 = Segment[i].linked_segment[j];
                    double bond_distance[3];
                    double distance = sqrt(Get_Distance2(i, num2, bond_distance));
                    //U_bond = 0.5*k*(r-r0)^2, so
                    //f_bond = -k(r-r0)*r_vector/r
                    if(distance < bond_length_FENE_0)
                    {
                        //FENE potential
                        //U = -0.5*k*R0*ln(1-(r/R0)^2))
                        //F = -k*r/(1-(r/R0)^2)
                        double distance_ratio = distance / bond_length_FENE_0;
                        double bond_force = -k_FENE/(1-pow(distance_ratio, 2.0));
                        bond_distance[0] *= bond_force;
                        bond_distance[1] *= bond_force;
                        bond_distance[2] *= bond_force;
                        //#pragma omp atomic
                        Segment[i].acceleration[0] += bond_distance[0];
                        //#pragma omp atomic
                        Segment[i].acceleration[1] += bond_distance[1];
                        //#pragma omp atomic
                        Segment[i].acceleration[2] += bond_distance[2];
                        //#pragma omp atomic
                        Segment[num2].acceleration[0] -= bond_distance[0];
                        //#pragma omp atomic
                        Segment[num2].acceleration[1] -= bond_distance[1];
                        //#pragma omp atomic
                        Segment[num2].acceleration[2] -= bond_distance[2];
                        
                        pot_temp += log(1-pow(distance_ratio, 2.0));
                    }
                }
            }
            //add random force for speed
            //ma = -m*zeta*v + F + random_force
            int ran_pos = i * 3;
            Segment[i].acceleration[0] += rand_deviation * random_num[ran_pos];
            Segment[i].acceleration[1] += rand_deviation * random_num[ran_pos+1];
            Segment[i].acceleration[2] += rand_deviation * random_num[ran_pos+2];
        }
    }
    return -0.5*k_FENE*pow(bond_length_FENE_0, 2.0)*pot_temp;
}
