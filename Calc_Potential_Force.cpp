#include "Molecular_Model.hpp"

double Model_Segment::Calc_Nonbonded_Potential()
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

double Model_Segment::Calc_Bonded_Potential_AND_Apply_Langevin()
{
    double pot_temp_FENE = 0, pot_temp_Angle=0;
    double random_num[1000000];
    for(int i=0;i<nParticle*3;i++) random_num[i] = Rand_Standard_Normal_Dist();
//#pragma omp parallel reduction(+:                                                                                                                    pot_temp)
    {
//#pragma omp for
        for(int i=0; i<nParticle; i++)
        {
            for(int j=0; j<Segment[i].linked_segment_num; j++)
            {
                int num_prev = Segment[i].linked_segment[j];
                if(i < num_prev)
                {
                    double bond_distance[3];
                    double distance = sqrt(Get_Distance2(i, num_prev, bond_distance));
                    if(distance < bond_length_FENE_0)
                    {
                        //FENE potential
                        //U = -0.5*k*R0^2*ln(1-(r/R0)^2)) (r < R0), else U = inf
                        //F = -k*r/(1-(r/R0)^2 )
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
                        Segment[num_prev].acceleration[0] -= bond_distance[0];
                        //#pragma omp atomic
                        Segment[num_prev].acceleration[1] -= bond_distance[1];
                        //#pragma omp atomic
                        Segment[num_prev].acceleration[2] -= bond_distance[2];
                        
                        pot_temp_FENE += log(1-pow(distance_ratio, 2.0));
                    }
                }
#if APPLY_BOND_ANGLE
                if(Segment[i].segment_type)
                {
                    for(int k=j+1; k<Segment[i].linked_segment_num; k++)
                    {
                        int num_next = Segment[i].linked_segment[k];
                        if(Segment[num_next].segment_type || Segment[num_prev].segment_type)
                        {
                            double bond_prev[3], bond_next[3];
                            //ba vector
                            double distance_prev = sqrt(Get_Distance2(num_prev, i, bond_prev));
                            //bc vector
                            double distance_next = sqrt(Get_Distance2(num_next, i, bond_next));
                            
                            double Cos_Value = scalar_multiply(bond_prev, bond_next) / (distance_prev * distance_next);
                            double bond_angle;
                            if(Cos_Value >= 1) bond_angle = 0;
                            else if(Cos_Value <= -1) bond_angle = PI;
                            else bond_angle = acos(Cos_Value);
                            double Force_A[3], Force_B[3], Force_C[3];
                            
                            double dU = 2*bond_angle_k*(bond_angle - bond_angle_0);
                            
                            double BA_BC[3];
                            Vector_Cross(bond_prev, bond_next, BA_BC);
                            Vector_Cross(bond_prev, BA_BC, Force_A);
                            Normalize_Vector(Force_A);
                            Vector_Constant(Force_A, -dU/distance_prev, Force_A);
                            Vector_Constant(bond_next, -1, bond_next);
                            Vector_Cross(bond_next, BA_BC, Force_C);
                            Normalize_Vector(Force_C);
                            Vector_Constant(Force_C, -dU/distance_next, Force_C);
                            
                            Force_B[0] = -Force_A[0] -Force_C[0];
                            Force_B[1] = -Force_A[1] -Force_C[1];
                            Force_B[2] = -Force_A[2] -Force_C[2];
                            
                            Segment[num_prev].acceleration[0] += Force_A[0];
                            Segment[num_prev].acceleration[1] += Force_A[1];
                            Segment[num_prev].acceleration[2] += Force_A[2];
                            Segment[i].acceleration[0] += Force_B[0];
                            Segment[i].acceleration[1] += Force_B[1];
                            Segment[i].acceleration[2] += Force_B[2];
                            Segment[num_next].acceleration[0] += Force_C[0];
                            Segment[num_next].acceleration[1] += Force_C[1];
                            Segment[num_next].acceleration[2] += Force_C[2];
                            
                            pot_temp_Angle += pow(bond_angle-bond_angle_0, 2.0);
                        }
                    }
                }
#endif
            }
            //add random force for speed
            //ma = -m*zeta*v + F + random_force
            int ran_pos = i * 3;
            Segment[i].acceleration[0] += rand_deviation * random_num[ran_pos];
            Segment[i].acceleration[1] += rand_deviation * random_num[ran_pos+1];
            Segment[i].acceleration[2] += rand_deviation * random_num[ran_pos+2];
        }
    }
    return -0.5*k_FENE*pow(bond_length_FENE_0, 2.0)*pot_temp_FENE + bond_angle_k * pot_temp_Angle;
}
