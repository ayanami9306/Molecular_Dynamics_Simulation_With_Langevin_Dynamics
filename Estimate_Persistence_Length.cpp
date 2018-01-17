#include "Molecular_Model.hpp"

#define scalar_multiply(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

void Model_Segment::Estimate_Lp(double * Lp)
{
    double bond_length = 0;
    for(int i=0; i<dp_backbone-1; i++)
    {
        double pair_temp[3];
        bond_length += sqrt(Get_Distance2(i, i+1, pair_temp));
    }
    bond_length /= (double)(dp_backbone-1);
    
    double av_pair_scalar[2000] = {0, };
    for(int s = 1; s <= dp_backbone/2; s++)
    {
        for(int i = 1; i<dp_backbone  - s; i++)
        {
              double pair_distance_pre[3];
            double pair_distance_cur[3];
            
            Get_Distance2(i, i-1, pair_distance_pre);
            Get_Distance2(i+s, i+s-1, pair_distance_cur);
            
            av_pair_scalar[s] += scalar_multiply(pair_distance_pre, pair_distance_cur);
        }
        Lp[s] = - bond_length / log(av_pair_scalar[s]/(bond_length*bond_length)) ;
    }
}
