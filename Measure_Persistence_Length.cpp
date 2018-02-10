#include "Molecular_Model.hpp"

#define scalar_multiply(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

void Model_Segment::Measure_Persistence_Length()
{
    double pair_distance_pre[3];
    double pair_distance_cur[3];
    double bond_length_pre, bond_length_cur;
    for(int i=0; i<dp_backbone-1; i++) bond_length += sqrt(Get_Distance2(i, i+1, pair_distance_pre))/((dp_backbone-1)*step_AVG);
    for(int s = 0; s <= dp_backbone/5; s++)
    {
        for(int i=1;i<dp_backbone-s;i++)
        {
            bond_length_pre = Get_Distance2(i, i-1, pair_distance_pre);
            bond_length_cur = Get_Distance2(i+s, i+s-1, pair_distance_cur);
            Lp[s] += scalar_multiply(pair_distance_pre, pair_distance_cur)/(sqrt(bond_length_pre*bond_length_cur)*(dp_backbone-s-1)*step_AVG);
        }
    }
}
