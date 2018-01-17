#include "Molecular_Model.hpp"

void Model_Segment::Build_NebrList()
{
    num_NebrList = 0;
    double r2_Nebr = (rcut + radius_NebrShell) * (rcut + radius_NebrShell);
    for(int i=0; i<nParticle-1; i++)
    {
        for(int j=i+1; j<nParticle; j++)
        {
            double temp[3];
            double dist2 = Get_Distance2(i, j, temp);
            if(dist2 < r2_Nebr)
            {
                NebrList_1[num_NebrList] = i;
                NebrList_2[num_NebrList++] = j;
            }
        }
    }
}

