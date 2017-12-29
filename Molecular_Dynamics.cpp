//
//  Molecular_Dynamics.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Molecular_Dynamics()
{
    pot_step = 0; kin_step = 0; vCM_Total = 0;
    Set_Params();
    Mol2_File_Write(true);
    MINI_MC(100000);
    Build_NebrList();
    Rg = 0;
    //PDB_File_Write(true);
    for(int num_Cycle = 0; num_Cycle <= Limit_Cycle; num_Cycle++)
    {
        Single_Step();
        if(num_Cycle>Limit_Cycle/2) Rg += CALCULATE_RADIUS_OF_GYRATION();
        if(!(num_Cycle % step_AVG) && num_Cycle)
        {
            char temp_filename[100];
            sprintf(temp_filename, "%s_%08d", write_filename, num_Cycle);
            //Write_State(temp_filename);
            //PDB_File_Write(false);
            Mol2_File_Write(false);
            Evaluate_Properties();
            //step average property?
            //write data
        }
    }
    FILE *fp = fopen("Rg.csv","a");
    fprintf(fp, "%d,%.4lf\n",backbone_num, Rg/(Limit_Cycle/2));
    fclose(fp);
}
