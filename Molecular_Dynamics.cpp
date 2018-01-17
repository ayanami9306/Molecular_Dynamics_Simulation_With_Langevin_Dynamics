#include "Molecular_Model.hpp"

void Model_Segment::Molecular_Dynamics()
{
    pot_step = 0; kin_step = 0; vCM_Total = 0;
    Set_Params();
    Mol2_File_Write(true);
    MINI_MC(100000);
    Build_NebrList();
    for(int i=0;i<dp_backbone;i++) Lp[i] = 0;
    for(int num_Cycle = 0; num_Cycle <= Limit_Cycle; num_Cycle++)
    {
        Single_Step();
        Estimate_Lp(Lp);
        if(!(num_Cycle % step_AVG) && num_Cycle)
        {
            char temp_filename[100];
            sprintf(temp_filename, "%s_%08d", write_filename, num_Cycle);
            Write_State(temp_filename);
            Mol2_File_Write(false);
            Evaluate_Properties();
            for(int i=0;i<dp_backbone;i++) Lp[i] = 0;
        }
    }
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
