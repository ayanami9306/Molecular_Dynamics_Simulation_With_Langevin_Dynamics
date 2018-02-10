#include "Molecular_Model.hpp"

void Model_Segment::Molecular_Dynamics()
{
    pot_step = 0; kin_step = 0; vCM_Total = 0;
    Set_Params();
    Mol2_File_Write(true);
    MINI_MC(50000);
    Build_NebrList();
    double Rg = 0;
    bond_length = 0;
    for(int i=0;i<=dp_backbone/2; i++) Lp[i] = 0;
    for(int num_Cycle = 1; num_Cycle <= Limit_Cycle; num_Cycle++)
    {
        Single_Step();
        if(num_Cycle>Limit_Cycle/2)
        {
            Measure_Persistence_Length();
            Rg += Measure_Radius_of_Gyration();
            if(!(num_Cycle % step_AVG) && num_Cycle)
            {
                char temp_filename[100];
                sprintf(temp_filename, "%s_%08d", write_filename, num_Cycle);
                Write_State(temp_filename);
                Mol2_File_Write(false);
                Evaluate_Properties();
                for(int i=0;i<=dp_backbone/2; i++) Lp[i] = 0;
                bond_length = 0;
            }
        }
    }
    FILE *fp = fopen(Rg_filename,"a");
    fprintf(fp,"%d,%lf\n",dp_backbone,Rg/(Limit_Cycle/2));
    fclose(fp);
}


