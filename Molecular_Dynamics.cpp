#include "Molecular_Model.hpp"

void Model_Segment::Molecular_Dynamics()
{
    pot_step = 0; kin_step = 0; vCM_Total = 0;
    Set_Params();
    Mol2_File_Write(true);
    MINI_MC(100000);
    Build_NebrList();
    Rx_sidechain = 0; Ry_sidechain = 0; Rz_sidechain = 0; Rg_sidechain = 0;
    Rx_backbone = 0; Ry_backbone = 0; Rz_backbone = 0; Rg_backbone = 0;
    Rx_total = 0; Ry_total = 0; Rz_total = 0; Rg_total = 0;
    bond_length = 0;
    for(int i=0;i<=dp_backbone/2; i++) Lp[i] = 0;
    for(int num_Cycle = 1; num_Cycle <= Limit_Cycle; num_Cycle++)
    {
        Single_Step();
        if(num_Cycle>Limit_Cycle/2)
        {
            Measure_Persistence_Length();
            Measure_Radius_of_Gyration();
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
        else
        {
            kin_step = 0;
            pot_step = 0;
        }
    }
    FILE *fp = fopen(Rg_filename,"a");
    fprintf(fp,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",dp_backbone,Rx_sidechain/(Limit_Cycle/2),Ry_sidechain/(Limit_Cycle/2),Rz_sidechain/(Limit_Cycle/2),Rg_sidechain/(Limit_Cycle/2),Rx_backbone/(Limit_Cycle/2),Ry_backbone/(Limit_Cycle/2),Rz_backbone/(Limit_Cycle/2),Rg_backbone/(Limit_Cycle/2),Rx_total/(Limit_Cycle/2),Ry_total/(Limit_Cycle/2),Rz_total/(Limit_Cycle/2),Rg_total/(Limit_Cycle/2));
    fclose(fp);
    //Print_LAMMPS();    
}
