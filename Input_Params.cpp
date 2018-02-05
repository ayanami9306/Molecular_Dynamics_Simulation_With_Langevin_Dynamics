#include "Molecular_Model.hpp"

void Model_Segment::Input_Params(char *filename)
{
    FILE *fp = fopen(filename, "r");
    char property_value[50];
    char property_name[50];
    while(fscanf(fp, "%s %s\n", property_name, property_value) > 0)
    {
        tolower(property_name);
        if(!strcmp(property_name, "deltat")) deltaT = atof(property_value);
        else if(!strcmp(property_name, "rcut")) rcut = atof(property_value);
        else if(!strcmp(property_name, "kt_0")) kT_0 = atof(property_value);
        else if(!strcmp(property_name, "step_avg")) step_AVG = atoi(property_value);
        else if(!strcmp(property_name, "segment_mass")) segment_mass = atof(property_value);
        else if(!strcmp(property_name, "limit_cycle")) Limit_Cycle = atof(property_value);
        else if(!strcmp(property_name, "rcut")) rcut = atof(property_value);
        else if(!strcmp(property_name, "zeta")) zeta = atof(property_value);
        else if(!strcmp(property_name, "bond_length_fene_0")) bond_length_FENE_0 = atof(property_value);
        else if(!strcmp(property_name, "write_filename")) strcpy(write_filename, property_value);
        else if(!strcmp(property_name, "time_now")) time_Now = atof(property_value);
        else if(!strcmp(property_name, "radius_nebrshell")) radius_NebrShell = atof(property_value);
        else if(!strcmp(property_name, "boundaryx")) BOUNDARY_SIZE_X = atof(property_value);
        else if(!strcmp(property_name, "boundaryy")) BOUNDARY_SIZE_Y = atof(property_value);
        else if(!strcmp(property_name, "boundaryz")) BOUNDARY_SIZE_Z = atof(property_value);
        else if(!strcmp(property_name, "corenum")) pragma_core = atoi(property_value);
        else if(!strcmp(property_name, "dp_backbone")) dp_backbone = atoi(property_value);
        else if(!strcmp(property_name, "dp_dendron")) dp_dendron = atoi(property_value);
        else if(!strcmp(property_name, "generation_dendron")) generation_dendron = atoi(property_value);
        else if(!strcmp(property_name, "number_branch")) number_branch   = atoi(property_value);
        else if(!strcmp(property_name, "space_sidechain")) space_sidechain   = atoi(property_value);
        else if(!strcmp(property_name, "epsilon")) epsilon = atof(property_value);
    }
}

void Model_Segment::Set_Params()
{
    //omp_set_num_threads(pragma_core);
    
    accu_movement = 0;
    Build_NebrList();
    recon_Nebrlist = false;
    deltaT_half = 0.5*deltaT;
    deltaT2_half = deltaT_half*deltaT;
    rcut2 = rcut*rcut; //rcut unit : sigma
    inv_nParticle = 1/(double)nParticle;
    //kT_0 = epsilon * k_b * T
    rand_deviation = sqrt(2*kT_0*segment_mass*zeta/deltaT);
    rMax = bond_length_FENE_0*0.67*0.03;
    double inverse_rcut6 = 1.0 / pow(rcut2, 3.0);
    potential_rcut = 4 * epsilon * inverse_rcut6 * (inverse_rcut6 - 1.0);
    k_FENE = kT_0/1.2*30;
    inv_segment_mass = 1.0 / (double)segment_mass;
    
    inv_step_AVG = 1.0 / (double)step_AVG;
    
    sprintf(dat_filename, "%s_DAT.csv", write_filename);
    sprintf(traj_filename, "%s.mol2", write_filename);
    sprintf(Lp_filename, "%s_Lp.csv", write_filename);
    if(time_Now == 0)
    {
        FILE *fp = fopen(dat_filename, "w");
        fprintf(fp, "TIME,vCM,HAM,POT,KIN\n");
        fclose(fp);
        fp = fopen(Lp_filename, "w");
        fprintf(fp, "s,value\n");
        fclose(fp);
    }
    
    //initialize force
    for(int i=0;i<nParticle;i++)
    {
        Segment[i].acceleration[0] = 0;
        Segment[i].acceleration[1] = 0;
        Segment[i].acceleration[2] = 0;
    }
}
