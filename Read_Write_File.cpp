#include "Molecular_Model.hpp"

void Model_Segment::Read_State(char *filename)
{
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d %d %d %d %lf %lf %lf %d %d %lf ", &BOUNDARY_SIZE_X, &BOUNDARY_SIZE_Y, &BOUNDARY_SIZE_Z, &nParticle, &deltaT, &rcut, &time_Now, &step_AVG, &Limit_Cycle, &radius_NebrShell);
    
    for(int i=0; i<nParticle; i++)
    {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &Segment[i].coordinate[0], &Segment[i].coordinate[1], &Segment[i].coordinate[2], &Segment[i].velocity[0], &Segment[i].velocity[1], &Segment[i].velocity[2], &Segment[i].acceleration[0], &Segment[i].acceleration[1], &Segment[i].acceleration[2], &Segment[i].segment_type, &Segment[i].linked_segment_num);
        for(int j=0;j<Segment[i].linked_segment_num; j++)
        {
            fscanf(fp, "%d ", &Segment[i].linked_segment[j]);
        }
    }
}

void Model_Segment::Write_State(char *filename)
{
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d %d %d %d %.5lf %.5lf %.5lf %d %d %.5lf ", BOUNDARY_SIZE_X, BOUNDARY_SIZE_Y, BOUNDARY_SIZE_Z, nParticle, deltaT, rcut, time_Now, step_AVG, Limit_Cycle, radius_NebrShell);
    
    for(int i=0; i<nParticle; i++)
    {
        fprintf(fp, "%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %d %d ", Segment[i].coordinate[0], Segment[i].coordinate[1], Segment[i].coordinate[2], Segment[i].velocity[0], Segment[i].velocity[1], Segment[i].velocity[2], Segment[i].acceleration[0], Segment[i].acceleration[1], Segment[i].acceleration[2], Segment[i].segment_type, Segment[i].linked_segment_num);
        for(int j=0; j<Segment[i].linked_segment_num; j++)
            fprintf(fp, "%d ", Segment[i].linked_segment[j]);
    }
    fclose(fp);
}

//don't use for dendronized polymer(some bugs about VMD or Material Studio)
void Model_Segment::PDB_File_Write(bool is_new_file)
{
    FILE *fp_pdb;
    if(is_new_file)
    {
        fp_pdb = fopen(traj_filename, "w");
    }
    else fp_pdb = fopen(traj_filename, "a");
    //fprintf(fp_pdb, "MODEL%9.4lf\n", time_Now);
    //print HETATM section
    for(int i=0; i<nParticle; i++)
    {
        fprintf(fp_pdb, "ATOM  %5d  C%-2d MOL  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00%12c\n\n",i+1, (i+1)%100, Segment[i].segment_type, Segment[i].coordinate[0], Segment[i].coordinate[1], Segment[i].coordinate[2],'C');
    }
    fprintf(fp_pdb,"TER   %5d\n",nParticle+1);
    //print CONECT section
    for(int i=0; i<nParticle; i++)
    {
        int check = 0;
        for(int j=0; j<Segment[i].linked_segment_num; j++)
        {
            if(!check)
            {
                fprintf(fp_pdb, "CONECT%5d%5d", i+1, Segment[i].linked_segment[j]+1);
                check=1;
            }
            else fprintf(fp_pdb, "%5d", Segment[i].linked_segment[j]+1);
        }
        fprintf(fp_pdb, "\n");
    }
    //print ENDMDL
    fprintf(fp_pdb, "END\n\n");
    fclose(fp_pdb);
}

void Model_Segment::Mol2_File_Write(bool is_new_file)
{
    FILE *fp_mol;
    if(is_new_file)
    {
        fp_mol = fopen(traj_filename, "w");
        
        time_t current_time;
        time(&current_time);
        struct tm t = *localtime(&current_time);
        
        //write file comment
        fprintf(fp_mol, "#\tName : %s\n#\tCreating time : %dyear %dmonth %dday %dhour %dmin %dsec\n\n", traj_filename, t.tm_year+1900, t.tm_mon+1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec);
    }
    else fp_mol = fopen(traj_filename, "a");
    
    //Get Total bond number
    int Total_Bond_Number = 0;
    for(int i=0; i<nParticle; i++) Total_Bond_Number += Segment[i].linked_segment_num;
    Total_Bond_Number /= 2;
    
    fprintf(fp_mol, "\n@<TRIPOS>MOLECULE\n%s\n%d\t%d\t0\t0\t0\nSMALL\nNO_CHARGES\n\n",
            traj_filename, nParticle, Total_Bond_Number);
    
    
    fprintf(fp_mol, "\n@<TRIPOS>ATOM\n");
    for(int i=0; i<nParticle; i++)
        fprintf(fp_mol, "\t%d C%d %g %g %g C.3 %d\n", i+1, i+1, Segment[i].coordinate[0], Segment[i].coordinate[1], Segment[i].coordinate[2], Segment[i].segment_type);
    
    
    fprintf(fp_mol, "\n@<TRIPOS>BOND\n");
    int bond_count = 1;
    for(int i=0; i<nParticle; i++)
        for(int j=0; j<Segment[i].linked_segment_num; j++)
            if(i < Segment[i].linked_segment[j])
                fprintf(fp_mol, "\t%d %d %d 1\n", bond_count++, i+1, Segment[i].linked_segment[j]+1);
    fclose(fp_mol);
}

//input file info to segment about off lattice by mol2 file
void Model_Segment::Mol2_File_Read(char *filename)
{
    int TOTAL_BOND_NUMBER = 0;      //bond number
    char buffer[2001];              //get line string
    
    FILE *fp_mol;
    fp_mol = fopen(filename, "r");
    
    while(!feof(fp_mol))
    {
        fgets(buffer, 2000, fp_mol);
        if(strstr(buffer, "@<TRIPOS>MOLECULE") != NULL)//if molecule section?
        {
            fgets(buffer, 2000, fp_mol);
            fscanf(fp_mol, "%d %d", &nParticle, &TOTAL_BOND_NUMBER);
            //get total segment number and total bond number
        }
        else if(strstr(buffer, "@<TRIPOS>ATOM") != NULL)//if atom section?
        {
            for(int i=0; i<nParticle; i++)
            {
                int atom_num;
                double coordinate[3];
                char atom_type[10];
                
                //get data
                fgets(buffer, 2000, fp_mol);
                sscanf(buffer, "%d %*s %lf %lf %lf %s", &atom_num, &coordinate[0], &coordinate[1], &coordinate[2], atom_type);
                //get coordinate, atom type
                for(int j=0; j<3; j++) Segment[atom_num-1].coordinate[j] = coordinate[j];
                Segment[atom_num-1].linked_segment_num = 0;
            }
        }
        else if(strstr(buffer, "@<TRIPOS>BOND") != NULL)//if bond section?
        {
            for(int i=0; i<TOTAL_BOND_NUMBER; i++)
            {
                int segment1, segment2;
                char bond_type[10];
                
                fgets(buffer, 2000, fp_mol);
                sscanf(buffer, "%*d %d %d %s", &segment1, &segment2, bond_type);
                //get bond data
                Segment[segment1-1].linked_segment[Segment[segment1-1].linked_segment_num++] = segment2-1;
                Segment[segment2-1].linked_segment[Segment[segment2-1].linked_segment_num++] = segment1-1;
            }
        }
    }
    fclose(fp_mol);
}
