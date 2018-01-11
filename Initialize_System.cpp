//
//  Initialize_System.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Initialize_System(int nTypes)
{
    //char mol2_file[50] = "Dendronized.mol2";
    //Mol2_File_Read(mol2_file);
    Initialize_Coordinate(nTypes);
    inv_nParticle = 1/(double)nParticle;
    Initialize_Velocity();
}

void Model_Segment::Initialize_Coordinate(int nTypes)
{
    nParticle = 0;
    //octane, for test
    if(nTypes == 1)
    {
        nParticle = backbone_num;
        double x = BOUNDARY_SIZE_X/2 - (nParticle/2) * (rcut/2), y = BOUNDARY_SIZE_Y/2, z = BOUNDARY_SIZE_Z/2;
        for(int i=0; i<nParticle; i++)
        {
            double torsional_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            double bond_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            Segment[i].coordinate[0] = x + bond_length_0 * sin(bond_angle)*cos(torsional_angle);
            Segment[i].coordinate[1] = y + bond_length_0 * sin(bond_angle)*sin(torsional_angle);
            Segment[i].coordinate[2] = z + bond_length_0 * cos(bond_angle);
            Segment[i].linked_segment_num = 0;
            x = Segment[i].coordinate[0];
            y = Segment[i].coordinate[1];
            z = Segment[i].coordinate[2];
            Segment[i].segment_type = 0;
            if(i != nParticle-1) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i+1;
            if(i != 0) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i-1;
        }
    }
    //dedronized polymer
    else if(nTypes == 2)
    {
        //int backbone_interval = 4;
        int side_generation = 0;
        //int branch_num = 2;
        backbone_num = 20;
        double x = BOUNDARY_SIZE_X/2 - (backbone_num/2) * (bond_length_0/2), y = BOUNDARY_SIZE_Y/2, z= BOUNDARY_SIZE_Z/2;
        for(int i=0;i<backbone_num;i++)
        {
            Segment[i].coordinate[0] = x+bond_length_0 * (1+(rand()*inv_RAND_MAX-0.5)*0.1);
            Segment[i].coordinate[1] = y + ((rand()*inv_RAND_MAX-0.5)*0.1);
            Segment[i].coordinate[2] = z + ((rand()*inv_RAND_MAX-0.5)*0.1);
            x = Segment[i].coordinate[0];
            Segment[i].linked_segment_num = 0;
            Segment[i].segment_type = 0;
            if(i != backbone_num-1) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i+1;
            if(i != 0) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i-1;
            nParticle++;
        }
        for(int i=0; i<backbone_num;i++)
        {
            recursive_branch(i, side_generation);
        }
    }
}

void Model_Segment::recursive_branch(int i, int generation)
{
   
    if(generation)
    {
        double torsional_angle = (double)(rand()%180)*PI/(180.0);
        double bond_angle = (double)(rand()%90)*PI/(180.0);
        Segment[nParticle].coordinate[0] = Segment[i].coordinate[0] + bond_length_0 * sin(bond_angle)*cos(torsional_angle);
        Segment[nParticle].coordinate[1] = Segment[i].coordinate[1] + bond_length_0 * sin(bond_angle)*sin(torsional_angle);
        Segment[nParticle].coordinate[2] = Segment[i].coordinate[2] + bond_length_0 * cos(bond_angle);
        Segment[i].linked_segment[Segment[i].linked_segment_num++] = nParticle;
        Segment[nParticle].segment_type = 1;
        Segment[nParticle].linked_segment[Segment[nParticle].linked_segment_num++] = i;
        nParticle++;
        int parent_node = nParticle -1;
        recursive_branch(parent_node, generation-1);
        recursive_branch(parent_node, generation-1);
    }
    else return;
}

void Model_Segment::Initialize_Velocity()
{
    double vCM[3] = {0, 0, 0};
    
    //initialize velocity
    for(int i=0; i<nParticle; i++)
    {
        double temp_v;
        for(int j=0; j<3; j++)
        {
            temp_v = rand()*inv_RAND_MAX - 0.5;
            Segment[i].velocity[j] = temp_v;
            vCM[j] += temp_v;
        }
    }
    
    // shift and rescale velocities
    // require v_CM(t=0) == 0 & rkT(t=0) == rkT0
    for(int i=0; i<3; i++) vCM[i] *= inv_nParticle;
    
    for(int i=0; i<nParticle; i++)
        for(int j=0; j<3; j++) Segment[i].velocity[j] -= vCM[j];
    =9
    double initial_Ek = 0;
    for(int i=0; i<nParticle; i++)
        for(int j=0; j<3; j++) initial_Ek += pow(Segment[i].velocity[j], 2.0);
    
    //2Ek = <mv^2> = 3kT
    initial_Ek *= (inv_nParticle * segment_mass);
    double multiply_ratio = sqrt(3 * kT_0 / initial_Ek);
    
//#pragma omp parallel for
    for(int i=0; i<nParticle; i++)
        for(int j=0; j<3; j++) Segment[i].velocity[j] *= multiply_ratio;
    
    initial_Ek = 0;
    
//#pragma omp parallel for reduction(+:initial_Ek)
    for(int i=0; i<nParticle; i++)
        for(int j=0; j<3; j++) initial_Ek += pow(Segment[i].velocity[j], 2.0);
    
    initial_Ek *= (0.5 * inv_nParticle * segment_mass);
    printf("%lf\n",initial_Ek);
}
