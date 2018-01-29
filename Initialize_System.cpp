#include "Molecular_Model.hpp"

void Model_Segment::Initialize_System(int nTypes)
{
    Initialize_Coordinate(nTypes);
    inv_nParticle = 1/(double)nParticle;
    Initialize_Velocity();
}

void Model_Segment::Initialize_Coordinate(int nTypes)
{
    //dedronized polymer
    if(nTypes == 1)
    {
        double x = 0, y = 0, z= 0;
        for(int i=0;i<dp_backbone;i++)
        {
            double torsional_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            double bond_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            Segment[i].coordinate[0] = x + distance_FENE_0*0.67 * sin(bond_angle)*cos(torsional_angle);
            Segment[i].coordinate[1] = y + distance_FENE_0*0.67 * sin(bond_angle)*sin(torsional_angle);
            Segment[i].coordinate[2] = z + distance_FENE_0*0.67 * cos(bond_angle);
            Segment[i].linked_segment_num = 0;
            Segment[i].segment_type = 0;
            x = Segment[i].coordinate[0];
            y = Segment[i].coordinate[1];
            z = Segment[i].coordinate[2];
            if(i != dp_backbone-1) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i+1;
            if(i != 0) Segment[i].linked_segment[Segment[i].linked_segment_num++] = i-1;
            nParticle++;
        }
        for(int i=0; i<dp_backbone;i+=space_sidechain) recursive_branch(i, generation_dendron);
    }
}

void Model_Segment::recursive_branch(int parent, int generation)
{
    if(generation)
    {
        int previous_node = parent;
        for(int i=0; i<dp_dendron; i++)
        {
            double torsional_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            double bond_angle = ((double)rand()*inv_RAND_MAX)*2*PI;
            Segment[nParticle].coordinate[0] = Segment[previous_node].coordinate[0] + distance_FENE_0*0.67 * sin(bond_angle)*cos(torsional_angle);
            Segment[nParticle].coordinate[1] = Segment[previous_node].coordinate[1] + distance_FENE_0*0.67 * sin(bond_angle)*sin(torsional_angle);
            Segment[nParticle].coordinate[2] = Segment[previous_node].coordinate[2] + distance_FENE_0*0.67 * cos(bond_angle);
            Segment[nParticle].segment_type = 1;
            Segment[previous_node].linked_segment[Segment[previous_node].linked_segment_num++] = nParticle;
            Segment[nParticle].linked_segment[Segment[nParticle].linked_segment_num++] = previous_node;
            previous_node = nParticle;
            nParticle++;
        }
        for(int i=0; i<number_branch; i++) recursive_branch(previous_node, generation-1);
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
