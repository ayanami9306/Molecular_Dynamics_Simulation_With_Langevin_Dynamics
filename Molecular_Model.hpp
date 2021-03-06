/***********************************************
 
 LENGTH IS PROPORTIONAL TO SIGMA
 ENERGY IS PROPORTIONAL TO EPSILON
 
 **********************************************/
#ifndef Molecular_Model_hpp
#define Molecular_Model_hpp

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "omp.h"
#include "WELL512a.hpp"

#define NUM_MOLECULE 100000
#define MAX_CORENUM 30
#define PI 3.141592

#define scalar_multiply(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define APPLY_BOND_ANGLE 1
using namespace std;

#endif /* Molecular_Model_hpp */

typedef struct Node
{
    double coordinate[3];
    double velocity[3];
    double acceleration[3];
    int linked_segment_num;
    int linked_segment[10];
    int segment_type;
}Node;

class Model_Segment
{
    //variables
private:
    Node Segment[NUM_MOLECULE];
    int nParticle, step_AVG, Limit_Cycle, BOUNDARY_SIZE_X, BOUNDARY_SIZE_Y, BOUNDARY_SIZE_Z, pragma_core, dp_dendron, generation_dendron, number_branch, num_NebrList, space_sidechain;
    double deltaT, rcut, segment_mass, zeta, rand_deviation, epsilon, rMax, inv_nParticle, inv_RAND_MAX;
    double deltaT_half, deltaT2_half, rcut2, inv_step_AVG, potential_rcut, kT_0, k_FENE, vvMax, vCM_Total, pot_step, kin_step, time_Now, inv_Rand_Max, accu_movement, radius_NebrShell, bond_length_FENE_0, Lp[10000], inv_segment_mass, bond_length, bond_angle_0, bond_angle_k;
    double Rx_backbone, Ry_backbone, Rz_backbone, Rg_backbone, Rx_sidechain, Ry_sidechain, Rz_sidechain, Rg_sidechain, Rx_total, Ry_total, Rz_total, Rg_total;
    unsigned int NebrList_1[100000000], NebrList_2[100000000];
    bool recon_Nebrlist;
    char write_filename[100], dat_filename[100], traj_filename[100], Lp_filename[100], Rg_filename[100];
    
    //methods
public:
    Model_Segment();
    ~Model_Segment();
    void Initialize_System(int nTypes);
    void Input_Params(char *filename);
    void Molecular_Dynamics();
    void Read_State(char *filename);
    void Write_State(char *filename);
    double Rand_Standard_Normal_Dist();
    int dp_backbone;
    
    
private:
    void INITIALIZE_RANDOM_SEED(int thread_num);
    void Initialize_Coordinate(int nTypes);
    void Initialize_Velocity();
    void tolower(char *data);
    void Single_Step();
    void Set_Params();
    void Compute_Forces();
    double Get_V2(int i);
    double Get_Distance2(int num1, int num2, double * dr);
    void Build_NebrList();
    void Periodic_Length (double * data);
    void Periodic_Boundary (double * data);
    void Velocity_Verlet_Step();
    void Velocity_Verlet_After_Step();
    double Calc_Nonbonded_Potential();
    double Calc_Bonded_Potential_AND_Apply_Langevin();
    void recursive_branch(int parent, int generation);
    void Evaluate_Properties();
    void PDB_File_Write(bool is_new_file);
    void Mol2_File_Write(bool is_new_file);
    void Mol2_File_Read(char *filename);
    void MINI_MC(int nCycle);
    void Measure_Persistence_Length();
    void Measure_Radius_of_Gyration();
    double Measure_Rotational_Viscosity();
    void Vector_Cross(double *V1, double *V2, double *VR);
    void Vector_Constant(double *V1, double constant, double*VR);
    void Normalize_Vector(double *Vector);
    void Print_LAMMPS();
    void MC_NEIGHBORLIST();
};
