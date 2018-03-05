#include "Molecular_Model.hpp"

void Model_Segment::Evaluate_Properties()
{
    pot_step *= inv_step_AVG;
    kin_step *= inv_step_AVG;
    vCM_Total *= inv_step_AVG;
    double Hamiltonian_Energy = pot_step + kin_step;
    
    FILE *fp = fopen(dat_filename, "a");
    fprintf(fp, "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n",time_Now, vCM_Total, Hamiltonian_Energy, pot_step, kin_step);
    fclose(fp);
    fp = fopen(Lp_filename, "a");
    for(int i=0;i<=dp_backbone/2; i++) fprintf(fp, "%d,%lf,%lf\n", i, bond_length, Lp[i]);
    fclose(fp);
    
    pot_step = 0; kin_step = 0; vCM_Total = 0;
}

double Model_Segment::Get_V2(int i)
{
    return Segment[i].velocity[0]*Segment[i].velocity[0] + Segment[i].velocity[1]*Segment[i].velocity[1] + Segment     [i].velocity[2]*Segment[i].velocity[2];
}

double Model_Segment::Get_Distance2(int num1, int num2, double * dr)
{
    dr[0] = Segment[num1].coordinate[0] - Segment[num2].coordinate[0];
    dr[1] = Segment[num1].coordinate[1] - Segment[num2].coordinate[1];
    dr[2] = Segment[num1].coordinate[2] - Segment[num2].coordinate[2];
    //Periodic_Length(dr);
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}

void Model_Segment::Vector_Cross(double *V1, double *V2, double *VR)
{
    VR[0] = V1[1]*V2[2] - V1[2]*V2[1];
    VR[1] = V1[2]*V2[0] - V1[0]*V2[2];
    VR[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

void Model_Segment::Normalize_Vector(double *Vector)
{
    double Vec_Size = sqrt(pow(Vector[0], 2.0) + pow(Vector[1], 2.0) + pow(Vector[2], 2.0));
    
    Vector[0] /= Vec_Size;
    Vector[1] /= Vec_Size;
    Vector[2] /= Vec_Size;
}
 
void Model_Segment::Vector_Constant(double *V1, double constant, double*VR)
{
    VR[0] = constant * V1[0];
    VR[1] = constant * V1[1];
    VR[2] = constant * V1[2];
}
