//
//  Eval_Property.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Evaluate_Properties()
{
    pot_step *= inv_step_AVG;
    kin_step *= inv_step_AVG;
    vCM_Total *= inv_step_AVG;
    double Hamiltonian_Energy = pot_step + kin_step;
    
    FILE *fp = fopen(dat_filename, "a");
    //double persistance_length = Estimate_Lp();
    fprintf(fp, "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n",time_Now, vCM_Total, Hamiltonian_Energy, pot_step, kin_step);
    fclose(fp);
    
    pot_step = 0; kin_step = 0; vCM_Total = 0;
}

double Model_Segment::Get_V2(int i)
{
    return Segment[i].velocity[0]*Segment[i].velocity[0] + Segment[i].velocity[1]*Segment[i].velocity[1] + Segment[i].velocity[2]*Segment[i].velocity[2];
}

double Model_Segment::Get_Distance2(int num1, int num2, double * dr)
{
    dr[0] = Segment[num1].coordinate[0] - Segment[num2].coordinate[0];
    dr[1] = Segment[num1].coordinate[1] - Segment[num2].coordinate[1];
    dr[2] = Segment[num1].coordinate[2] - Segment[num2].coordinate[2];
    //Periodic_Length(dr);
    return dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
}
