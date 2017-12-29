//
//  main.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 9. 27..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include <iostream>
#include "Molecular_Model.hpp"
#include <sys/time.h>
#include <unistd.h>

Model_Segment System;

int main() {
    struct timeval start, end;
    char Input_Param_File[50] = "Parameters.dat";
    FILE *fp = fopen("Rg.csv","w");
    fclose(fp);
    int num[11] = {50, 100, 150, 200, 600, 1000, 2000, 4000, 6000, 8000, 10000};
    for(int i = 0; i<11; i++)
    { 
        for(int j=0; j<50; j++)
            {
        System.Input_Params(Input_Param_File);
        System.backbone_num = num[i];
    System.Initialize_System(1);
    gettimeofday(&start, NULL);
    System.Molecular_Dynamics();
    gettimeofday(&end, NULL);
    double time_operating = (double)(end.tv_sec)+(double)(end.tv_usec)/1000000.0 - (double)(start.tv_sec)+(double)(start.tv_usec)/1000000.0;
    printf("%f\n",time_operating);
            }
    }
}
