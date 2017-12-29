//
//  tolower.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::tolower(char *data)
{
    for(int i=0; i<strlen(data); i++)
    {
        if(data[i] <= 90 && data[i] >= 65) data[i] += 32;
    }
}
