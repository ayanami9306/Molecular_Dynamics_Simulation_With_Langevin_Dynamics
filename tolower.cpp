#include "Molecular_Model.hpp"

void Model_Segment::tolower(char *data)
{
    for(int i=0; i<strlen(data); i++)
    {
        if(data[i] <= 90 && data[i] >= 65) data[i] += 32;
    }
}
