#include "Molecular_Model.hpp"
#include <sys/time.h>
#include <unistd.h>

Model_Segment::Model_Segment()
{
    inv_RAND_MAX = 1.0 / RAND_MAX;
    srand((unsigned int)time(NULL));
    unsigned int seed[16];
    for(int i=0;i<16;i++) seed[i] = (rand()+getpid())%RAND_MAX;
    InitWELLRNG512a(seed);
}

Model_Segment::~Model_Segment()
{
    
}

