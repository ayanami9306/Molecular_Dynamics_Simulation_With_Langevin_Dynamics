#include "Molecular_Model.hpp"

double Model_Segment::Rand_Standard_Normal_Dist()
{
    static double n2 = 0;
    static bool n2_cached = false;
    //Box Muller Distribution
    if(!n2_cached)
    {
        double x1, x2, abs_x;
        do
        {
            //x = 2.0*rand()/Rand_Max - 1
            x1 = 2.0*rand()*inv_RAND_MAX - 1;
            x2 = 2.0*rand()*inv_RAND_MAX - 1;
            
            abs_x = x1*x1 + x2*x2;
        }
        while (abs_x == 0 || abs_x >= 1.0);
        double factor = sqrt((-2.0*log(abs_x))/abs_x);
        n2 = factor*x2;
        n2_cached = true;
        
        //return factor*x1*stddev + mean
        return factor*x1;
    }
    else
    {
        n2_cached = false;
        //return factor*x2*stddev + mean
        return n2;
    }
}
