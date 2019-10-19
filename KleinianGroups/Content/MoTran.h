//                      
//                      
//               file MoTran.h        
//               
//               Moebius Transformations and SL2(R)
//               
//               Garrett Vance June 19, 2018   
//                      
//                      


#pragma once 

#include <complex>



#define GHV_OPTION_GOLD


#define GHV_OPTION_SPHERICAL_RENDER



#define GHV_OPTION_MESH_POINT  // must use GHV_OPTION_MESH_POINT;


#undef GHV_OPTION_WEIRD_FRACTAL  // active;


#undef GHV_OPTION_HYPERBOLOID  // active; 



namespace VHG 
{


    class MoTran
    {
    public:
        MoTran(
            std::complex<double> p_a,
            std::complex<double> p_b,
            std::complex<double> p_c,
            std::complex<double> p_d
        ) : e_a(p_a), e_b(p_b), e_c(p_c), e_d(p_d) {}

        std::complex<double> MatrixDeterminant(void) { return e_a * e_d - e_b * e_c; }

    public:
        std::complex<double> e_a;
        std::complex<double> e_b;
        std::complex<double> e_c;
        std::complex<double> e_d;
    };



    class EineKleine
    {
    public:
        static unsigned int             k_arc_density;
        static float                    k_sprite_scale;
        static unsigned int             k_depth;
        static double                   k_eps;

    };


    //       
    //   Values for "depth": 
    //   ========================================
    //   formerly : depth = 60;
    //   depth = 4;     // special case to demo animator;
    //   depth = 8;     // tested and working : depth = 8;
    //   depth = 15;    // tested and working if disable render of GeometricPrimitive sphere;
    //   depth = 21, eps = 0001;    // tested and working : depth = 21 will max-out my GPU; 
    //   depth = 26;    // depth = 26 will max-out my GPU and stutter; 
    //               
    //   depth = 31, eps = 0.00005;     ???? 
    //               
    //               

    //   important : depth = 21; eps = 0.0001;  //  GOLD : eps = 0.0001; 

    //   important : depth = 31; eps = 0.00007;













}
//  Closes namespace VHG; 


