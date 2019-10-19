

#include "pch.h"

#include "..\Common\DirectXHelper.h"
#include "gv_Scene3D.h"
#include "MoTran.h"
#include "cx.h"





#ifdef GHV_OPTION_GOLD
unsigned int        VHG::EineKleine::k_depth = 31;     //  GOLD : k_depth = 31 and k_eps = 0.00007;
double              VHG::EineKleine::k_eps = 0.00007;  //  GOLD : k_eps = 0.00007 and k_depth = 31;  
unsigned int        VHG::EineKleine::k_arc_density = 512;       //  k_arc_density = 512 for cameo;  
float               VHG::EineKleine::k_sprite_scale = 0.002f;   //  k_sprite_scale = 0.002f; 

#else
//  experimental : 
unsigned int        VHG::EineKleine::k_depth = 80;    //  
// double              VHG::EineKleine::k_eps = 0.00001;  //  too small for k_eps;
// double              VHG::EineKleine::k_eps = 0.00005;  //  
   double              VHG::EineKleine::k_eps = 0.00007;  //  
// double              VHG::EineKleine::k_eps = 0.0001;  //  largest viable k_eps value = 0.0001; 
unsigned int        VHG::EineKleine::k_arc_density = 512;       //  k_arc_density = 8;  

// float               VHG::EineKleine::k_sprite_scale = 0.002f;   //  k_sprite_scale = 0.002f; 
float               VHG::EineKleine::k_sprite_scale = 0.004f;   //  k_sprite_scale = 0.002f; 
#endif




#define MAXC  300000 	/* Circle stack */
#define MAXG  1000	/* Generator stack */
static double  huge = 100.00; 
static double agree = AGREE;  
static int nfinal = 0; 
static int  nstack = 0;
static size_t circle_bytes;
static double circle_rad_hwm = 0.00; 

VHG::VHG_circle final[MAXC]; 
VHG::VHG_circle stack[MAXC];


VHG::VHG_matrix   g0_generators[MAXG];  //  formerly named "gens"; 
static int  ngens = 0; 



using namespace KG;
using namespace DirectX;
using namespace Windows::Foundation;



int push_circle(VHG::VHG_circle p_a_circle); // forward declaration; 




void evolve_a_circle(VHG::VHG_circle c)
{
    VHG::VHG_circle d;
    int i;

    //     Loop over all the matrices in "g0_generators", 
    //     and call  image_circle() using that current loop's matrix: 
    //   
    //     For the "Limit Sets of Kleinian Groups" aka "cusp", 
    //     the g0_generators array has cardinality four (4); 
    //   

    for (i = 0; i < ngens; i++)
    {
        //     The image_circle() method takes a circle 
        //     and a matrix from the "g0_generators[]" array. 


        d = image_circle(g0_generators[i], c);

        push_circle(d);
    }
}




void evolve_stack()
{
    int nold; nold = nstack;


    for (int i = 0; i < nold; i++) 
    { 
        evolve_a_circle(stack[i]); 
    }
}





/* Cut circle of slopes at an irrational number such as PI */

double slope(double x)
{
    x = x + PI;
    while (x < 0)    x = x + 1.0;
    while (x >= 1.0) x = x - 1.0;
    return(x);
}




// orig : int compare_circles(circle * c, circle * d)

int compare_circles(const void * pc, const void * pd)
{
    const VHG::VHG_circle * c = (const VHG::VHG_circle *)pc;
    const VHG::VHG_circle * d = (const VHG::VHG_circle *)pd;


    double dr, dx, dy;
    double dblsmall;

    /* Lines precede circles */

    if (c->r <= 0 && d->r > 0) return(-1);


    if (d->r <= 0 && c->r > 0) return(1);


    dblsmall = VHG::EineKleine::k_eps / agree;


    dx = c->c.x - d->c.x;
    if (dx < -dblsmall) return(-1);
    if (dx > dblsmall) return(1);

    dy = c->c.y - d->c.y;
    if (dy < -dblsmall) return(-1);
    if (dy > dblsmall) return(1);

    if (c->r <= 0) 		/* Lines */
    {
        dr = slope(c->r) - slope(d->r);
        if (dr < -PARAFUZZ) return(-1);
        if (dr >  PARAFUZZ) return(1);
        return(0);
    }
    else {		/* Circles */
        dr = c->r - d->r;
        if (dr < -dblsmall) return(-1);
        if (dr > dblsmall) return(1);
        return(0);
    }
}





void sort_circles(VHG::VHG_circle stack[], int nstack)
{
    qsort(stack, nstack, circle_bytes, compare_circles);
}





/* Sort a list of circles and eliminate dups */

void unique_circles(VHG::VHG_circle stack[], int * nstack)
{
    int i, n;

    if (*nstack<2) return;
    sort_circles(stack, *nstack);
    i = n = 1;
    while (i<*nstack)
    {
        if (0 == compare_circles(&stack[i], &stack[i - 1]))	i++;
        else { stack[n] = stack[i]; n++; i++; }
    }
    *nstack = n;
}






/* Move circles from stack to final, skipping dups */
/* Returns 1 if run out of room */

int merge_circles()
{
    int c, i, ifinal, istack, nold;

    unique_circles(stack, &nstack);

    if (nstack + nfinal > MAXC) return(1);



    ifinal = istack = 0;
    nold = nfinal;
    while (ifinal < nold && istack < nstack)
    {
        c = compare_circles(&final[ifinal], &stack[istack]);
        if (c<0) ifinal++;
        if (c == 0) istack++;
        if (c>0)
        {
            final[nfinal] = stack[istack];
            nfinal++;
            istack++;
        }
    }
    while (istack < nstack)
    {
        final[nfinal] = stack[istack];
        nfinal++;
        istack++;
    }
    for (i = nold; i<nfinal; i++)
        stack[i - nold] = final[i];
    nstack = nfinal - nold;
    sort_circles(final, nfinal);
    return(0);
}




//   end qsort mystery 








int push_circle(VHG::VHG_circle p_a_circle)
{
    //  Return code: 
    //            
    //   0) OK 
    //   1) too small or far 
    //   2) full


    if (p_a_circle.r < VHG::EineKleine::k_eps && p_a_circle.r > 0) return(1);

    if (p_a_circle.r > huge || cx_abs(p_a_circle.c) > huge) return(1);

    if (nstack == MAXC) return(2);

    stack[nstack] = p_a_circle;

    nstack++;

    return(0);
}
//  Closes push_circle(); 






    





void gv_moebius_to_matrix(
    int             p_ngens,
    VHG::MoTran&    p_moebius
    )
{
    g0_generators[p_ngens].a.x = p_moebius.e_a.real(); g0_generators[p_ngens].a.y = p_moebius.e_a.imag();
    g0_generators[p_ngens].b.x = p_moebius.e_b.real(); g0_generators[p_ngens].b.y = p_moebius.e_b.imag();
    g0_generators[p_ngens].c.x = p_moebius.e_c.real(); g0_generators[p_ngens].c.y = p_moebius.e_c.imag();
    g0_generators[p_ngens].d.x = p_moebius.e_d.real(); g0_generators[p_ngens].d.y = p_moebius.e_d.imag();
}











void read_stack()
{
    nstack = ngens = 0;

    VHG::VHG_circle circ; circ.c.x = 0.f; circ.c.y = 0.f; circ.r = 1.f;
    push_circle(circ);


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //   
    //    Note on the four (4) canonical Moebius Transformations: 
    //   
    //    Observation 1: 
    //              m1_C = complex_conjugate( m1_B ); 
    //              m1_D = complex_conjugate( m1_A ); 
    //    The property above also holds true for m2. 
    //              
    //      Consequently, Moebius Transformations motr1 and motr2 
    //      are fully specified BY ONLY TWO COMPLEX NUMBERS 
    //      to be named canon_A and canon_B. 
    //    
    //    
    //  
    //    Observation 2: 
    //    On a larger scale, the two Moebius Transformations m1 and m2 
    //    are related by 
    //                         m1 = g(A,B,C,D); 
    //                         m2 = g( cc(A),  cc(B),  cc(C),  cc(D) );
    //    In other words, Moebius Transformation m2 is built from 
    //    the four complex conjugates of the complex numbers used to build m1. 
    //    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //              
    //    Specification of Moebius transformations motr1 and motr2 
    //    requires only two complex numbers canon_A and canon_B: 
    //              
    std::complex<double> canon_A(+1.0, +1.0);
    std::complex<double> canon_B( 0.0, +1.0);
    //              
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    VHG::MoTran motr1( 
        canon_A, 
        canon_B, 
        std::conj(canon_B), 
        std::conj(canon_A)
    ); 
    std::complex<double> determinant_motr1 = motr1.MatrixDeterminant(); 
    gv_moebius_to_matrix(ngens, motr1); 
    ngens++;

    VHG::MoTran motr2( 
        std::conj(canon_A), 
        std::conj(canon_B), 
        canon_B,   //  conjugate(conjugate(canon_B)) = canon_B; 
        canon_A    //  conjugate(conjugate(canon_A)) = canon_A; 
    ); 
    std::complex<double> determinant_motr2 = motr2.MatrixDeterminant(); 
    gv_moebius_to_matrix(ngens, motr2); 
    ngens++;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //              
    //    Specification of Moebius transforms motr3 and motr4 
    //    requires only three complex numbers:
    //              
    std::complex<double> canon_P(+0.955, -0.025);
    std::complex<double> canon_Q(+0.045, +0.025); 
    std::complex<double> canon_R(-1.955, +0.025); 
    //              
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    VHG::MoTran motr3( 
        canon_P, 
        canon_Q, 
        canon_R, 
        canon_P 
    ); 
    std::complex<double> determinant_motr3 = motr3.MatrixDeterminant(); 
    std::complex<double> fixedpoint_motr3 = sqrt(
        canon_Q / canon_R
    );
    gv_moebius_to_matrix(ngens, motr3); 
    ngens++;


#ifdef GHV_OPTION_WEIRD_FRACTAL
    VHG::MoTran motr4( 
        canon_P,  
        -canon_R,   //  weird : swap positions of canon_R and canon_Q...
        -canon_Q,   //  weird : swap positions of canon_R and canon_Q...
        canon_P 
    ); 
#else
    VHG::MoTran motr4( 
        canon_P,  
        -canon_Q,   //  std::complex has the unary operator- 
        -canon_R,   //  std::complex has the unary operator- 
        canon_P 
    ); 


#endif


    std::complex<double> determinant_motr4 = motr4.MatrixDeterminant(); 
    gv_moebius_to_matrix(ngens, motr4); 
    ngens++;
}
//  Closes read_stack(); 






void generate_to_final()
{
    int nevo;


    //   Count of Circles = nstack, 
    //   Count of Generators = ngens;


    /* Put original circles on final list */

    for (int idx_a = 0; idx_a < nstack; idx_a++)
    {
        final[idx_a] = stack[idx_a];
    }
    nfinal = nstack;

    unique_circles(final, &nfinal);

    for (UINT idx_b = 1; idx_b <= VHG::EineKleine::k_depth; idx_b++)
    {
        if (nstack == 0)
        {
            break;  //  Generation is now complete.
        }
        nevo = nstack;


        evolve_stack();


        //   fprintf(stderr, "Gen %3d: final %6d evolving %6d stack %6d\n", i, nfinal, nevo, nstack);


        if (merge_circles())
        {
            break;  //  fprintf(stderr, "Stack full\n");
        }
    }

    //  if (verbose && i > VHG::EineKleine::k_depth) fprintf(stderr, "Depth attained\n");
}
//  Closes generate_to_final(); 




DirectX::XMFLOAT3 gv_FindDiskCenter(DirectX::XMFLOAT3 p_circle_normal, float p_circle_radius)
{
    float sphere_radius = Hvy3DScene::g_S2SphereRadius;

    float h_distance = sqrt(sphere_radius * sphere_radius - p_circle_radius * p_circle_radius); 

    float length_of_normal = VHG::torvec_vabs(p_circle_normal); 

    DirectX::XMFLOAT3 disk_center = DirectX::XMFLOAT3(
        h_distance * p_circle_normal.x / length_of_normal,
        h_distance * p_circle_normal.y / length_of_normal,
        h_distance * p_circle_normal.z / length_of_normal
    ); 
    return disk_center;
}
//  Closes gv_FindDiskCenter(); 



void gv_ComputeCircleFromNormal(
    DirectX::XMFLOAT3           p_circle_normal,
    float                       p_circle_radius,
    DirectX::XMFLOAT3&          p_R,
    DirectX::XMFLOAT3&          p_S
)
{
    //   Choose a point NOT on the line defined by p_circle_normal: 

    float dominant_f = max(
        abs(p_circle_normal.x), 
        max( abs(p_circle_normal.y), abs(p_circle_normal.z) )
    ); 
    XMFLOAT3 pt_any_fl3 = XMFLOAT3(dominant_f, 3.f * dominant_f, 5.f * dominant_f); 

    XMFLOAT3 pt_disk_center = gv_FindDiskCenter(p_circle_normal, p_circle_radius); 

    XMFLOAT3 r_any_fl3 = VHG::torvec_vsubtract(pt_any_fl3, pt_disk_center);
    XMVECTOR r_any_xmv = XMLoadFloat3(&r_any_fl3);

    XMVECTOR circle_normal_xmv = XMLoadFloat3(&p_circle_normal);

    XMVECTOR R_xmv = XMVector3Cross(r_any_xmv, circle_normal_xmv); 
    XMVECTOR S_xmv = XMVector3Cross(R_xmv, circle_normal_xmv); 

    XMVECTOR unit_R_xmv = XMVector3Normalize(R_xmv);
    XMVECTOR unit_S_xmv = XMVector3Normalize(S_xmv);

    XMStoreFloat3(&p_R, unit_R_xmv);
    XMStoreFloat3(&p_S, unit_S_xmv);
}
//  Closes gv_ComputeCircleFromNormal(); 














void gv_CircleSolid(
    std::vector<VHG_Instance>          *p_vect_Instances, 
    DirectX::XMFLOAT3                   p_circle_normal, 
    float                               p_circle_radius
)
{
    VHG_Instance            current_instance;
    XMFLOAT3                disk_point; 
    XMFLOAT3                basis_R; 
    XMFLOAT3                basis_S;

    DirectX::XMFLOAT3 disk_center = gv_FindDiskCenter(p_circle_normal, p_circle_radius); 

    gv_ComputeCircleFromNormal(p_circle_normal, p_circle_radius, basis_R, basis_S); 

    UINT const n_radial_steps = 256;  //  128;

    for (UINT idx_radius = 0; idx_radius < n_radial_steps; idx_radius++)
    {
        float radius_now = p_circle_radius - (idx_radius * p_circle_radius / (float)n_radial_steps); 

        float arc_density_now_fl = VHG::EineKleine::k_arc_density * radius_now / p_circle_radius;

        unsigned int arc_density_now_ui = lrint(arc_density_now_fl);



        for (UINT idx_azimuth = 0; idx_azimuth < arc_density_now_ui; idx_azimuth++)
        {
            float angle_azimuth = idx_azimuth * XM_2PI / arc_density_now_fl;

            disk_point.x = disk_center.x
                + radius_now * cosf(angle_azimuth) * basis_R.x
                + radius_now * sinf(angle_azimuth) * basis_S.x;

            disk_point.y = disk_center.y
                + radius_now * cosf(angle_azimuth) * basis_R.y
                + radius_now * sinf(angle_azimuth) * basis_S.y;

            disk_point.z = disk_center.z
                + radius_now * cosf(angle_azimuth) * basis_R.z
                + radius_now * sinf(angle_azimuth) * basis_S.z;

            current_instance.inst_pos = disk_point;
            current_instance.inst_attributes.x = radius_now;
            current_instance.inst_attributes.y = -2.f;
            p_vect_Instances->push_back(current_instance);
        }
    }
}
//  Closes gv_CircleSolid();  








void gv_CircleTo3DProper(
    std::vector<VHG_Instance>          *p_vect_Instances, 
    DirectX::XMFLOAT3                   p_circle_normal, 
    float                               p_circle_radius
)
{
    VHG_Instance            current_instance;
    XMFLOAT3                disk_point; 
    XMFLOAT3                basis_R; 
    XMFLOAT3                basis_S;

    DirectX::XMFLOAT3 disk_center = gv_FindDiskCenter(p_circle_normal, p_circle_radius); 

    gv_ComputeCircleFromNormal(p_circle_normal, p_circle_radius, basis_R, basis_S); 

    for (UINT idx_azimuth = 0; idx_azimuth < VHG::EineKleine::k_arc_density; idx_azimuth++)
    {
        float angle_azimuth = idx_azimuth * XM_2PI / (float)VHG::EineKleine::k_arc_density; 

        disk_point.x = disk_center.x
            + p_circle_radius * cosf(angle_azimuth) * basis_R.x
            + p_circle_radius * sinf(angle_azimuth) * basis_S.x; 

        disk_point.y = disk_center.y
            + p_circle_radius * cosf(angle_azimuth) * basis_R.y
            + p_circle_radius * sinf(angle_azimuth) * basis_S.y; 

        disk_point.z = disk_center.z
            + p_circle_radius * cosf(angle_azimuth) * basis_R.z
            + p_circle_radius * sinf(angle_azimuth) * basis_S.z; 

        current_instance.inst_pos = disk_point;
        current_instance.inst_attributes.x = p_circle_radius;
        current_instance.inst_attributes.y = +2.f;
        p_vect_Instances->push_back(current_instance);

        //  Rotate the entire frame to render the back face of S^2 sphere: 

        XMVECTOR disk_point_xmv = XMLoadFloat3(&disk_point); 
        XMMATRIX rotation_mx = XMMatrixRotationY(DirectX::XM_PI); 
        XMVECTOR antipode_xmv = XMVector3TransformCoord(disk_point_xmv, rotation_mx); 
        XMFLOAT3 antipode_fl3; 
        XMStoreFloat3(&antipode_fl3, antipode_xmv); 
        current_instance.inst_pos = antipode_fl3;
        current_instance.inst_attributes.x = p_circle_radius;
        current_instance.inst_attributes.y = +2.f;
        p_vect_Instances->push_back(current_instance);
    }
}
//  Closes gv_CircleTo3DProper()










void gv_CircleTo3DEconomy(
    std::vector<VHG_Instance>          *p_vect_Instances, 
    DirectX::XMFLOAT3                   p_circle_normal, 
    float                               p_circle_radius
)
{
    VHG_Instance            current_instance;
    XMFLOAT3 disk_point; 
    DirectX::XMFLOAT3 disk_center = gv_FindDiskCenter(p_circle_normal, p_circle_radius); 


#ifdef GHV_OPTION_SPHERICAL_RENDER
    disk_point.x = disk_center.x; 
    disk_point.y = disk_center.y; 
    disk_point.z = disk_center.z; 

#else

    float s2u = disk_center.x;
    float s2v = disk_center.y;
    float s2w =
        (abs(1.00f - disk_center.z) < 0.000001f) ? disk_center.z + 0.0001f : disk_center.z;

    disk_point.x = s2u / (1.f - s2w); 

    disk_point.y = s2v / (1.f - s2w); 

    disk_point.z = 3.f;

#endif


    float obfuscated_limit = 0.05f;  

    // obfuscated_limit = 0.9f;  //  undo : experimental 

    if (p_circle_radius < obfuscated_limit)  //  Most circles have p_circle_radius < 0.1f; 
    {
        current_instance.inst_pos = disk_point;
        current_instance.inst_attributes.x = p_circle_radius;
        current_instance.inst_attributes.y = +2.f;
        p_vect_Instances->push_back(current_instance);

        //  Rotate the entire frame to render the back face of S^2 sphere: 

        XMVECTOR disk_point_xmv = XMLoadFloat3(&disk_point); 
        XMMATRIX rotation_mx = XMMatrixRotationY(DirectX::XM_PI); 
        XMVECTOR antipode_xmv = XMVector3TransformCoord(disk_point_xmv, rotation_mx); 
        XMFLOAT3 antipode_fl3; 
        XMStoreFloat3(&antipode_fl3, antipode_xmv); 
        current_instance.inst_pos = antipode_fl3;
        current_instance.inst_attributes.x = p_circle_radius;
        current_instance.inst_attributes.y = +2.f;
        p_vect_Instances->push_back(current_instance);
    }
}
//  Closes gv_CircleTo3DEconomy()











void ps_plot_circles_s2(
    std::vector<VHG_Instance>    *p_vect_Instances, 
    VHG::VHG_circle p_the_stack[], 
    int nstack
)
{
    //   Kleinian group limit set on sphere  

    for (int i = 0; i < nstack; i++)
    {
        DirectX::XMFLOAT3   v_3DCircleNormal;
        double              r_3DCircleRadius;

        v_3DCircleNormal = circle_to_s2(p_the_stack[i], &r_3DCircleRadius);


        if (r_3DCircleRadius > circle_rad_hwm)
        {
            circle_rad_hwm = r_3DCircleRadius;
        }



        //   Choose one of the following two rendering methods: 


#ifdef GHV_OPTION_MESH_POINT

        gv_CircleTo3DEconomy(p_vect_Instances, v_3DCircleNormal, (float)r_3DCircleRadius);

#else

        gv_CircleTo3DProper(p_vect_Instances, v_3DCircleNormal, (float)r_3DCircleRadius);
#endif

    }
}




uint32_t Hvy3DScene::gv_MeshMcMullenCusp( std::vector<VHG_Instance>    *p_vect_Instances)
{
    //      ghv : set variables to emulate "cusp.run" with mode = SPHERE;

    circle_bytes = sizeof(VHG::VHG_circle);

    UINT gg = 7;


    read_stack();  //  initializes four (4) Moebius transformations;
    generate_to_final();
    ps_plot_circles_s2(p_vect_Instances, final, nfinal);


    // GOLD : gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.f, 0.f, 0.f), 0.5f);  // target aligned with x-axis;

    // EXCELLENT : gv_CircleSolid(p_vect_Instances, XMFLOAT3(0.f, 0.f, 1.f), 0.5f); 

    //  gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.f, 0.f, 1.f), 0.5f); 
    //  gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.7f, 0.f, 1.f), 0.5f); 
    //  gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.7f, 0.f, 1.f), 0.4f); 
    //  gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.866f, 0.f, 1.f), 0.43f); 
    //  gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.877f, 0.f, 1.f), 0.45f); 
    //    Let's try an intentionally smaller circle: 

    //  Normal is VERY close : gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.877f, 0.f, 1.f), 0.43f); 

    gv_CircleSolid(p_vect_Instances, XMFLOAT3(1.881f, 0.f, 1.f), 0.43f); 

    UINT card_instance_cubes = (uint32_t)p_vect_Instances->size();

    return card_instance_cubes;
}
//  Closes VHG_Scene3D::gv_MeshMcMullenCusp(); 





