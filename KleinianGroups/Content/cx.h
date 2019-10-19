
#pragma once 


namespace VHG 
{

#define PI 3.1415926535897932384  

#define TWO30  1073741824


#define VEPS 0.0000000001     /* Recognizing great circles */


#define INF 1e20

#define PARAFUZZ  1e-3  /* Nearly parallel slopes */

#define LINEFUZZ  1e-5  /* Line vs. Circle */

#define AGREE 100	/* 1/AGREE to match circles */







typedef struct VHG_complex
{
    double x; 
    double y;
} VHG_complex; 


typedef struct VHG_circle
{
    VHG_complex c;
    double r;
} VHG_circle;



typedef struct VHG_matrix
{
    VHG_complex a;
    VHG_complex b;
    VHG_complex c;
    VHG_complex d;
} VHG_matrix;





VHG::VHG_complex add (VHG_complex z, VHG_complex w);

VHG_complex sub (VHG_complex z, VHG_complex w); 
VHG_complex mult (VHG_complex z, VHG_complex w); 


VHG_complex divide (VHG_complex z, VHG_complex w); 
VHG_complex recip (VHG_complex z); 
VHG_complex cx_conj(VHG_complex z); 






VHG::VHG_complex cx_sqrt(VHG_complex z);

/* Compute sqrt(z) in the half-plane perpendicular to w. */

VHG_complex contsqrt(VHG_complex z, VHG_complex w); 






double cx_abs (VHG_complex z);       /* L 2 norm of z */
double infnorm (VHG_complex z);   /* L infinity norm of z */



VHG_complex polar (double radius, double angle);  /*Convert to VHG_complex. */



/* Values in [-pi,pi]. */

double arg(VHG_complex z);



VHG_complex cx_exp(VHG_complex z);
VHG_complex cx_log(VHG_complex z);
VHG_complex cx_sin(VHG_complex z); 
VHG_complex cx_cos(VHG_complex z); 
VHG_complex cx_sinh(VHG_complex z);
VHG_complex cx_cosh(VHG_complex z);


VHG_complex power(VHG_complex z, double t);   /* Raise z to a real power t */



/*  Map points in the unit disk onto the lower hemisphere of the
    Riemann sphere by inverse stereographic projection.
    Projecting, r -> s = 2r/(r^2 + 1);  inverting this,
    s -> r = (1 - sqrt(1-s^2))/s.   */
 
VHG_complex disk_to_sphere(VHG_complex z);



VHG_complex mobius(
    VHG_complex a, 
    VHG_complex b, 
    VHG_complex c, 
    VHG_complex d, 
    VHG_complex z
);



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



VHG_complex det(VHG_matrix x);
VHG_matrix make_sl2(VHG_matrix x);
VHG_matrix conjmat(VHG_matrix x);
VHG_matrix inverse(VHG_matrix x);
VHG_matrix product(VHG_matrix x, VHG_matrix y);
VHG_matrix circle_to_matrix(VHG_circle circ);
VHG_matrix line_to_matrix(VHG_circle circ);
VHG_circle matrix_to_circle(VHG_matrix x);
VHG_circle matrix_to_line(VHG_matrix x);
VHG_circle image_circle(VHG_matrix x, VHG_circle c);









/* Convert circle in plane to circle on S^2 */
/* Method is to find a pair of opposite points on circle */
/* Return norm vector; put radius in R^3 in r*/

//   formerly : VHG_vector circle_to_s2(VHG_circle circ, double * r);

DirectX::XMFLOAT3 circle_to_s2(VHG_circle circ, double * r);




DirectX::XMFLOAT3   torvec_vadd(DirectX::XMFLOAT3 p_a, DirectX::XMFLOAT3 p_b); 
DirectX::XMFLOAT3   torvec_vsubtract(DirectX::XMFLOAT3 p_a, DirectX::XMFLOAT3 p_b); 
float               torvec_vabs(DirectX::XMFLOAT3 p_a);




}
//  Closes namespace VHG; 


