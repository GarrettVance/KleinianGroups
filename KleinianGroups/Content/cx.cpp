

#include "pch.h"

#include <math.h>
#include "cx.h"

using namespace VHG; 




/* We consider SL_2(C) adjoined orientation-reversing
maps.  The latter are distinguished by det=-1
and represent f(z) = (aw+b)/(cw+d) where w=z-bar.
We compute the action on circles.

Lines a represented by circles (c,r) with r < 0 and
with c = closest point to the origin, and with
r = angle in revolutions [-1,0].
*/












VHG_complex VHG::add (VHG_complex z, VHG_complex w)
{
	VHG_complex t;
	t.x = z.x + w.x;
	t.y = z.y + w.y;
	return(t);
}





VHG_complex VHG::sub (VHG_complex z, VHG_complex w)
{
	VHG_complex t;
	t.x = z.x - w.x;
	t.y = z.y - w.y;
	return(t);
}








VHG_complex VHG::mult (VHG_complex z, VHG_complex w)
{
	VHG_complex t;
	t.x = z.x*w.x - z.y*w.y;
	t.y = z.x*w.y + z.y*w.x;
	return(t);
}











VHG_complex VHG::divide (VHG_complex z, VHG_complex w)
{
	return(
        mult(z, 
            VHG::recip(w)
        )
    );
}





VHG_complex VHG::recip (VHG_complex z)
{
	VHG_complex w;
	double r;
	r = z.x*z.x + z.y*z.y;
	w.x = z.x / r;
	w.y = -z.y / r;
	return(w);
}






VHG_complex VHG::cx_conj(VHG_complex z)
{
	VHG_complex t;
	t = z;
	t.y = -t.y;
	return(t);
}







VHG_complex VHG::cx_sqrt(VHG_complex z)
{
	VHG_complex w;
	
    // double fabs(), sqrt();



    /* Worry about numerical stability */

	if (z.x == 0.0 && z.y == 0.0) return(z); 
	else
	if (z.x > fabs(z.y))
		{
		w.x = sqrt((z.x+sqrt(z.x*z.x+z.y*z.y))/2);
		w.y = z.y/(2*w.x);
		}
	else
		{
		w.y = sqrt((-z.x+sqrt(z.x*z.x+z.y*z.y))/2);
		w.x = z.y/(2*w.y);
		}
	return(w);
}











/* Compute sqrt(z) in the half-plane perpendicular to w. */
VHG_complex VHG::contsqrt(VHG_complex z, VHG_complex w)
{
    //  VHG_complex cx_sqrt();


    VHG_complex t;

	t = VHG::cx_sqrt(z);


	if (0 > (t.x*w.x + t.y*w.y)) 
    {
        t.x = -t.x; 
        t.y = -t.y;
    }

    return(t);
}











double VHG::cx_abs (VHG_complex z)       /* L 2 norm of z */
{
	return (sqrt (z.x*z.x + z.y*z.y));
}











double VHG::infnorm (VHG_complex z)   /* L infinity norm of z */
{
	double a,b;
	a = (z.x > 0) ? z.x : -z.x;
	b = (z.y > 0) ? z.y : -z.y;
	return ((a>b) ? a : b);
}










VHG_complex VHG::polar (double radius, double angle) /*Convert to VHG_complex. */
{
	VHG_complex z;
	z.x = cos (angle) * radius;
	z.y = sin (angle) * radius;
	return(z);
}












/* Values in [-pi,pi]. */
double VHG::arg(VHG_complex z)  
{
	return(atan2(z.y,z.x));  // namely the std atan2...
}










VHG_complex VHG::cx_exp(VHG_complex z)
{
	VHG_complex w;
	double m;

	m = exp(z.x);
	w.x = m * cos(z.y);
	w.y = m * sin(z.y);
	return(w);
}










VHG_complex VHG::cx_log(VHG_complex z)
{
	VHG_complex w;

	w.x = log(cx_abs(z));
	w.y = arg(z);
	return(w);
}









VHG_complex VHG::cx_sin(VHG_complex z)
{
	VHG_complex w;

	w.x = sin(z.x) * cosh(z.y);
	w.y = cos(z.x) * sinh(z.y);
	return(w);
}










VHG_complex VHG::cx_cos(VHG_complex z)
{
	VHG_complex w;

	w.x = cos(z.x) * cosh(z.y);
	w.y = -sin(z.x) * sinh(z.y);
	return(w);
}








VHG_complex VHG::cx_sinh(VHG_complex z)
{
	VHG_complex w;

	w.x = sinh(z.x) * cos(z.y);
	w.y = cosh(z.x) * sin(z.y);
	return(w);
}





VHG_complex VHG::cx_cosh(VHG_complex z)
{
	VHG_complex w;

	w.x = cosh(z.x) * cos(z.y);
	w.y = sinh(z.x) * sin(z.y);
	return(w);
}












VHG_complex VHG::power(VHG_complex z, double t)   /* Raise z to a real power t */
{	
    // double arg(), cx_abs(), pow();

	// VHG_complex polar();

	return(
        VHG::polar(
            pow(
                VHG::cx_abs(z), 
                t
            ), 
            t * arg(z)
        )
        );
}














/*  Map points in the unit disk onto the lower hemisphere of the
    Riemann sphere by inverse stereographic projection.
    Projecting, r -> s = 2r/(r^2 + 1);  inverting this,
    s -> r = (1 - sqrt(1-s^2))/s.   */
 
VHG_complex VHG::disk_to_sphere(VHG_complex z)
{       
        VHG_complex w;
        double r, s; 


        // double cx_abs(), 
        // double sqrt();

 
        s = VHG::cx_abs(z);
        if (s == 0) return(z);
                else r = (1 - sqrt(1-s*s))/s;
        w.x = (r/s)*z.x;
        w.y = (r/s)*z.y;
        return(w);
}









VHG_complex VHG::mobius(
    VHG_complex a, 
    VHG_complex b, 
    VHG_complex c, 
    VHG_complex d, 
    VHG_complex z
)
{       
    return(divide(
            add(mult(a,z),b),
            add(mult(c,z),d)
            )
    );
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





VHG_complex VHG::det(VHG_matrix x)
{
    return(sub(mult(x.a, x.d), mult(x.b, x.c)));
}








VHG_matrix VHG::make_sl2(VHG_matrix x)
{
    VHG_complex d;

    d = recip(cx_sqrt(det(x)));

    x.a = mult(x.a, d);
    x.b = mult(x.b, d);
    x.c = mult(x.c, d);
    x.d = mult(x.d, d);

    return(x);
}










//    begin f.cpp



VHG_matrix VHG::conjmat(VHG_matrix x)
{
	x.a = cx_conj(x.a);
	x.b = cx_conj(x.b);
	x.c = cx_conj(x.c);
	x.d = cx_conj(x.d);
	return(x);
}




VHG_matrix VHG::inverse(VHG_matrix x)
{
	VHG_matrix z;
	VHG_complex d;
	d = det(x);
	z.a = x.d;
	z.b.x = -x.b.x;
	z.b.y = -x.b.y;
	z.c.x = -x.c.x;
	z.c.y = -x.c.y;
	z.d = x.a;
	if(d.x > 0) return(z);
	z.a.x = -z.a.x;
	z.b.x = -z.b.x;
	z.c.x = -z.c.x;
	z.d.x = -z.d.x;
	return(z);
}






VHG_matrix VHG::product(VHG_matrix x, VHG_matrix y)
{
	VHG_matrix z;
	VHG_complex d;

	d = det(x);
	if(d.x < 0) y = conjmat(y);
	z.a = add(mult(x.a,y.a),mult(x.b,y.c));
	z.b = add(mult(x.a,y.b),mult(x.b,y.d));
	z.c = add(mult(x.c,y.a),mult(x.d,y.c));
	z.d = add(mult(x.c,y.b),mult(x.d,y.d));
	return(z);
}







VHG_matrix VHG::circle_to_matrix(VHG_circle circ)
{
	VHG_matrix z;
	if(circ.r <= 0) return(line_to_matrix(circ));
	z.a.x = circ.c.x/circ.r; 
	z.a.y = circ.c.y/circ.r;
	z.b.x = circ.r-(circ.c.x*circ.c.x+circ.c.y*circ.c.y)/circ.r;
	z.b.y = 0;
	z.c.x = 1/circ.r;
	z.c.y = 0;
	z.d.x = -circ.c.x/circ.r;
	z.d.y =  circ.c.y/circ.r;
	return(z);
}







VHG_matrix VHG::line_to_matrix(VHG_circle circ)
{
	VHG_matrix z;
	VHG_complex unit;
	static VHG_complex i = {0.0,1.0};

	unit.x = cos(-PI*circ.r);
	unit.y = sin(-PI*circ.r);
	z.a = mult(i,unit);
	z.b = mult(i,sub(  mult(cx_conj(unit),circ.c),
			   mult(unit,cx_conj(circ.c))   ));
	z.c.x = 0.0; z.c.y = 0.0; 
	z.d = mult(i,cx_conj(unit));
	return(z);
}











VHG_circle VHG::matrix_to_circle(VHG_matrix x)
{
	VHG_circle circ;

	if(x.c.x < LINEFUZZ && x.c.x > -LINEFUZZ) return(matrix_to_line(x));
	circ.c = divide(x.a,x.c);
	circ.r = 1/x.c.x;
	if(circ.r < 0) circ.r = -circ.r;
	return(circ);
}










VHG_circle VHG::matrix_to_line(VHG_matrix x)
{
	VHG_circle circ;

	circ.r =   -arg(x.a)/PI - 1.5;
	circ.c.x = -x.b.x/2; 
	circ.c.y = -x.b.y/2; 
	circ.c = mult(circ.c,x.a);
	return(circ);
}











VHG_circle VHG::image_circle(VHG_matrix x, VHG_circle p_input_circle)
{


    VHG_matrix  y; y = VHG::circle_to_matrix(p_input_circle);


    VHG_matrix  xi; xi = inverse(x); 


	y = product(x,product(y,xi));



	// undo : p_input_circle = VHG::matrix_to_circle(y);
	// undo : return(p_input_circle);

	VHG_circle other_circle = VHG::matrix_to_circle(y);
    return other_circle;

}








/* Sphere subroutines */


DirectX::XMFLOAT3 complex_to_s2(VHG_complex z)
{

    double r;

    r = z.x * z.x + z.y * z.y;

    //   ghv : obfuscate : TODO: looks like stereographic projection... 

    double x_dbl = 2 * z.x / (1 + r);
    double y_dbl = 2 * z.y / (1 + r);
    double z_dbl = (1 - r) / (1 + r);

    DirectX::XMFLOAT3 v = DirectX::XMFLOAT3(
        (float)x_dbl, 
        (float)y_dbl,
        (float)z_dbl
        );

    return(v);
}










DirectX::XMFLOAT3 VHG::torvec_vadd(DirectX::XMFLOAT3 p_a, DirectX::XMFLOAT3 p_b)
{
    DirectX::XMFLOAT3 retval = DirectX::XMFLOAT3(p_a.x + p_b.x, p_a.y + p_b.y, p_a.z + p_b.z); 
    return retval;
}


DirectX::XMFLOAT3 VHG::torvec_vsubtract(DirectX::XMFLOAT3 p_a, DirectX::XMFLOAT3 p_b)
{
    DirectX::XMFLOAT3 retval = DirectX::XMFLOAT3(p_a.x - p_b.x, p_a.y - p_b.y, p_a.z - p_b.z); 
    return retval;
}


float VHG::torvec_vabs(DirectX::XMFLOAT3 p_a)
{
    float f_norm_squared = p_a.x * p_a.x
        + p_a.y * p_a.y
        + p_a.z * p_a.z; 

    return sqrt(f_norm_squared); 
}







/* Convert line in plane to circle on S^2 */

DirectX::XMFLOAT3 line_to_s2(VHG_circle circ, double * r)
{
    DirectX::XMFLOAT3 p, v;


    static DirectX::XMFLOAT3 north_pole = { 0.0,  0.0,   -1.0 };


    // double norm;

    p = complex_to_s2(circ.c);



    if (p.x == 0.0 && p.y == 0.0)
    {
        //   Downgrade from double-precision to float:  

        v.x = (float)cos(-PI * circ.r);
        v.y = (float)sin(-PI * circ.r);
        v.z = (float)0.0;
        *r = 1.0;
    }
    else
    {
        v = torvec_vadd(p, north_pole);

        *r = torvec_vabs(
            torvec_vsubtract(p, north_pole)
        ) / 2;

    }
    return(v);
}

















/* Convert circle in plane to circle on S^2 */
/* Method is to find a pair of opposite points on circle */
/* Return norm vector; put radius in R^3 in r*/

DirectX::XMFLOAT3 VHG::circle_to_s2(VHG_circle circ, double * p_little_r)
{
    //          ghv : TODO:   here's where the switch from double to float begins...


    double ac, s;

    VHG_complex w1, w2, c, cn;


    //  formerly VHG_vector v, v1, v2;

    DirectX::XMFLOAT3    v, v1, v2;



    if (circ.r <= 0) return(line_to_s2(circ, p_little_r));

    c = circ.c;
    s = circ.r;
    ac = cx_abs(c);

    /* Circle centered at origin */

    if (ac == 0.0)
    {
        v.x = v.y = 0; v.z = 1;
        if (s>1.0) v.z = -1;
        *p_little_r = 2 * s / (1 + s * s);
        return(v);
    }


    /* Else find points closest and farthest from origin */

    cn.x = c.x / ac; cn.y = c.y / ac;
    w1.x = c.x + s * cn.x; w1.y = c.y + s * cn.y;
    w2.x = c.x - s * cn.x; w2.y = c.y - s * cn.y;

    v1 = complex_to_s2(w1);


    v2 = complex_to_s2(w2);


    v = torvec_vadd(v1, v2);



    *p_little_r = torvec_vabs(torvec_vsubtract(v1, v2)) / 2;



    /* Take care if we have a great circle */

    if (torvec_vabs(v) < VEPS) 
    {
        s = sqrt(v1.x*v1.x + v1.y*v1.y);

        //   Downgrade from double-precision to float: 

        v.x = -v1.x * v1.z / (float)s;
        v.y = -v1.y * v1.z / (float)s;
        v.z = (float)s;
    }
    return(v);
}











//   end of f.cpp





