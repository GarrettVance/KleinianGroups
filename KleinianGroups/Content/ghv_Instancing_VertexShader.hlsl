//         
// 
//                 ghv Instancing Vertex Shader 
// 
            
cbuffer WVP_ConstantBufferStruct : register(b0)
{
	matrix   model;
	matrix   view;
	matrix   projection;
    uint4    animator_count;
};



struct VertexShaderInput
{
	float3    pos                   : POSITION;
	float2    texco                 : TEXCOORD0; 
    float3    inst_pos              : TEXCOORD1; 
    float3    inst_attributes       : TEXCOORD2; 
    uint      inst_id               : SV_InstanceID; 
    uint       vert_id              : SV_VertexID;
};



struct PixelShaderInput
{
	float4     pos         : SV_POSITION;
	float2     texco       : TEXCOORD0;
    float4     anim_color  : COLOR; 
    float      circle_radius : TANGENT;
};



PixelShaderInput main(VertexShaderInput input)
{
     PixelShaderInput output;
     
     float4 vertex_location = float4(input.pos, 1.0f);

     float4 tmp_anim_color = float4(0.f, 0.f, 0.f, 0.f); 
     

    //   Update the position of the vertices 
    //   based on the data for this particular instance.

    vertex_location.x += input.inst_pos.x;
    vertex_location.y += input.inst_pos.y;
    vertex_location.z += input.inst_pos.z;



    //     
    //     obtain spherical coordinates 
    //     from the Cartesian coords in "vertex_location"
    //     

    float r_radial = sqrt(
        vertex_location.x * vertex_location.x +
        vertex_location.y * vertex_location.y +
        vertex_location.z * vertex_location.z);

    float phi_inclination = acos(vertex_location.y / r_radial); 

    float theta_azimuth = atan2(vertex_location.z, vertex_location.x);

    float k_pi = 3.1415926535; 


    //  GOLD : float phi1 = 45.f * k_pi / 180.f; 
    //  GOLD : float phi2 = 135.f * k_pi / 180.f; 

    //  very close : float phi1 = 60.f * k_pi / 180.f; 
    //  very close : float phi2 = 120.f * k_pi / 180.f; 

    // too small : float phi1 = 70.f * k_pi / 180.f; 
    // too small : float phi2 = 110.f * k_pi / 180.f; 

    float phi1 = 62.f * k_pi / 180.f; 
    float phi2 = 118.f * k_pi / 180.f; 



    // float az0 = 0.f; 
    // float az0 = 30.f * k_pi / 180.f;
    //  float az0 = -30.f * k_pi / 180.f;
    // float az0 = -45.f * k_pi / 180.f;
    // float az0 = -59.f * k_pi / 180.f;
    //  float az0 = -60.f * k_pi / 180.f;
    float az0 = -62.f * k_pi / 180.f;


    // GOLD : float az1 = az0 + 45.f * k_pi / 180.f; 
    // GOLD : float az2 = az0 + 135.f * k_pi / 180.f; 

    //  very close : float az1 = az0 + 60.f * k_pi / 180.f; 
    //  very close : float az2 = az0 + 120.f * k_pi / 180.f; 

    //  too small : float az1 = az0 + 70.f * k_pi / 180.f; 
    //  too small : float az2 = az0 + 110.f * k_pi / 180.f; 


    float az1 = az0 + 62.f * k_pi / 180.f; 
    float az2 = az0 + 118.f * k_pi / 180.f; 



    float tex_u = 0.f;
    float tex_v = 0.f;

    if ( (theta_azimuth > az1) && (theta_azimuth < az2) && 
         (phi_inclination > phi1) && (phi_inclination < phi2) )
    {
        //   tex_u = 2.f * theta_azimuth / 3.1415926535;

        tex_u = (theta_azimuth - az1) / (az2 - az1); 

        tex_v = (phi_inclination - phi1) / (phi2 - phi1); 
    }




    if(input.inst_id == animator_count.x)
    {
         // vertex_location.x *= 2.f;

         tmp_anim_color.x = 2.f;
    }


    tmp_anim_color.y = (float)input.inst_id;
    tmp_anim_color.z = (float)input.vert_id;

     // Transform the vertex position "vertex_location" into projected space.
     
     vertex_location = mul(vertex_location, model);
     vertex_location = mul(vertex_location, view);
     vertex_location = mul(vertex_location, projection);
     output.pos = vertex_location;
     
     //   not  output.texco = input.texco;

     output.texco = float2(tex_u, tex_v);


     output.anim_color = tmp_anim_color;


     //   inst_attributes.y encodes whether or not 
     //   the vertex belongs to the solid disk used 
     //   only to paint cameo of Doctor Moebius. 
     //  
     //   inst_attributes.y < 0.f ===>  vertex belongs to solid disk;
     // 
     //   inst_attributes.y > 0.f ===>  its just a regular Julia Set vertex;


     //  output.circle_radius = input.inst_attributes.x;

     output.circle_radius = input.inst_attributes.y;



     return output;
}




