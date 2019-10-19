//               
//              
//                     ghv Instancing Pixel Shader  
//     


Texture2D ObjTexture;
SamplerState ObjSamplerState;



struct PixelShaderInput
{
	float4     pos         : SV_POSITION;
	float2     texco       : TEXCOORD0;
    float4     anim_color  : COLOR; 
    uint       prim_id     : SV_PrimitiveID; 
    float      circle_radius : TANGENT;
};




float4       main(PixelShaderInput      input) : SV_TARGET
{
     float4 v_rgba;  

     
     if(input.anim_color.x > 0.5)
     {
         v_rgba.r = 0.9f;
         v_rgba.g = 0.9f;
         v_rgba.b = 0.9f;
     }

     //  
     //    Using SV_InstanceID: 
     //   
     float arbiter_f = input.anim_color.y; 
     float magnitude_f = 1e5;   //  for anim_color.y;

     magnitude_f = 1.2e5;

     //  
     //    Using SV_PrimitiveID: 
     //   
     //  float arbiter_f = (float)input.prim_id;
     //  float magnitude_f = 3;   //  for primitive_id;



     //  
     //    Using SV_VertexID: 
     //   
     //  float arbiter_f = (float)input.anim_color.z;
     //  float magnitude_f = 3;
     //  magnitude_f = 4; 


     if (arbiter_f > 2.04f * magnitude_f)  //  Using 2.04 is 99 percent perfect; 
     {
         v_rgba = float4(1.f, 0.f, 0.f, 1.f); // red; 
     }
     else if (arbiter_f > 1.65 * magnitude_f)  //   1.65 is better than 1.5;

     {
         v_rgba = float4(0.f, 1.f, 0.f, 1.f);  // green;
     }
     else if (arbiter_f > 1.2f * magnitude_f)  //  Using 1.2f is essentially perfect to 99.9 percent;
     {
         v_rgba = float4(0.f, 1.f, 0.f, 1.f);  // green; 
     }
     else 
     {
         v_rgba = float4(0.f, 0.f, 1.f, 1.f);  // blue; 
     }



     if (input.circle_radius < 0.f)
     {
         //   The solid disk used to paint the cameo portrait: 

         v_rgba = ObjTexture.Sample(ObjSamplerState, input.texco);
     }


     return v_rgba;
 
}




