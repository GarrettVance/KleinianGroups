#pragma once

#include "..\Common\DeviceResources.h"
#include "..\Common\StepTimer.h"

#include "GeometricPrimitive.h"
#include "WICTextureLoader.h" 
#include "Keyboard.h"   //  ghv : for DirectX::Keyboard;
#include "Mouse.h"      //  ghv : for DirectX::Mouse;


namespace KG
{
    struct ModelViewProjectionConstantBuffer
    {
        DirectX::XMFLOAT4X4     model;
        DirectX::XMFLOAT4X4     view;
        DirectX::XMFLOAT4X4     projection;
        DirectX::XMUINT4        animator_count;
    };

    struct VHG_Vertex_PosTex
    {
        DirectX::XMFLOAT3       e_pos;
        DirectX::XMFLOAT2       e_texco;
    };

    struct VHG_Instance
    {
        DirectX::XMFLOAT3       inst_pos;
        DirectX::XMFLOAT3       inst_attributes;
    };

    class VHG_Scale
    {
    public:
        VHG_Scale(float p_scale_factor) : e_scale_factor(p_scale_factor) {}


        void posApply(std::vector<VHG_Vertex_PosTex> & p_vectpostex)
        {
            for_each(
                p_vectpostex.begin(),
                p_vectpostex.end(),

                //      Important : pass p_postex BY REFERENCE 
                //      otherwise won't be able to alter it!!! 

                [this](VHG_Vertex_PosTex & p_postex) {
                    p_postex.e_pos.x *= this->e_scale_factor;
                    p_postex.e_pos.y *= this->e_scale_factor;
                    p_postex.e_pos.z *= this->e_scale_factor;
                }
            );
        }

    private:
        float e_scale_factor;
    };




	class Hvy3DScene
	{
	public:
		Hvy3DScene(const std::shared_ptr<DX::DeviceResources>& deviceResources);

        static float g_S2SphereRadius; 

		void CreateDeviceDependentResources();
		void CreateWindowSizeDependentResources();
		void ReleaseDeviceDependentResources();


		void Update(DX::StepTimer const& timer);
		void Render();


        UINT gv_get_instance_idx(void) { return m_constantBufferData.animator_count.x;  }

	private:

        uint32_t gv_MeshHyperboloid(std::vector<VHG_Instance>    *p_vect_Instances);

        uint32_t gv_MeshMcMullenCusp(std::vector<VHG_Instance>    *p_vect_Instances); 

        void gv_mesh_for_point();

        void gv_mesh_for_cube();

        void gv_rotate_cloud(double p_total_seconds);

	private:
		// Cached pointer to device resources.
		std::shared_ptr<DX::DeviceResources> m_deviceResources;


		Microsoft::WRL::ComPtr<ID3D11InputLayout>	m_inputLayout;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_vertexBuffer;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_indexBuffer;
		Microsoft::WRL::ComPtr<ID3D11VertexShader>	m_vertexShader;
		Microsoft::WRL::ComPtr<ID3D11PixelShader>	m_pixelShader;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_constantBuffer;
		ModelViewProjectionConstantBuffer	                    m_constantBufferData;
        std::unique_ptr<DirectX::Keyboard>                      m_keyboard;
        std::unique_ptr<DirectX::Mouse>                         m_mouse;
        bool                                                    e_option_rotate_world;
        bool                                                    e_option_wireframe;
		bool	                                                m_loadingComplete;
		float	                                                m_degreesPerSecond;
        Microsoft::WRL::ComPtr<ID3D11SamplerState>              e_texture_sampler_state;
        Microsoft::WRL::ComPtr<ID3D11ShaderResourceView>        e_texture_srv_1;


        Microsoft::WRL::ComPtr<ID3D11Buffer>                    m_instanceBuffer; 
        uint32_t                                                m_instanceCount;
        uint32_t                                                m_vertexCount;
        uint32	                                                m_indexCount;
        std::unique_ptr<DirectX::GeometricPrimitive>            e_primitive;
	};
}

