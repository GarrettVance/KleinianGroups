#pragma once

#include "Common\StepTimer.h"
#include "Common\DeviceResources.h"
#include "Content\gv_Scene3D.h"
#include "D2D1\gv_HUD.h"



namespace KleinianGroups
{
	class KleinianGroupsMain : public DX::IDeviceNotify
	{
	public:
		KleinianGroupsMain(const std::shared_ptr<DX::DeviceResources>& deviceResources);
		~KleinianGroupsMain();
		void CreateWindowSizeDependentResources();
		void Update();
		bool Render();

		// IDeviceNotify
		virtual void OnDeviceLost();
		virtual void OnDeviceRestored();

	private:
		// Cached pointer to device resources.
		std::shared_ptr<DX::DeviceResources> m_deviceResources;


		std::unique_ptr<KG::Hvy3DScene> m_sceneRenderer;

		std::unique_ptr<KG::HUD> m_fpsTextRenderer;


		// Rendering loop timer.
		DX::StepTimer m_timer;
	};
}
