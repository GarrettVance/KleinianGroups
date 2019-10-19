# KleinianGroups
A port of Curtis T. McMullen's "Lim (Version 2.2)" (Limit Sets of Kleinian Groups) application to C++ / DirectX D3D11.

This is a complete solution for Visual Studio 2017 which builds a DirectX D3D11 application for Windows 10 UWP. You can build for either x64 debug or x64 release.

Dependencies: The KleinianGroups project requires an external library, namely Microsoft DirectXTK. The easiest way to satisfy this dependency is to use the nuget console built into VisualStudio 2017. DirectXTK can be installed for the KleinianGroups solution by running the following command in the VisualStudio Package Manager Console: 

     PM> Install-Package directxtk_uwp -Version 2019.10.17.1

Alternately, if nuget is not available: DirectXTK is available on github as a VisualStudio solution targeting Windows 10 UWP. Build DirecXTK for x64, after which you can build the KleinianGroups project.
