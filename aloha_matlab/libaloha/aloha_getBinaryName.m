function binary_name = aloha_getBinaryName(arch, version)
%  Determine which is the architecture of the build platform
%  
%  binary_name = aloha_getBinaryName(arch, version)
%  
%  INPUTS
%   - arch [str] : architecture of the platform who is running the code
%   - version [integer] : version of the binary code 
%
%  OUTPUTS
%   - binary_name [str] : binary name to be launched
%  
%  AUTHOR(S) : JH
%  LAST UPDATE : 
%   - 17/06/2008 [creation]
%  

    binary_name = ['coupl_plasma_version', num2str(version), '_', arch];