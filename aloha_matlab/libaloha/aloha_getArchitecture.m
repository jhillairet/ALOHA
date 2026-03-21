function arch = aloha_getArchitecture
%  Determine which is the architecture of the build platform
%  
%  arch = aloha_getArchitecture
%  
%  INPUT : none
%  
%  OUTPUTS
%  - str : string which may be :
%      'glnx86' : GNU on x86
%      'glnxa64': GNU on x86_64
%      'alpha'  : DEC Alpha 
%      'pcwin'  : Windows 32 et 64 bits 
%  
%  COMMENTS : this programm could be use to choose which binary 
%  version should be used.
%  
%  AUTHOR(S) : JH
%  
%  LAST UPDATE : 
%   - 17/06/2008 [creation]
%  

    arch = lower(computer);