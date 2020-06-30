function [depth,dz] = taodepth
%
% set a common depth grid for all TaoTritonPirataRama temperature and
% velocity data. 
%
% Note: chipods are too sparce in depth to interpolate to
% an even depth grid, and numerous variables such as winds and surface
% fluxes do not need a depth grid.)
%
% Note: The best depth grid is slightly unusual. I want the top-most grid
% to be at 1m depth because that's where the top-most CTD is located on
% almost all TaoTritonPirataRama moorings. If the top-most bin is set to be
% at 0m, then the data from the 1m sensor is often not included in the final
% gridded product. Below that, I want an even 5-meter depth grid at 5m, 
% 10m, 15m, etc., down to 300m. 

dz = 5;
depth = (0:dz:300)';
end

