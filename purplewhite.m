function c=purplewhite(m)
 
% PURPLE-TO-WHITE Linear green to red colormap
%
% usage:
%  redgreen
%  redgreen(M)
%
% generates an Mx3 matrix representing the colormap.  If no input argument
% is supplied, M is set to the length of the current colormap.
%
% See also COLORMAP, JET, YELLOWBLUE
 
% by Alex Rosenberg
% 2-FEB-2016
% Copyright (c) 2011. All rights reserved.
% This software is offered with no guarantees of any kind.
% afr@bioinfx.com
 
 
if nargin<1; m=size(get(gcf,'colormap'),1); end         % get size of current colormap
if mod(m,2); m=m-1; end                                 % size must be even number
c2=((1:m)/m)';
c13=(.4:.6*(1/m):1)';
c13 = c13(2:end);
c = [c13 c2 c13];
 
return
 
 
