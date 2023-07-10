function rg=usa(m)
 
% USA Linear blue to white to red colormap
%
% usage:
%  usa
%  usa(M)
%
% generates an Mx3 matrix representing the colormap.  If no input argument
% is supplied, M is set to the length of the current colormap.
%
% See also COLORMAP, JET, YELLOWBLUE, REDGREEN
 
% by Alex Rosenberg
% 11-APR-2013
% Copyright (c) 2013. All rights reserved.
% This software is offered with no guarantees of any kind.
 
 
if nargin<1; m=size(get(gcf,'colormap'),1); end         % get size of current colormap
if mod(m,2); m=m-1; end                                 % size must be even number
blue=flipud(((1:(m/2))/(m/2))'*[1 1 0])+(ones(m/2,1)*[0 0 1]);
red=((1:(m/2))/(m/2))'*[0 1 1]+(ones(m/2,1)*[1 0 0]);
rg=flipud([red; blue]);
 
return
 
 
