%
% Filename: $RCSfile: franzmap.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Function FM = FRANZMAP(M)
% Franz's modified jet color map
% This function returns the colormap and should
% be used like any other colormaps, e.g.
% colormap(franzmap(128));
% or
% imwrite(uint8(255*myarray/max(max(myarray))),franzmap(256),'myarray.jpg');
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% none
%
% history:
%
% May 9th 2008, Oliver Bunk: add CVS header
%
% Franz Pfeiffer: 1st version
%
function fm = franzmap(m)

if (nargin < 1)
    m = size(get(gcf,'colormap'),1);
end

maxval = ceil(m/8);
fm = jet(m);
fm(1:maxval,3) = [1:maxval]/maxval;
