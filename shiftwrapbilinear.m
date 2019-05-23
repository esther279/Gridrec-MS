% out = shiftwraplinear(in,dy,dx)
% Shifts an image with wrap around and bilinear interpolation
function out = shiftwrapbilinear(in,dy,dx)
% 
% in = flipud(im2double(imread('cameraman.tif')));
% in = in(1,:);
% %in = in*0;
% %in(100,50) = 1;
% dy = -12;
% dx = -500;

[ny nx] = size(in);
dyfloor = floor(dy);
dxfloor = floor(dx);

x = [1:nx]+dxfloor;
y = [1:ny]+dyfloor;

% Shift integer step (floor of desired)
y = mod(y-1,ny)+1;
x = mod(x-1,nx)+1;

out = in(y,x);


% Subpixel (bilinear)
taux = dx-dxfloor;
tauy = dy-dyfloor;
if (taux~=0)||(tauy~=0)
    indx    = [1:nx];
    indxp1  = [2:nx 1]; 
    indy    = [1:ny];
    indyp1  = [2:ny 1];
    out =   out(indy,indx)*(1-tauy)*(1-taux) + ...
            out(indyp1,indx)*tauy*(1-taux) + ...
            out(indy,indxp1)*(1-tauy)*taux + ...
            out(indyp1,indxp1)*tauy*taux;
end



% figure(100);
% imagesc(out);
% axis xy equal tight;
% colorbar
% figure(101);
% imagesc(real(shiftpp2(in,dy,dx)));
% axis xy equal tight;
% colormap gray
% colorbar