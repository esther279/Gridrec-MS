% taperfunc = radtap(X,Y,tappix,zerorad), Creates a central cosine tapering for
% beamstop. Recieves the X and Y coordinates, tappix is the extent of
% tapering, zerorad is the radius with no data (zeros).

function taperfunc = radtap(X,Y,tappix,zerorad),

tau = 2*tappix; % period of cosine function (only half a period is used)

R = sqrt(X.^2+Y.^2);
% taperfunc = (R)>=10;
taperfunc = 0.5*(1+cos(  2*pi*(R-zerorad-tau/2)/tau  ));
taperfunc = (R>zerorad+tau/2)*1.0 + taperfunc.*(R<=zerorad+tau/2);
taperfunc = taperfunc.*(R>=zerorad);