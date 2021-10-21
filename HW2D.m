%{
% 
%  DESCRIPTION:   
%       A solver for Darcy flow and particle tracking
% 
%  ASSUMPTIONS:  
%       isotropic conductivy: kx = ky = kz
% 
%  GEOMETRIC INDEXING:                                                           
%    - 3 X 3 matrix indexing example: 
%           
%       | 1  4  7 |              -----> +x 
%       | 2  5  8 |              |
%       | 3  6  9 |             \|/
%                                +y 
%     
%    - Symbols for cell interfaces (Left,Right,Up,Down):
%    
%               4                      
%           
%           --- U ---               
%           |       |                                              
%      2    L   5   R    8     
%           |       |            
%           --- D ---                                             
%                                 
%               6                                                      
% 
%}

clear all; clc; 
%% ------------------------------------------------------------------------
% Grid structure
% -------------------------------------------------------------------------
nx = 100; ny = 50;      % number of grid cells in each direction
Lx = 1; Ly = 0.5;       % [m] fracture size
dx = Lx/nx; dy = Ly/ny; % [m] cell size


% --- coordinates of cell centers
xCells = [dx/2:dx:Lx]'; 
yCells = [dy/2:dy:Ly]'; 


%% ------------------------------------------------------------------------
% Boundary Conditions
% -------------------------------------------------------------------------
Q_ = 1e-3; % [m/s] flux boundary at inlet 
h_ = 1;    % [m]   head boundary at outlet


%% ------------------------------------------------------------------------
% Physical Constants
% -------------------------------------------------------------------------
mu = 10^-3;     % viscosity [kg/m/s]
g = 9.8;        % gravitational acceleration [m/s^2]. 
rho = 1000;     % fresh water density [kg/m^3]
phi = 0.3;      % porosity

% heterogeneous k field ----------------------------------
load(['Kmat_var1.mat'],'Kmat')
Kmat = reshape(Kmat(:,5),50,500);
kMat = Kmat(1:ny,1:nx);

kMatMean = log(1e-2);%mean(log(kMat(:)));
kMatStd = std(log(kMat(:)));
kStd = 1;
kMat = exp((log(kMat)-kMatMean) / kMatStd * kStd + kMatMean);

%% ------------------------------------------------------------------------
% transmissibility matrix
% -------------------------------------------------------------------------
Ty = 2./( 1./kMat(2:ny,:) + 1./kMat(1:ny-1,:) ) /dy;
Tx = 2./( 1./kMat(:,2:nx) + 1./kMat(:,1:nx-1) ) /dx;

TU = cat(1, zeros(1,nx), Ty); % no flow B.C. at the top 
TD = cat(1, Ty, zeros(1,nx)); % no flow B.C. at the bottom 

TL = cat(2, zeros(ny,1), Tx);  
TR = cat(2, Tx, 2*dy/dx*kMat(:,nx,:));  % hydrostatic B.C. at the outlet, assuming inf Tx at outlet

TC = TL(:) + TR(:) + TU(:) + TD(:) ;
Tmat = spdiags([-TR(:) -TD(:) TC -TU(:) -TL(:)], [-ny,-1,0,1,ny], ny*nx,ny*nx)';

%% ------------------------------------------------------------------------
% Flow Simulation: head distribution computation
% -------------------------------------------------------------------------
b = zeros(ny,nx);

b(:,1) = b(:,1) + Q_*dy;                   % flux boundary at inlet
b(:,end) = b(:,end) + TR(:,end).* (h_);    % head boundary at outlet

pvec = Tmat\b(:);

pmat=reshape(pvec,ny,nx);

%% ------------------------------------------------------------------------
% Flow Simulation: velocity filed computation
% -------------------------------------------------------------------------
clc
Pzdif=[pmat(2:end,:)-pmat(1:end-1,:)];
Pxdif=[pmat(:,2:end)-pmat(:,1:end-1)];


Pxdif_L = cat(2, zeros(ny,1), Pxdif);        % flux inlet B.C.

Pxdif_R = cat(2, Pxdif, h_-pmat(:,end)); % pressure outlet B.C.
Pzdif_U = cat(1, zeros(1,nx), Pzdif);
Pzdif_D = cat(1, Pzdif, zeros(1,nx));

U_L = -TL/dy.*Pxdif_L;     
if Q_ ~= 0 % if flux inlet B.C.
    U_L(:,1) = Q_;
end
U_R = -TR/dy.*Pxdif_R;
U_U = -TU/dx.*Pzdif_U;
U_D = -TD/dx.*Pzdif_D;

U = (U_L + U_R)/2;
V = (U_U + U_D)/2;
    
% ==== visualizaiton
clf;
colormap('jet')

subplot(3,4,1:2); imagesc(xCells,yCells,log10(kMat)); axis equal tight;
colorbar;
title('log(k)'); xlabel('x'); ylabel('z');

subplot(3,4,3:4); imagesc(xCells,yCells,pmat); axis equal tight;
colorbar;
title('pressure field'); xlabel('x'); ylabel('z')

subplot(3,4,5:12); 
quiver(xCells,yCells,U,V,4); axis equal tight;
set(gca,'ydir','reverse')
xlim([0 Lx]); ylim([0 Ly])
xlabel('x'); ylabel('z')
title('flow field')

set(gcf,'color','w')

%% ------------------------------------------------------------------------
% Particle transport 
% -------------------------------------------------------------------------
% --- physical constants
Disp = 10^-7;       % dispersion coefficient

% --- coordinates of cell faces
xFaces = [xCells-dx/2; xCells(end)+dx/2];
yFaces = [yCells-dy/2; yCells(end)+dy/2];

nptcl = 100000;

% inlet particles (randomly distributed) 
xp = zeros(nptcl,1);
yp = linspace(dy/2,Ly-dy/2,nptcl)';

t = 0;
tObs = [1:0.5:50]'*60; it = 1;

% cl = [-16.44 -1.32];
while t<tObs(end)
    
    % advection
    [ixC, iyC] = sub_findNearest(xp,yp,xFaces,yFaces); % identify the nearest cell & face

    indP = (ixC-1)*ny + iyC;
    ptclU = U(indP);
    ptclV = V(indP);

    doPlot = false;
    dtAdv = 0.5*min(dx / max(abs(ptclU(:))), dy / max(abs(ptclV(:))));
    if t+dtAdv >= tObs(it)
        dtAdv = tObs(it) - t;
        doPlot = true;
        it = it+1;
    end
    dtDiff = 0.5*dy*dx/Disp;
    
    dt = min(dtAdv,dtDiff);
    
    
    t = t+dt;
    xp = xp + dt*ptclU + sqrt(2*Disp*dt)*randn(length(xp),1);
    yp = yp + dt*ptclV + sqrt(2*Disp*dt)*randn(length(yp),1);
    
    ipOutL = xp <= xFaces(1);
    ipOutR = xp >= xFaces(end);
    ipOutB = yp <= yFaces(1);
    ipOutU = yp >= yFaces(end);
    
    yp(ipOutB) = yFaces(1)   + abs(yp(ipOutB));
    yp(ipOutU) = 2*yFaces(end) - abs(yp(ipOutU));
    xp(ipOutL) = xFaces(1);
    xp(ipOutR) = []; yp(ipOutR) = []; 
    
    if doPlot
        clf
        imagesc(xCells,yCells,log10(kMat)); axis equal tight; hold on;
        colormap('jet')
        colorbar; 
        if it-1 == 1
            figName = [sprintf('kStd%.1f_it%d.jpg', kStd,0)];
            title(sprintf('time = %.1f [day]', 0))
            print( figName, '-djpeg','-r0');
        end
        scatter(xp,yp,'.k')
        title(sprintf('time = %.1f [min]', t/60))
        drawnow;
        doPlot = false;
        
%         figName = [sprintf('kStd%.1f_it%d.jpg', kStd,it-1)];
%         print( figName, '-djpeg','-r0');
%         endIt = it - 1;
    end
end

% videoName = sprintf('kStd%.1f.avi', kStd);
% outputVideo = VideoWriter(videoName);
% outputVideo.FrameRate = 5;
% open(outputVideo)
% for it = 0:endIt
%     figName = [sprintf('kStd%.1f_it%d.jpg', kStd,it)];
%     img = imread(figName);
%     writeVideo(outputVideo,img);
% end
% close(outputVideo);
