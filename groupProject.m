%{
%  Developed by:        
%       Seonkyoo Yoon < yoonx213@umn.edu >                     
% 
%  DESCRIPTION:   
%       A solver for Darcy flow and particle tracking
% 
%  ASSUMPTIONS:  
%       isotropic permeability: kx = ky = kz
%       anisotropic permeability: kx = ky = 10*kz
% 
%  GEOMETRIC INDEXING:                                                           
%    % The indice are arranged like this:                             
%    %         IndU                   +y (j)                             
%    %           | IndB               /                               
%    %           |/                  /                                
%    %  IndL -- Ind -- IndR         /---- +x (i)                        
%    %          /|                  |                                 
%    %      IndF |                  |                                 
%    %         IndD                +z (k)                              
%    %                                                                
%    % priority: z > x > y   
%    %                                                               
%    % example:                                                      
%    %              ____19 22 25                                    
%    %       ____10 13 16 |23 26                                    
%    %     1  4  7 |14 17 |24 27                                    
%    %     2  5  8 |15 18 |                                         
%    %     3  6  9 |                                                
% 
%}

clear all; clc; 


%% ------------------------------------------------------------------------
% Grid structure
% -------------------------------------------------------------------------
dx = 1; dy = 1; dz = 1; %[m]

nx = 100; ny = 1; nz = 50;   % 2d

Lx = nx*dx; Ly = ny*dy; Lz = nz*dz;

nc = nz*nx*ny;  % number of grid cells 

xCells = [dx/2:dx:Lx]'; xFaces = [xCells-dx/2; xCells(end)+dx/2];
yCells = [dy/2:dy:Ly]'; yFaces = [yCells-dy/2; yCells(end)+dy/2];
zCells = [dz/2:dz:Lz]'; zFaces = [zCells-dz/2; zCells(end)+dz/2];

[xxCells,zzCells,yyCells]=meshgrid(xCells,zCells,yCells);

%% ------------------------------------------------------------------------
% B.C
% -------------------------------------------------------------------------

dP = 1e+4; Uin = 0;   % pressure inlet & outlet
% dP = 0; Uin = 3*1e-5; % fluxt inlet & pressure outlet

%% ------------------------------------------------------------------------
% Physical Constants
% -------------------------------------------------------------------------
D = 1;              % dispersivity [m]

visco = 10^-3;         % viscosity [kg/m/s]
g = 0;%9.8;         % gravitational acceleration [m/s^2]. 
                       % 0: horizontal | 9.8: vertical
rho0 = 1000;            % fresh water density [kg/m^3]
DelRho = 30;             % TCE solubility into water [kg/m^3]
rho1 = rho0+DelRho;     % saline water density [kg/m^3]
phi = 0.3;              % porosity


% stratigraphic k field ----------------------------------
MN = 1:13; BPP = 14:15; HF = 16:30; MF = 31:50;

zIndFracMN = MN; zIndFracMN(1)=[]; zIndFracMN(end)=[];
zIndFracHF = HF; zIndFracHF(1)=[]; zIndFracHF(end)=[];
zIndFracMF = HF; zIndFracMF(1)=[]; zIndFracMF(end)=[];

anisotropicRateMN = 3;
anisotropicRateBPP = 1;
anisotropicRateHF = 5;
anisotropicRateMF = 10;

kTrueMean = exp(-23);  % mean permeability
kMat = kTrueMean*ones(nz,nx);
kMat(MN,:) =  kTrueMean * 10/10;
kMat(BPP,:) = kTrueMean * 1000/10;
kMat(HF,:) = kTrueMean * 20/10;
kMat(MF,:) = kTrueMean * 5/10;


kMat(zIndFracMN,[1:3 11:13 21:23 31:33 41:43 51:53 61:63 71:73 81:83 91:93]+2) = ...
                kTrueMean * 1000/10;

% kMat(HF,[1:3 11:13 21:23 31:33 41:43 51:53 61:63 71:73 81:83 91:93]+3) = ...
%                 kTrueMean * 1000/10;

kMat(zIndFracHF,[ 11:13  31:33  51:53  71:73  91:93]+5) = ...
                kTrueMean * 1000/10;

% kMat(MF,[11:13 21:23 31:33 41:43 51:53 61:63 71:73 81:83 91:93]) = ...
%                 kTrueMean * 1000/10;

kMat(MF,[11:13 31:33 51:53 71:73 91:93]) = ...
                kTrueMean * 1000/10;

% kMat(10:20,[11:13 21:23 31:33 41:43 51:53 61:63 71:73 81:83 91:93]) = ...
%                 kTrueMean * 1000/10;


% kMat(11:40,41:60) = kTrueMean * 1000/10;
% kMat(20:30,20:80) = kTrueMean * 1000/10;


% heterogeneous k field ----------------------------------
load(['Kmat_var1.mat'],'Kmat')
Kmat = reshape(Kmat(:,1),50,500);
kMat = Kmat(1:nz,1:nx,1:ny);

kMatMean = mean(log(kMat(:)));
kMatStd = std(log(kMat(:)));
kMat = exp((log(kMat)-kMatMean) / kMatStd * 0.25 + kMatMean);

anisotropicRateMN = 1;
anisotropicRateBPP = 1;
anisotropicRateHF = 1;
anisotropicRateMF = 1;

% imagesc(zCells,xCells,log(kMat))
% clf; imagesc(log(kMat))
% colorbar;
% 
% 
% axis equal tight
%% ------------------------------------------------------------------------
% transmissibility matrix
% -------------------------------------------------------------------------
lMat = kMat./visco./phi; % mobility matrix

Tz = 2./( 1./lMat(2:nz,:,:) + 1./lMat(1:nz-1,:,:) ) * dx*dy/dz;
Tx = 2./( 1./lMat(:,2:nx,:) + 1./lMat(:,1:nx-1,:) ) * dy*dz/dx;
Ty = 2./( 1./lMat(:,:,2:ny) + 1./lMat(:,:,1:ny-1) ) * dz*dx/dy;

T.U = cat(1, zeros(1,nx,ny), Tz); % no flow B.C. at the top 
T.D = cat(1, Tz, zeros(1,nx,ny)); % no flow B.C. at the bottom 

T.U(MN,:) = T.U(MN,:) / anisotropicRateMN;
T.U(BPP,:) = T.U(BPP,:) / anisotropicRateBPP;
T.U(HF,:) = T.U(HF,:) / anisotropicRateHF;
T.U(MF,:) = T.U(MF,:) / anisotropicRateMF;

T.D(MN,:) = T.D(MN,:) / anisotropicRateMN;
T.D(BPP,:) = T.D(BPP,:) / anisotropicRateBPP;
T.D(HF,:) = T.D(HF,:) / anisotropicRateHF;
T.D(MF,:) = T.D(MF,:) / anisotropicRateMF;

T.L = cat(2, zeros(nz,1,ny), Tx);  
T.R = cat(2, Tx, zeros(nz,1,ny));  

T.F = cat(3, zeros(nz,nx,1), Ty); % no flow B.C. at the front 
T.B = cat(3, Ty, zeros(nz,nx,1)); % no flow B.C. at the back 

TU=T.U(:); TD=T.D(:); TL=T.L(:); TR=T.R(:); TF=T.F(:); TB=T.B(:); 

if Uin == 0 && dP ~= 0 % pressure inlet B.C., assuming inf Tx at inlet
    T.L = cat(2, 2*dz*dy/dx*lMat(:,1,:), Tx);
end
T.R = cat(2, Tx, 2*dz*dy/dx*lMat(:,nx,:));  % hydrostatic B.C. at the outlet, assuming inf Tx at outlet

TC = T.L(:) + T.R(:) + T.U(:) + T.D(:) + T.F(:) + T.B(:);
Tmat = spdiags([-TB -TR -TD TC -TU -TL -TF],...
                          [-nz*nx,-nz,-1,0,1,nz,nz*nx],nz*nx*ny,nz*nx*ny)';
clear TB TR TD TC TU TL TF

[Lmat, Dmat, perm] = ldl(Tmat,'vector');
invDmat = 1./diag(Dmat);


%% ------------------------------------------------------------------------
% Flow and Transport Simulation
% -------------------------------------------------------------------------
tEnd = 60*24*60*60; %[sec]
tObs = linspace(0,tEnd,10);

it = 1; 
doPlot = false;
% c0 = zeros(nz,nx,ny);  

[xx,~,~] = meshgrid(dx/2:dx:Lx-dx/2,dz/2:dz:Lz-dz/2,dy/2:dy:Ly-dy/2);
c0= 1-(1/2*( tanh(1*(xx-0.01*Lx)) )+0.5); % smooth initial concentration profile

t = 0;
%%  pressure field
c0_avg = cat(1, zeros(1,nx,ny), (c0(1:end-1,:,:)+c0(2:end,:,:))/2, zeros(1,nx,ny));

b = + T.U.*(rho0+DelRho*c0_avg(1:end-1,:,:))*g*dz ... vertical flux
    - T.D.*(rho0+DelRho*c0_avg(2:end,:,:))*g*dz;

Phydro = repmat(rho0*g*dz*[1:nz]',1,ny); % hydrostatic pressure at the right boundary
Phydro = reshape(Phydro, nz,1,ny);

if Uin ~= 0 && dP == 0
    b(:,1,:) = b(:,1,:) + Uin*dz*dy  ;
elseif Uin == 0 && dP ~= 0
    b(:,1,:) = b(:,1,:) + T.L(:,1,:).* (Phydro+dP);
end
b(:,end,:) = b(:,end,:) + T.R(:,end,:).* (Phydro);

bvec = b(:);

pvec(perm,:)=Lmat'\(invDmat.*(Lmat\(bvec(perm,:))));

pmat=reshape(pvec,nz,nx,ny);

%% velocity filed
Pzdif=[pmat(2:end,:,:)-pmat(1:end-1,:,:)]-(rho0+DelRho*c0_avg(2:end-1,:,:))*g*dz;
Pxdif=[pmat(:,2:end,:)-pmat(:,1:end-1,:)];
Pydif=[pmat(:,:,2:end)-pmat(:,:,1:end-1)];

if Uin ~= 0 % if flux inlet B.C.
    Pxdif_L = cat(2, zeros(nz,1,ny), Pxdif);        % flux inlet B.C.
else
    Pxdif_L = cat(2, pmat(:,1,:)-(Phydro+dP), Pxdif);% pressure inlet B.C.
end
Pxdif_R = cat(2, Pxdif, Phydro-pmat(:,end,:)); % pressure outlet B.C.
Pzdif_U = cat(1, zeros(1,nx,ny), Pzdif);
Pzdif_D = cat(1, Pzdif, zeros(1,nx,ny));
Pydif_F = cat(3, zeros(nz,nx,1), Pydif);
Pydif_B = cat(3, Pydif, zeros(nz,nx,1));

U_L = -T.L/dz/dy.*Pxdif_L;     
if Uin ~= 0 % if flux inlet B.C.
    U_L(:,1) = Uin;
end
U_R = -T.R/dz/dy.*Pxdif_R;
U_U = -T.U/dx/dy.*Pzdif_U;
U_D = -T.D/dx/dy.*Pzdif_D;
U_F = -T.F/dz/dx.*Pydif_F;
U_B = -T.B/dz/dx.*Pydif_B;

U = (U_L + U_R)/2;
V = (U_U + U_D)/2;


    
%% visualizaiton
clf;
subplot(4,5,[1 2 6 7]); imagesc(xCells,zCells,log10(kMat)); axis equal tight;
colorbar;
title('log(k)'); xlabel('x'); ylabel('z')
subplot(4,5,[11 12 16 17]); imagesc(xCells,zCells,pmat); axis equal tight;
colorbar;
title('pressure field'); xlabel('x'); ylabel('z')

subplot(4,5,[3 8 13 18]);
plot(pmat(:,50),1:50);
set(gca,'ydir','reverse')
xl = xlim;
xlim([xl(1)-10 xl(2)+10])
title('x=50')
xlabel('head')
ylabel('z')

% subplot(4,,4,[3 4 7 8 11 12 15 16 ]); 
% quiver(xCells,zCells,U,V,1.5); axis equal tight;
% set(gca,'ydir','reverse')
% xlim([25 75]); ylim([0 50])

subplot(4,5,[4 9 14 19]); 
quiver(xCells,zCells,U,V,1.5); axis equal tight;
set(gca,'ydir','reverse')
xlim([45 55]); ylim([0 50])
xlabel('x'); ylabel('z')
title('flow field')

subplot(4,5,[5 10 15 20]); 
Unorm = sqrt(U.^2 + V.^2);
quiver(xCells,zCells,U./Unorm,V./Unorm); axis equal tight;
set(gca,'ydir','reverse')
xlim([45 55]); ylim([0 50])
xlabel('x'); ylabel('z')
title('flow direction')

set(gcf,'color','w')

%% particle tracking

nptcl = 100000;

% inlet particles (randomly distributed) 
xp = zeros(nptcl,1);
zp = (Lz-1)*rand(nptcl,1) + dz/2;

t = 0;
tObs = [1:40]'*24*60*60; it = 1;

% cl = [-16.44 -1.32];
while t<tObs(end)
    [ixC, izC, ixFce, izFce] = sub_findNearest(xp,zp,xFaces,zFaces); % identify the nearest cell & face

    indP = (ixC-1)*nz + izC;
    ptclU = U(indP);
    ptclV = V(indP);

    doPlot = false;
    dt = min(dx / max(abs(ptclU(:))), dz / max(abs(ptclV(:))));
    if t+dt >= tObs(it)
        dt = tObs(it) - t;
        doPlot = true;
        it = it+1;
    end
    t = t+dt;
    xp = xp + dt*ptclU;
    zp = zp + dt*ptclV;
    

    if doPlot
        clf
        imagesc(xCells,zCells,log10(kMat)); axis equal tight; hold on;
        colormap('jet')
        colorbar; 
        scatter(xp,zp,'.k')
        title(sprintf('time = %.1f [day]', t/60/60/24))
        drawnow;
        doPlot = false;
    end
end
