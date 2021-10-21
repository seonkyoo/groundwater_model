%{
%  Developed by:        
%       Seonkyoo Yoon < yoonx213@umn.edu >                     
% 
%  DESCRIPTION:   
%       A solver for groundwater flow and contaminant transport
% 
%  ASSUMPTIONS:  
%       isotropic permeability: kx = ky = kz
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

xCells = [dx/2:dx:Lx]';
yCells = [dy/2:dy:Ly]';
zCells = [dz/2:dz:Lz]';

[xxCells,zzCells,yyCells]=meshgrid(xCells,zCells,yCells);

%% ------------------------------------------------------------------------
% B.C
% -------------------------------------------------------------------------

dP = 1e+4;
Uin = 0;%3*1e-5;

%% ------------------------------------------------------------------------
% Physical Constants
% -------------------------------------------------------------------------
D = 1;              % dispersivity [m]

visco = 10^-3;         % viscosity [kg/m/s]
g = 0%9.8;         % gravitational acceleration [m/s^2]. 
                       % 0: horizontal | 9.8: vertical
rho0 = 1000;            % fresh water density [kg/m^3]
DelRho = 30;             % TCE solubility into water [kg/m^3]
rho1 = rho0+DelRho;     % saline water density [kg/m^3]
phi = 0.3;              % porosity

kTrueMean = exp(-23);  % mean permeability

load(['Kmat_var1.mat'],'Kmat')
Kmat = reshape(Kmat(:,1),50,500);
kMat = Kmat(1:nz,1:nx,1:ny);

% kMat(1:13,:) = kTrueMean * 10/10;
% % kMat(1:13,[11:13 21:23 31:33 41:43 51:53 61:63 71:73 81:83 91:93]) = kTrueMean * 1000/10;
% kMat(14:16,:) = kTrueMean * 1000/10;
% kMat(17:30,:) = kTrueMean * 20/10;
% kMat(31:50,:) = kTrueMean * 5/10;

% kMat(:) = kTrueMean ;
% kMat(11:40,41:60) = kTrueMean * 1000/10;
% kMat(21:30,21:80) = kTrueMean * 1000/10;
% kMat(14:16,:) = kTrueMean * 1000/10;
% kMat(17:30,:) = kTrueMean * 20/10;
% kMat(31:50,:) = kTrueMean * 5/10;

imagesc(zCells,xCells,log(kMat))
clf; imagesc(log(kMat))

axis equal tight
%% ------------------------------------------------------------------------
% transmissibility matrix
% -------------------------------------------------------------------------
lMat = kMat./visco./phi; % mobility matrix

Tz = 2./( 1./lMat(2:nz,:,:) + 1./lMat(1:nz-1,:,:) ) * dx*dy/dz;
Tx = 2./( 1./lMat(:,2:nx,:) + 1./lMat(:,1:nx-1,:) ) * dy*dz/dx;
Ty = 2./( 1./lMat(:,:,2:ny) + 1./lMat(:,:,1:ny-1) ) * dz*dx/dy;

T.U = cat(1, zeros(1,nx,ny), Tz); TU=T.U(:); % no flow B.C. at the top 
T.D = cat(1, Tz, zeros(1,nx,ny)); TD=T.D(:); % no flow B.C. at the bottom 

T.L = cat(2, zeros(nz,1,ny), Tx);  TL=T.L(:); 
T.R = cat(2, Tx, zeros(nz,1,ny));  TR=T.R(:); 

T.F = cat(3, zeros(nz,nx,1), Ty); TF=T.F(:); % no flow B.C. at the front 
T.B = cat(3, Ty, zeros(nz,nx,1)); TB=T.B(:); % no flow B.C. at the back 

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
tObs = linspace(0,tEnd,100);

it = 1; 
doPlot = false;
% c0 = zeros(nz,nx,ny);  

[xx,~,~] = meshgrid(dx/2:dx:Lx-dx/2,dz/2:dz:Lz-dz/2,dy/2:dy:Ly-dy/2);
c0= 1-(1/2*( tanh(1*(xx-0.01*Lx)) )+0.5); % smooth initial concentration profile

t = 0;
while t < tEnd
%     c0(:,1,:) = 1;
%     c0(1020) = 1;
    %%  pressure field
    if g ~= 0 || ~exist('pmat','var')
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
        
    end
        
    %% Transport (advection)
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

    c0_L = cat(2, ones(nz,1,ny), c0(:,1:end-1,:));  
    c0_R = cat(2, c0(:,2:end,:), zeros(nz,1,ny));    
    c0_U = cat(1, zeros(1,nx,ny), c0(1:end-1,:,:));       
    c0_D = cat(1, c0(2:end,:,:), zeros(1,nx,ny));
    c0_F = cat(3, zeros(nz,nx,1), c0(:,:,1:end-1));
    c0_B = cat(3, c0(:,:,2:end), zeros(nz,nx,1));

    e0=zeros(nz,nx,ny); 
    Fadv_L = e0; Fadv_R = e0; 
    Fadv_D = e0; Fadv_U = e0; 
    Fadv_F = e0; Fadv_B = e0; 

    Fadv_L(U_L>=0) = U_L(U_L>=0) .* c0_L(U_L>=0)*dz*dy; 
    Fadv_L(U_L<0)  = U_L(U_L<0)  .* c0(U_L<0)   *dz*dy; 
    Fadv_R(U_R>=0) = U_R(U_R>=0) .* c0(U_R>=0)  *dz*dy;
    Fadv_R(U_R<0)  = U_R(U_R<0)  .* c0_R(U_R<0) *dz*dy;
    Fadv_U(U_U>=0) = U_U(U_U>=0) .* c0_U(U_U>=0)*dx*dy; 
    Fadv_U(U_U<0)  = U_U(U_U<0)  .* c0(U_U<0)   *dx*dy; 
    Fadv_D(U_D>=0) = U_D(U_D>=0) .* c0(U_D>=0)  *dx*dy;
    Fadv_D(U_D<0)  = U_D(U_D<0)  .* c0_D(U_D<0) *dx*dy;
    Fadv_F(U_F>=0) = U_F(U_F>=0) .* c0_F(U_F>=0)*dz*dx;
    Fadv_F(U_F<0)  = U_F(U_F<0)  .* c0(U_F<0)   *dz*dx;
    Fadv_B(U_B>=0) = U_B(U_B>=0) .* c0(U_B>=0)  *dz*dx;
    Fadv_B(U_B<0)  = U_B(U_B<0)  .* c0_B(U_B<0) *dz*dx;

    Fadv = (-Fadv_F -Fadv_L - Fadv_U + Fadv_D + Fadv_R + Fadv_B);    

    %% Transport (Dispersion)
    Fdif_L = (10^-9+D*abs(U_L)).*(c0 - c0_L)/dx*dz*dy;
    Fdif_R = (10^-9+D*abs(U_R)).*(c0_R - c0)/dx*dz*dy;
    Fdif_U = (10^-9+D*abs(U_U)).*(c0 - c0_U)/dz*dx*dy;
    Fdif_D = (10^-9+D*abs(U_D)).*(c0_D - c0)/dz*dx*dy; 
    Fdif_F = (10^-9+D*abs(U_F)).*(c0 - c0_F)/dy*dz*dx; 
    Fdif_B = (10^-9+D*abs(U_B)).*(c0_B - c0)/dy*dz*dx; 

    Fdif = (-Fdif_F -Fdif_L - Fdif_U + Fdif_D + Fdif_R + Fdif_B); 

    %% Advance in Time
    umax = max([...
        max(max(max(abs(cat(3, U_F, U_B(:,:,end)))))), ...
        max(max(max(abs(cat(2, U_L, U_R(:,end,:)))))), ...
        max(max(max(abs(cat(1, U_U, U_D(end,:,:))))))  ...
        ]);
    dt_adv=.1*dy/umax;
    
    Dmax = max([...
        max(max(max(10^-9+D*abs(cat(3, U_F, U_B(:,:,end)))))), ...
        max(max(max(10^-9+D*abs(cat(2, U_L, U_R(:,end,:)))))), ...
        max(max(max(10^-9+D*abs(cat(1, U_U, U_D(end,:,:))))))  ...
        ]);
    dt_diff=.1*dx*dx/Dmax; % diffusive time scale

    dt=min(dt_adv,dt_diff);

    if t+dt >= tObs(it)
        dt = tObs(it) - t;
        doPlot = true;
        it = it+1;
    end

    t=t+dt;
    c0 = c0 - dt/(dx*dz*dy)*(Fadv-Fdif);

    if doPlot 
        clf;
        subplot(3,1,1); imagesc(xCells,zCells,log10(kMat)); axis equal tight;
        subplot(3,1,2); imagesc(xCells,zCells,pmat); axis equal tight;
        colorbar;
        subplot(3,1,3); 
        imagesc(xCells,zCells,c0); axis equal tight; colorbar;
        caxis([0 1])
        title(sprintf('time: %.1f day',t/60/60/24))
        drawnow;

        doPlot = false;
    end
    
end
