function [Q,Om,Nx,Ny,biHarm] = biharmonicEigs(rho,E,nu,ldim,h,BCs,Nmodes,plot_type)
% BIHARMONICEIGS What does this do?
%   [Q,Om,Nx,Ny,biHarm] = BIHARMONICEIGS (rho,E,nu,Lx,Ly,Lz,h,K0y,R0y,Kx0,Rx0,KLy,RLy,KxL,RxL,Nmodes)
%   A function that returns:
%           Q       : Eigen vector(s)
%           Om      : Angular modal frequencies
%           Nx      : Grid points along the x-axis
%           Ny      : Grid points along the y-axis
%           biHarm  : Biharmonic Matrix for the plate
%
%   Arguments:
%       rho          %-- density [kg/m^3]
%       E            %-- Young's mod [Pa]
%       nu           %-- poisson's ratio
%       
%       3 element array  representing  x,yz dimensions of plate
%       ldim = [Lx,  %-- length along x [m]
%               Ly,  %-- length along y [m]
%               Lz]  %-- thickness [m]
%
%       h            %-- grid spacing
%       
%       2 column array of eleastic boundary constants around each edge of
%       the plate.
%
%       Column 1 repesents  ...
%       Column 2 repesents  ...
%
%       BCs = [K0y, R0y;
%              Kx0, Rx0;
%              KLy, RLy;
%              KxL, RxL];
%
%       Nmodes      %-- number of modes to compute
%
%
%       Example:
%
%           %% physical and elastic parameters
%           Lx = 0.10; Ly = 0.08; Lz = 0.81e-3;
%           ldim = [Lx Ly Lz];   % plate dimensions [x, y, z] in metres
%           E       = 1.01e+11 ;        %-- Young's mod [Pa]
%           rho     = 8765 ;            %-- density [kg/m^3]
%           nu      = 0.3 ;             %-- poisson's ratio
%           Nmodes  = 16;               %-- number of modes to compute
%           h       = sqrt(Lx*Ly)*0.01; %--           
%           BCs = ones(4,2) * 1e15      %-- elastic constants around the edges
%           
%           [Q,Om,Nx,Ny,biHarm] = biharmonicEigs(rho, E, nu, ldim, h, BCs, Nmodes);
    %% Validation        
    validateattributes(rho,      {'double'}, {'nonempty'});
    validateattributes(E,        {'double'}, {'nonempty'});
    validateattributes(nu,       {'double'}, {'nonempty'});
    validateattributes(ldim,     {'double'}, {'numel', 3});
    validateattributes(h,        {'double'}, {'nonempty'});    
    validateattributes(BCs,      {'double'}, {'size', [4,2]});
    validateattributes(Nmodes,   {'numeric'}, {'integer','positive'});
    validatestring(plot_type,["chladni","3D","none"]);
              
    %% Unpack array variables
    pack_ldim = num2cell(ldim);
    pack_BCs = num2cell(BCs);
    [Lx, Ly, Lz] = pack_ldim{:};
    [K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};

    %%--- derived parameters (don't change here)    
    D = E * Lz^3 / 12 / (1-nu^2);
    Nx      = floor(Lx/h) ;
    Ny      = floor(Ly/h) ;
    %%----------------------------
    %% Build BiHarmonic
    biHarm = bhmat(BCs,[Nx Ny], h, Lz, E, nu);
    
    %% EIGENVALUES

    [Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
    [~,indSort] = sort(diag((Dm))) ;
    Q = Q(:,indSort) ;


    freqs = sqrt(abs(diag(Dm)))*sqrt(D/rho/Lz)/2/pi ;
    Om    = 2*pi*freqs ;


    switch plot_type
        case 'chladni'
            subs = ceil(sqrt(Nmodes));
            
            colormap('copper') ;
            cmp = colormap;
            cmp = flipud(cmp);
            colormap(cmp);    
                
            for m = 1 : Nmodes
                mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
                subplot(subs,subs,m)
                mesh(3e3*real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ; 
                view(2); 
                axis equal; 
                axis tight;        
                axis off;
                clim([0.00005 0.002]);        
            end

        case '3D'
        
            subs = ceil(sqrt(Nmodes));
            colormap('parula') ;

            for m = 1 : Nmodes
                mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
                subplot(subs,subs,m)
                mesh(3000*(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
            end 
    end
    
end
