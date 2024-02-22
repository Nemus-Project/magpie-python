import scipy
def magpie(rho,E,nu,ldim,h,BCs,Nmodes,plot_type):
	##----------------------------
    D = E * (Lz ** 3) / 12 / (1-(nu ** 2));
    Nx      = floor(Lx/h) ;
    Ny      = floor(Ly/h) ;
    ##----------------------------
    ## Build BiHarmonic
    biHarm = bhmat(BCs,[Nx Ny], h, Lz, E, nu);
    
    ## EIGENVALUES    
    [Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
    _,indSort = sort(diag((Dm))) ;
    Q = Q[:,indSort];
    
    
    freqs = np.sqrt(np.abs(diag(Dm)))*np.sqrt(D/rho/Lz)/2/pi ;
    Om    = 2*np.pi*freqs ;
    
	if plot_type == 'chladni'
		subs = ceil(sqrt(Nmodes));
		
		colormap('copper') ;
		cmp = colormap;
		cmp = flipud(cmp);
		colormap(cmp);    
			
		for m in np.r_[1 : Nmodes]:
			mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
			subplot(subs,subs,m)
			mesh(3e3*real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ; 
			view(2); 
			axis equal; 
			axis tight;        
			axis off;
			clim([0.00005 0.002]);        
		

	elif plot_type == '3D':
        
            subs = ceil(sqrt(Nmodes));
            colormap('parula') ;

            for m in np.r_[1 : Nmodes]:
                mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
                subplot(subs,subs,m)
                mesh(3000*(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;

