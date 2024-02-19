function FDMat = fidimat(arg1,arg2,arg3,arg4)
  % FIDIMAT Generate Finite Difference Spatial Sparcity Matrices
  %   FDMat = FIDIMAT(l,ord) Generate stencil FDMat with order given by ord 
  %           for a 1D system length l. Boundary condition defaults to simply supported
  %
  %   FDMat = FIDIMAT(l,m,ord) Generate stencil FDMat with order given by ord 
  %           for a 2D system of size [l m]. Boundary condition defaults to simply supported
  %
  %   FDMat = FIDIMAT(l,m,ord,bctype) Generate stencil FDMat with order given 
  %           by ord for a 2D system of size [l m] with specified boundary conditions bctype.
  %
  %
  %         m       % number of total grid points X axis
  %         l       % number of total grid points Y axis
  %         ord     % order of the matrix (string)
  %         bctype  % boundary condition type: 1: simply supported, 2: clamped
  %
  %                 % Valid order inputs
  %                   ['x-','x+','x.','xx','xxxx',
  %                   'y-','y+','y.','yy','yyyy',                  
  %                   'grad','xy','xxyy','laplace','biharm','I'];

  %% Variable Arguement Length Check
  if nargin<2
    error('Not enough input arguments')

  elseif nargin==2     
      ord = convertCharsToStrings(arg2);   
      l = arg1;
      m = 1;
      bctype = 1;
  elseif nargin==3

    if ischar(arg2) || isstring(arg2)
      ord = convertCharsToStrings(arg2);
      l = arg1;
      m = 1;
      bctype = arg3;
    else
      ord = convertCharsToStrings(arg3);
      l = arg1;
      m = arg2;
      bctype = 1;  
    end
  elseif nargin==4
    l = arg1;
    m = arg2;
    ord = convertCharsToStrings(arg3);
    bctype = arg4;

  elseif nargin>1
    error('Too many input arguments')
  end

  %% Validate Arguments
  valid_orders = ["x-","x+","x.","xx","xxxx", ... 
                  "y-","y+","y.","yy","yyyy", ...
                  "grad","xy","xxyy","laplace","biharm","I"];

  validateattributes(l,      {'numeric'}, {'integer','positive'});
  validateattributes(m,      {'numeric'}, {'integer','positive'});
  ord = validatestring    (ord,    valid_orders);
  validateattributes(bctype, {'numeric'}, {'integer','positive', '<', 3});

  %% Set Identity Matrices

  Iy = speye(l); % identity matrix for y axis
  Ix = speye(m); % identity matrix for x axis
  ss = l*m;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 1D Case
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % if only one dimension is stipulated the function will return only
  % a 1D FD matrix. In this case input is by the letter x

  if m == 1

    switch ord

    case 'x-'
      FDMat = spdiags([-ones(l,1),ones(l,1)],-1:0,speye(l));

    case 'x+'
      FDMat = spdiags([-ones(l,1),ones(l,1)],0:1,speye(l));

    case 'x.'
      FDMat = spdiags([-ones(l,1),zeros(l,1),ones(l,1)],[-1:1],speye(l));

    case 'xx'

      XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
      % if bctype == 0
      %   XX([1 end],:) = 0;
      % end
      FDMat = XX;

    case 'xxxx'
      XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
      FDMat = XX;
      FDMat = XX^2;

    case 'I'
      I = speye(ss);
      I([1 end],:) = 0;
      FDMat = I;

    case {"grad","xy","xxyy","laplace","biharm"}
        error('%s is not a valid order for a 1D scheme', ord);
    otherwise
      error('something went wrong, check your arguments');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2D Case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    % if two dimension arguements are given then the matri will out put the
    % corresponding 2D Matrix.
    %
    % y is for columns of a matrix whilst x is for rows.
    % x matrices in this case will contain 'off centre' diagonals corresponding to
    % the number of elements in the row.

    switch ord

    case 'y-'
      Y = spdiags([-ones(l,1),ones(l,1)],-1:0,speye(l));
      FDMat = kron(Ix,Y);

    case 'y+'
      Y = spdiags([ones(l,1),-ones(l,1)],0:1,speye(l));
      FDMat = kron(Ix,Y);

    case 'y.'
      Y = spdiags([-ones(l,1),zeros(l,1),ones(l,1)], -1:1 ,speye(l));
      FDMat = kron(Ix,Y);

    case 'yy'
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
      FDMat = kron(Ix,YY);

    case 'yyyy'
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
      FDMat = kron(Ix,YY)^2;

    case 'x-'
      X = spdiags([-ones(m,1),ones(m,1)],-1:0,speye(m));
      FDMat = kron(X,Iy);

    case 'x+'
      X = spdiags([ones(m,1),-ones(m,1)],0:1,speye(m));
      FDMat = kron(X,Iy);

    case 'x.'
      X = spdiags([-ones(m,1),zeros(m,1),ones(m,1)],[-1:1],speye(m));
      FDMat = kron(X,Iy);

    case 'xx'
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
      FDMat = kron(XX,Iy);

    case 'xxxx'
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
      FDMat = kron(XX,Iy)^2;
    
    case {'grad', 'xy'}
      Y = spdiags([-ones(l,1),zeros(l,1),ones(l,1)],[-1:1],speye(l));
      X = spdiags([-ones(m,1),zeros(m,1),ones(m,1)],[-1:1],speye(m));
      FDMat = kron(X,Iy) + kron(Ix,Y);

    case 'xxyy'
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
      YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;

      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,speye(m));
      XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;

      FDMat = kron(Ix,YY)*kron(XX,Iy);

    case 'laplace'

      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
      XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
      YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
      LA = kron(XX,Iy) + kron(Ix,YY);
      FDMat = LA;

    case 'biharm'

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Building Bi-Harmonic
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
      XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
      YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
      YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
      LA = kron(XX,Iy) + kron(Ix,YY);
      BH = LA*LA;
      FDMat = BH;

    case 'I'
      I = speye(ss);
      FDMat = I;

    otherwise
      error('something went wrong, check your arguements');

    end

    % set the correct points to 0
    FDMat(:,[1:l+1, ss-l:ss]) = 0;FDMat([1:l+1, ss-l:ss],:) = 0;
    FDMat(:,[2*l:l:ss, (2*l)+1:l:ss]) = 0; FDMat([2*l:l:ss, (2*l)+1:l:ss],:) = 0;

  end
