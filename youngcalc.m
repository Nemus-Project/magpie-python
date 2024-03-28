function E  = youngcalc(rho,ldim,h,BCs,ExpFreq,Ntrain)
% YOUNGCALC: what does this do ?
%      E  = YOUNGCALC(rho,ldim,h,BCs,ExpFreq,Ntrain)
%
% this is a function that returns Young's modulus (E) of an experimental plate
% starting from a batch of experimentally measured frequencies, leveraging MAGPIE
%
% Input parameters
%           rho     : the experimental plate density
%           ldim    : a 3X1 array containing the Lx Ly Lz dimensions
%           h       : the grid spacing of the FD scheme
%           BCs     : a 4X2 array containing the rigidities of the boundary supports of the experimental plate
%           ExpFreq : an array contaning the first Nmodes measured modal frequencies
%           Ntrain  : an integer. Must be <= Nmodes. It is the number of training modes out of the available batch
%
% example usage
%
%           ExpFreq = [73.2; 148; 376; 431; 559; 910] ;  %-- these are measured from a plate 
%           rho     = 8765 ;            %-- density [kg/m^3]
%           Lx      = 0.1 ;
%           Ly      = 0.08 ;
%           Lz      = 0.00081 ;
%
% cantilever BCs. Clamped edge along x:
%           BCs = [0,    0;
%                  1e15, 1e15;
%                  0,    0;
%                  0,    0];
%
%           ldim    = [Lx Ly Lz] ;
%           h       = sqrt(Lx*Ly)*0.01 ;  %-- grid spacing [m]
%
%           E  = youngcalc(rho,ldim,h,BCs,ExpFreq,3) ;
%
%--------------------------------------------------------------------------
%

Nmodes = length(ExpFreq) ;

if Ntrain > Nmodes
    disp('Choose Ntrain < total experimental Freqs')
    return
end

TrainFreq = ExpFreq(1:Ntrain) ;
if Ntrain == Nmodes
    TestFreq  = [] ;
else
    TestFreq  = ExpFreq(Ntrain+1:end) ;
end

%-- zero parameters
E0      = 2e11 ;                 %-- Young's modulus [Pa] (just a number here, results shouldnt change if this changes)
nu      = 0.3 ;                  %-- poisson's ratio (average value for metals)


Lx      = ldim(1) ;
Ly      = ldim(2) ;
Lz      = ldim(3) ;
A       = Lx*Ly ;                %-- area

D0      = E0*Lz^3/12/(1-nu^2) ;  %-- zero-rigidity



%-- run magpie and get non-dimensional freqs
Om        = magpie(rho,E0,nu,ldim,h,BCs,Ntrain,"none") ;
OmNDim    = Om ./ sqrt(D0) * sqrt(rho * A^2 * Lz) ;
OmNDimsq  = OmNDim.^2 ;

%-- least-square (LS) optimisation
psi       = (TrainFreq*2*pi).^2 * rho * Lz * A^2 ;
DLS       = (OmNDimsq.'*psi) / (OmNDimsq.' * OmNDimsq) ;
ELS       = DLS / (Lz^3/12/(1-nu^2)) ;
% OmSqLS    = OmNDimsq*DLS ;

%-- launch a numerical simulation to get the frequencies of the numerical model
%-- using the estimated Youngs Mod
%-- and compare against the experimental freqs
NumOm   = magpie(rho,ELS,nu,ldim,h,BCs,Nmodes,"none") ;
NumFreq = NumOm/2/pi ;



if Ntrain < Nmodes

    subplot(2,2,1)
    Y = [TrainFreq,NumFreq(1:Ntrain)] ;
    X = 1:Ntrain ;
    bar(X,Y)
    xlabel('Mode Number') ;
    ylabel('f (Hz)')
    legend('Exp','Num')
    title('Training Set')

    subplot(2,2,3)
    errTrain = (1-Y(:,2)./Y(:,1))*100 ;
    bar(X,errTrain)
    xlabel('Mode Number') ;
    ylabel('rel err (%)')

    subplot(2,2,2)
    Y = [TestFreq,NumFreq(Ntrain+1:end)] ;
    X = Ntrain+1:Nmodes ;
    bar(X,Y)
    xlabel('Mode Number') ;
    ylabel('f (Hz)')
    legend('Exp','Num')
    title('Testing Set')

    subplot(2,2,4)
    errTest = (1-Y(:,2)./Y(:,1))*100 ;
    bar(X,errTest)
    xlabel('Mode Number') ;
    ylabel('rel err (%)')

else

    subplot(2,1,1)
    Y = [TrainFreq,NumFreq(1:Ntrain)] ;
    X = 1:Ntrain ;
    bar(X,Y)
    xlabel('Mode Number') ;
    ylabel('f (Hz)')
    legend('Exp','Num')
    title('Training Set')

    subplot(2,1,2)
    errTrain = (1-Y(:,2)./Y(:,1))*100 ;
    bar(X,errTrain)
    xlabel('Mode Number') ;
    ylabel('rel err (%)')

end



E = ELS ;
