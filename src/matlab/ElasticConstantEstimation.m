clear all
close all
clc

%------------- Training Phase: Obtain the linear modal paramters a b c
%---- global plate parameters
rho      = 439 ;
Lx       = 0.215 ;
Ly       = 0.108 ;
Lz       = 0.001 ;

Ex0      = 10.7e9 ;
Ey0      = 716e6 ;
Gxy0     = 500e6 ;
nux0     = 0.51 ;

scaleVec = [0.8,0.9,1,1.1,1.2] ;
scaleMat = zeros(2,5^2) ;
ind = 0 ;
for p = 1 : 5
    for q = 1 : 5
        ind = ind + 1 ;
        scaleMat(:,ind) = [scaleVec(p);scaleVec(q)] ;
    end
end

load('Eigenfreqs_training_Spruce_S3_1mm_CFFF.csv') ;
OmMat = Eigenfreqs_training_Spruce_S3_1mm_CFFF(:,4) ;
OmMat = reshape(OmMat,[6,25]) ;

temp              = OmMat(4,[3,4,5,9,10,15]) ;
OmMat(4,[3,4,5,9,10,15]) = OmMat(5,[3,4,5,9,10,15]) ;
OmMat(5,[3,4,5,9,10,15]) = temp ;
% return
OmMat(:,21) = [] ;
scaleMat(:,21) = [] ;

%-- init
DxVec = zeros(24,1) ;
DyVec = zeros(24,1) ;
DsVec = zeros(24,1) ;
om1sq = zeros(24,1) ;
om2sq = zeros(24,1) ;
om3sq = zeros(24,1) ;
om4sq = zeros(24,1) ;
om5sq = zeros(24,1) ;
om6sq = zeros(24,1) ;
rxVec = zeros(24,1) ;
ryVec = zeros(24,1) ;

for n = 1 : 24


    Ex      = Ex0*scaleMat(1,n) ;
    Ey      = Ey0*scaleMat(2,n) ;
    Gxy     = Gxy0 ;

    nux     = nux0 ;
    nuy     = Ey/Ex*nux ;

    Dx            = Ex*Lz^3/12/(1-nux*nuy) ;
    Dy            = Ey*Lz^3/12/(1-nux*nuy) ;
    Ds            = Gxy*Lz^3/3 ;

    DxVec(n)      = Dx ;
    DyVec(n)      = Dy ;
    DsVec(n)      = Ds ;
    rxVec(n)      = Dx/Ds ;
    ryVec(n)      = Dy/Ds ;
    om1sq(n)      = OmMat(1,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om2sq(n)      = OmMat(2,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om3sq(n)      = OmMat(3,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om4sq(n)      = OmMat(4,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om5sq(n)      = OmMat(5,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om6sq(n)      = OmMat(6,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;

end

omsq = [om1sq,om2sq,om3sq,om4sq,om5sq,om6sq] ;


FitMat = zeros(4,6) ; %-- this collects the best fit paramters for the six modes
%-- first row: a ; second row: b ; third row: c.

for n = 1 : 6

    linfit = fit([rxVec,ryVec],omsq(:,n),'poly11','Lower', [0 0 0]) ;
    % linfit = fit([rxVec,ryVec],omsq(:,n),'poly11') ;

    aCur = linfit.p10 ;
    bCur = linfit.p01 ;
    cCur = linfit.p00 ;


    FitMat(1,n) = aCur ;
    FitMat(2,n) = bCur ;
    FitMat(3,n) = cCur ;

    %-- compute Rsquared
    SStot = 0 ; SSres = 0 ; meany = mean(omsq(:,n)) ;

    for nFit = 1 : 23
        rx        = rxVec(nFit) ;
        ry        = ryVec(nFit) ;
        yfit      = aCur*rx + bCur*ry + cCur ;
        y         = omsq(nFit,n) ;
        SStot     = SStot + (y-meany).^2;                    % Total Sum-Of-Squares
        SSres     = SSres + (y-yfit).^2 ;                    % Residual Sum-Of-Squares

    end
    Rsq       = 1-SSres/SStot  ;
    FitMat(4,n) = Rsq ;

end



%------------- End Training Phase ........................................
%-------------------------------------------------------------------------

%--------------- Plate: Measured spruce plate

indCell = {
    [1 2 3]
    [1 2 4]
    [1 2 5 ]
    [1 2 6 ]
    [1 3 4 ]
    [1 3 5 ]
    [1 3 6 ]
    [1 4 5 ]
    [1 4 6 ]
    [1 5 6 ]
    [2 3 4]
    [2 3 5]
    [2 3 6]
    [2 4 5]
    [2 4 6]
    [2 5 6]
    [3 4 5]
    [3 4 6]
    [3 5 6]
    [4 5 6]
    [1 2 3 4]
    [1 2 3 5]
    [1 2 3 6]
    [1 2 4 5]
    [1 2 4 6]
    [1 2 5 6]
    [1 3 4 5]
    [1 3 4 6]
    [1 3 5 6]
    [1 4 5 6]
    [2 3 4 5]
    [2 3 4 6]
    [2 3 5 6]
    [2 4 5 6]
    [3 4 5 6]
    [1 2 3 4 5]
    [1 2 3 4 6]
    [1 2 3 5 6]
    [1 2 4 5 6]
    [1 3 4 5 6]
    [2 3 4 5 6]
    [1 2 3 4 5 6]
    } ;
Ntot = 42;

errDx = zeros(Ntot,1) ;
errDy = zeros(Ntot,1) ;
errDs = zeros(Ntot,1) ;
DxLS  = zeros(Ntot,1) ;
DyLS  = zeros(Ntot,1) ;
DsLS  = zeros(Ntot,1) ;
condNumb = zeros(Ntot,1) ;

rho      = 439;
Lx       = 0.215 ;
Ly       = 0.108 ;
Lz       = 0.00461 ;
% Lz       = 0.001;

% Ex0      = 10.7e9 ;
% Ey0      = 716e6 ;
% Gxy0     = 500e6 ;
% nux0     = 0.51 ;

nux      = nux0 ;
nuy      = Ey0/Ex0*nux ;

% 



 ExpFreqs = [82.8
 160.7
 517.2
 635.1
 657.3
 1059.5] ;

% % from same measurements repeated on 29 01 24
%   ExpFreqs = [80.8
%  156.9
%  507.2
%  629.4
%  670.9
%  1069.9] ;

 %  % from same measurements repeated on 29 01 24_2
 %  ExpFreqs = [82.4
 % 157.2
 % 510.25
 % 631
 % 669.6
 % 1070.1] ;


Om0=2*pi*ExpFreqs;

for nCase = 1: Ntot


    Nselect = indCell(nCase,:) ;
    Nselect = cell2mat(Nselect) ;


    a        = (FitMat(1,Nselect)).' 
    b        = (FitMat(2,Nselect)).' ;
    c        = (FitMat(3,Nselect)).' ;
    rsquared = (FitMat(4,Nselect)).' ;



    Om = Om0(Nselect) ;


    psi  = Om.^2 ;


    eta = (rho*Lz*Lx^2*Ly^2)^(-1) 

    X = [eta*a eta*b eta*c] 
     LSmat = (X.' * X) \ (X.')
     condNumb(nCase) = cond(LSmat) 


    temp = (X.' * X) \ (X.'  * psi) ;
% return

    DxLS(nCase)     = temp(1) ;
    DyLS(nCase)     = temp(2) ;
    DsLS(nCase)     = temp(3) ;


end
% return

disp('Test Experimental --------------------------------')


%Obtained Ex, Ey, Gxy



digits(3) ;
disp('Ex Ey Gxy')

ExLSall           = (DxLS*(12*(1-nux*nuy)))/Lz^3 ;
EyLSall           = (DyLS*(12*(1-nux*nuy)))/Lz^3 ;
GxyLSall          = 3/Lz^3*(DsLS) ;

ExLS              = [(1:Ntot)',ExLSall] ;
EyLS              = [(1:Ntot)',EyLSall] ;
GxyLS             = [(1:Ntot)',GxyLSall] ;

EMat   = [ExLSall EyLSall GxyLSall] ;
vpa(EMat) 


ExLSNeg              = find(ExLS(:,2)<0) ;   ExLS(ExLSNeg,:) = [] ;
EyLSNeg              = find(EyLS(:,2)<0) ;   EyLS(EyLSNeg,:) = [] ;
GxyLSNeg             = find(GxyLS(:,2)<0) ;  GxyLS(GxyLSNeg,:) = [] ;


gigiEx   = isoutlier(ExLS(:,2),"percentiles",[20 80]) ;
%gigiEx   = isoutlier(ExLS(:,2),"quartiles") ;
gigiEy   = isoutlier(EyLS(:,2),"percentiles",[20 80]) ;
%gigiEy   = isoutlier(EyLS(:,2),"quartiles") ;
gigiGxy  = isoutlier(GxyLS(:,2),"percentiles",[20 80]) ;
%gigiGxy  = isoutlier(GxyLS(:,2),"quartiles") ;

rEx      = find(gigiEx) ;
indOutEx = ExLS(rEx,1) ;
ExLS(rEx,:)  = [] ;

rEy      = find(gigiEy) ;
indOutEy = EyLS(rEy,1) ;
EyLS(rEy,:)  = [] ;

rGxy      = find(gigiGxy) ;
indOutGxy = GxyLS(rGxy,1) ;
GxyLS(rGxy,:)  = [] ;

meanEx = mean(ExLS(:,2)) ; stdEx = std(ExLS(:,2)) ;
meanEy = mean(EyLS(:,2)) ; stdEy = std(EyLS(:,2)) ;
meanGxy = mean(GxyLS(:,2)) ; stdGxy = std(GxyLS(:,2)) ;


meanVec = vpa([meanEx meanEy meanGxy])
stdVec  = vpa([stdEx/meanEx*100 stdEy/meanEy*100 stdGxy/meanGxy*100])


Nstd = 4 ;

subplot(1,3,1)
plot((0:Ntot+1),ones(Ntot+2,1)*meanEx/1e9,'--k') ; xlim([0,Ntot+1]); ylim([meanEx-Nstd*stdEx,meanEx+Nstd*stdEx]/1e9);
x = [0 Ntot+1 Ntot+1 0];
y = [meanEx-stdEx meanEx-stdEx meanEx+stdEx meanEx+stdEx]/1e9;

patch(x,y,'red','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(ExLS(:,1),ExLS(:,2)/1e9,'linestyle','none','marker','o','color','k'); 
% return
%plot(ExLSNeg,ExLSall(ExLSNeg),'linestyle','none','marker','*','color','k') ; 
plot(indOutEx,ExLSall(indOutEx)/1e9,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$e_x$ (GPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;

subplot(1,3,2)
plot((0:Ntot+1),ones(Ntot+2,1)*meanEy/1e6,'--k') ; xlim([0,Ntot+1]); ylim([meanEy-Nstd*stdEy,meanEy+Nstd*stdEy]/1e6)
x = [0 Ntot+1 Ntot+1 0];
y = [meanEy-stdEy meanEy-stdEy meanEy+stdEy meanEy+stdEy]/1e6;
patch(x,y,'blue','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(EyLS(:,1),EyLS(:,2)/1e6,'linestyle','none','marker','o','color','k'); 
%plot(EyLSNeg,EyLSall(EyLSNeg),'linestyle','none','marker','*','color','k') ; 
plot(indOutEy,EyLSall(indOutEy)/1e6,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$e_y$ (MPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;


subplot(1,3,3)
plot((0:Ntot+1),ones(Ntot+2,1)*meanGxy/1e6,'--k') ; xlim([0,Ntot+1]); ylim([meanGxy-Nstd*stdGxy,meanGxy+Nstd*stdGxy]/1e6)
x = [0 Ntot+1 Ntot+1 0];
y = [meanGxy-stdGxy meanGxy-stdGxy meanGxy+stdGxy meanGxy+stdGxy]/1e6;
patch(x,y,'green','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(GxyLS(:,1),GxyLS(:,2)/1e6,'linestyle','none','marker','o','color','k'); 
%plot(GxyLSNeg,GxyLSall(GxyLSNeg),'linestyle','none','marker','*','color','k') ; 
plot(indOutGxy,GxyLSall(indOutGxy)/1e6,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$g_{xy}$ (MPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;

sgtitle('Tonewood 3','interpreter','latex') 

return

%[1.38e+10, 8.31e+8, 6.05e+8]
NumFreqs = [51.240
99.849
315.05
335.54
394.33
629.48] ;



err     = NumFreqs - ExpFreqs;
errPct  = 100 * (NumFreqs - ExpFreqs)./ExpFreqs ;
errCent = 1200*log2(NumFreqs./ExpFreqs) ;

digits(3)
[vpa(ExpFreqs) vpa(NumFreqs) vpa(err) vpa(errPct) vpa(errCent)]

