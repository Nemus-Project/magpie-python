%---- modal synthesis of the experimentally tuned plate

clear all
close all
clc

%--------------------------------------------------------------------------
% CUSTOM PARAMETERS (can change below)

%--- Plate Parameters
% physical and elastic parameters
Lx      = 0.1 ;             %-- length along x [m]
Ly      = 0.08 ;            %-- length along y [m]
Lz      = 0.81e-3 ;         %-- thickness [m]
E       = 101e9 ;           %-- Young's mod [Pa]
rho     = 8765 ;            %-- density [kg/m^3]
nu      = 0.3 ;             %-- poisson's ratio

for nPl = 1 : 3


    for nAtt = 1 : 3

        if nPl == 1
            % input and output locs
            in       = [0.052 ,0.040] ;   % values in m
            out      = [0.072,0.042] ;    % values in m
            if nAtt == 1
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [72.1 150 382 460 578 890 920].' * 2 * pi ;
            elseif nAtt == 2
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [73 147 383 430 560 905 920].' * 2 * pi ;
            elseif nAtt == 3
                T60s     = [0.8, 0.4, 0.6, 0.3, 0.4, 0.5, 0.5].' ;
                ExpFreqs = [73 147 383 430 560 905 920].' * 2 * pi ;
            end
        elseif nPl == 2
            % input and output locs
            in      = [0.045, 0.035] ;   % values in m
            out     = [0.072,0.042] ;    % values in m
            if nAtt == 1
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [72.1 150 382 460 578 890 920].' * 2 * pi ;
            elseif nAtt == 2
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [73 147 384 430 560 905 920].' * 2 * pi ;
            elseif nAtt == 3
                T60s     = [0.8, 1.2, 0.6, 0.3, 0.4, 0.5, 0.5].' ;
                ExpFreqs = [73 147 384 430 560 905 920].' * 2 * pi ;
            end

        elseif nPl == 3
            % input and output locs
            in      = [0.035, 0.064] ;   % values in m
            out     = [0.072,0.042] ;    % values in m
            if nAtt == 1
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [72.1 150 382 460 578 890 920].' * 2 * pi ;
            elseif nAtt == 2
                T60s     = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].' ;
                ExpFreqs = [73 147 383 430 560 905 920].' * 2 * pi ;
            elseif nAtt == 3
                T60s     = [0.8, 1.2, 0.6, 0.3, 0.4, 0.5, 0.5].' ;
                ExpFreqs = [73 147 383 430 560 905 920].' * 2 * pi ;
            end
        end

        %-- damp and stiff matrix
        C        = 6*log(10)./T60s ;
        Omsq     = ExpFreqs.^2 ;

        % elastic constants around the edges (this allows to set the various bcs)
        K0y     = 0 ;
        R0y     = 0 ;
        Kx0     = 1e15 ;
        Rx0     = 1e15 ;
        KLy     = 0 ;
        RLy     = 0 ;
        KxL     = 0 ;
        RxL     = 0 ;

        % mesh size REALLY IMPORTANT
        h       =  sqrt(Lx*Ly)*0.01 ;

        % total number of modes
        Nmodes   = 7 ;

        % bottom and top freqs for plots and normalisation
        flow     = 20 ;
        fhigh    = 1000 ;
        %----------------------------------

        %--- Simulation Parameters
        T        = 1 ;              % sim length [s]
        fs       = 44100 ;          % sample rate [Hz]
        Fin      = 1 ;              % force amplitude [N]
        twid     = 0.003 ;          % forcing temporal width [s]
        %----------------------------------

        %--------------------------------------------------------------------------



        %--------------------------------------------------------------------------
        % DERIVED PARAMETERS (DON'T TOUCH BELOW)

        Ts           = round(T*fs) ;
        k            = 1 / fs ;
        tvec         = (0:Ts-1)*k ;
        fvec         = (0:Ts-1)*fs/Ts ;
        [~,indLow]   = sort(abs(fvec-flow)) ;
        [~,indHigh]  = sort(abs(fvec-fhigh)) ;
        fvec         = fvec(indLow:indHigh) ;

        % eigenvalue problem
        [Q,Om,Nx,Ny,biHarm] = biharmonicEigs(rho,E,nu,Lx,Ly,Lz,h,K0y,R0y,Kx0,Rx0,KLy,RLy,KxL,RxL,Nmodes) ;
        Q6 = Q(:,6) ;
        Q7 = Q(:,7) ;
        Q(:,6) = Q7 ;
        Q(:,7) = Q6 ;
        %Omsq = Om.^2 ;

        % build input and output

        Jin      = zeros((Nx+1)*(Ny+1),1) ;

        tempinx  = in(1)/Lx*Nx ;
        Min      = floor(tempinx) ;
        alx      = tempinx - Min ;

        tempiny  = in(2)/Ly*Ny ;
        Nin      = floor(tempiny) ;
        aly      = tempiny - Nin ;

        Jin((Ny+1)*Min+Nin+1)     = alx*aly ;
        Jin((Ny+1)*(Min+1)+Nin+1) = (1-alx)*aly ;
        Jin((Ny+1)*Min+Nin+2)      = alx*(1-aly) ;
        Jin((Ny+1)*(Min+1)+Nin+2) = (1-alx)*(1-aly) ;


        Jin = Jin / h^2 / rho / Lz ;
        Jin = pinv(Q)*Jin ;


        Jout      = zeros(1,(Nx+1)*(Ny+1)) ;

        tempinx  = out(1)/Lx*Nx ;
        Min      = floor(tempinx) ;
        alx      = tempinx - Min ;

        tempiny  = out(2)/Ly*Ny ;
        Nin      = floor(tempiny) ;
        aly      = tempiny - Nin ;

        Jout((Ny+1)*Min+Nin+1)      = alx*aly ;
        Jout((Ny+1)*(Min+1)+Nin+1)  = (1-alx)*aly ;
        Jout((Ny+1)*Min+Nin+2)      = alx*(1-aly) ;
        Jout((Ny+1)*(Min+1)+Nin+2)  = (1-alx)*(1-aly) ;

        %-- input forcing

        Nfin  = floor(twid*fs) ;
        fin   = zeros(Ts,1) ;
        fin(1:Nfin) = 1 - cos(2*pi*(0:Nfin-1)/Nfin) ;


        %--- init
        vm = zeros(Nmodes,1) ;
        v0 = zeros(Nmodes,1) ;
        out = zeros(Ts,1) ;

        for n = 1 : Ts

            fin = 0 ;
            if n == 10
                fin = 1 ;
            end
            vp     = (2*v0-vm - k^2*Omsq.*v0 + k*C.*vm + k^2*Jin*fin)./(1+k*C) ;
            out(n) = Jout * (Q*v0) ;
            vm     = v0 ;
            v0     = vp ;

        end

        acc             = [0;diff(out,2);0] ;
        acc             = acc / max(abs(acc)) ;
        fftacc          = abs(fft(acc.*hamming(Ts))) ;
        fftacc          = fftacc(indLow:indHigh) ;
        fftacc          = fftacc/max(fftacc) ;

        % -------------------------------------------------------------------------
        %-- load and postprocess experimental signals
        %-- DONT TOUCH BELOW

        %-- load signals
        [H1ACC1,fs]    = audioread('audio/300623_COPPER PLATE_H1ACC1.wav') ;
        H1ACC1         = H1ACC1(:,1) ;
        TsExp          = length(H1ACC1) ;
        hann2          = hamming(TsExp) ;
        H1ACC1         = H1ACC1.* hamming(TsExp) ;

        [H2ACC1,~]     = audioread('audio/300623_COPPER PLATE_H2ACC1.wav') ;
        H2ACC1         = H2ACC1(:,1) ;
        H2ACC1         = H2ACC1.* hamming(TsExp) ;

        [H3ACC1,~]     = audioread('audio/300623_COPPER PLATE_H3ACC1.wav') ;
        H3ACC1         = H3ACC1(:,1) ;
        H3ACC1         = H3ACC1.* hamming(TsExp) ;


        %-- compute spectra and smooth
        spc1         = abs(fft(H1ACC1)) ;
        spc1         = fastsmooth(spc1,14,3,1) ;

        spc2         = abs(fft(H2ACC1)) ;
        spc2         = fastsmooth(spc2,14,3,1) ;

        spc3         = abs(fft(H3ACC1)) ;
        spc3         = fastsmooth(spc3,14,3,1) ;


        fvExp        = (0:TsExp-1)*fs/TsExp ;               %-- frequency vector
        tvExp        = (0:TsExp-1)/fs ;                   %-- frequency vector

        [~,indLow]   = sort(abs(fvExp-flow)) ;
        [~,indHigh]  = sort(abs(fvExp-fhigh)) ;
        fvExp        = fvExp(indLow:indHigh) ;
        spc1         = spc1(indLow:indHigh) ;
        spc2         = spc2(indLow:indHigh) ;
        spc3         = spc3(indLow:indHigh) ;
        spc1         = spc1/max(spc1) ;
        spc2         = spc2/max(spc2) ;
        spc3         = spc3/max(spc3) ;

        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------

        % PLOT STUFF
        subplot(3,1,nPl)
        plot(fvec,20*log10(fftacc),'linewidth',1.0); xlim([0,910]) ; ylim([-70,10]);
        hold on ; drawnow ;
        if nAtt == 3
            if nPl == 1
                plot(fvExp,20*log10(spc1),'k','linewidth',1.0);
            elseif nPl == 2
                plot(fvExp,20*log10(spc2),'k','linewidth',1.0);
            elseif nPl == 3
                plot(fvExp,20*log10(spc3),'k','linewidth',1.0);
            end
        end

        xlabel('f (Hz)','interpreter','latex') ; ylabel('$\hat u$ (dB)','interpreter','latex') ;
        set(gca, 'TickLabelInterpreter', 'latex','fontsize',13);

    end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% SOME DATA FROM THE EXCEL ....
ResonantFreqsHz = [73.0, 430.0, 518.0, 617.0, 909.0, 933.0]' ;
vpa(ResonantFreqsHz)


T60sSecs = [0.579, 0.143, 0.0468, 0.0467, 0.263, 0.0455] ;


ResonantFreqsHz = [73.3, 147.0, 424.0, 560.0, 905.0, 922.0]' ;
vpa(ResonantFreqsHz)


T60sSecs = [0.951, 1.48, 0.0971, 0.159, 0.163, 0.183] ;


ResonantFreqsHz = [73.2, 148.0, 433.0, 564.0, 907.0, 924.0]' ;
vpa(ResonantFreqsHz)

T60sSecs = [1.01, 1.5, 0.165, 0.214, 0.109, 0.548] ;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

