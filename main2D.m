%%
%   Efficient Point-Mass Filter for 2D.
%   author: pesslovany@gmail.com
%   Papers reference will be added after publication

%% Parameters and system simulation
clc
clear
close all

format shortG

load("mapTAN.mat")

% System parameters
nx = 2; % state dimension
nz = 1; % measurement dimension
kf = 15; % final time

Q = diag([100 100]); % system noise covariance
invQ = inv(Q);
R = 2^2; % measurement noise covariance for both modes
invR = inv(R);

% PMF parameters
sFactor = 4; % scaling factor (number of sigmas covered by the grid)

Ts = 1; % time step
turn = deg2rad(15); % turn rate

% Dynamics
F = eye(2); % Discrete dynamics

A = logm(F); % Continuous dynamics (zeros)
Qc = Q; %Continuous Q

u = 50*ones(2,1); % known input

ffunct = @(x,w,k) F*x + u + w; % state equation for trajectory generation
ffunctSTD = @(x,w,k) F*x + w; % state equation for filters (u incorporated done by moving grid)

hfunct = @(x,v,k) interp2(map_m.x,map_m.y,map_m.z,x(1,:),x(2,:)) + v; % measurement equation

%% Gauss mix noise v
meanV = [0 20]; % Mean values of components of measurement noise
wV = [0.5 0.5]; % Weights of measurement noise

V.pdf = gmdistribution(meanV',R,wV); % Measurement PDF

MC = 1;
for mc = 1:1:MC
    mc
    Npa = 41;
    N = Npa^nx;
    % Initial condition - Gaussian
    meanX0 = [mean(map_m.x,"all"); mean(map_m.y,"all")]; % initial cond mean value
    varX0 = [160 20;
        20 90]; %initial cond covariance

    % PF parameters
    noPart = 1681; % number of particles fof PF
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF


    x = zeros(nx,kf);
    x(:,1) = mvnrnd(meanX0,varX0); % initial state

    % System simulation
    w = mvnrnd(zeros(nx,1),Q,kf-1)'; % system noise

    v = random(V.pdf,kf)'; % measurement noise

    for i=1:kf-1
        % trajectory generation
        x(:,i+1) = ffunct(x(:,i),w(:,i),i);
    end
    % measurement generation
    z = hfunct(x,v,1:length(x));


    %% PMF init

    % Initial grid
    [predGrid, GridDelta, gridDimOld] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator);% Adding probabilities to points

    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); %Denominator for convolution in predictive step


    for k = 1:1:kf
        %% STD
        tic
        %% Mesurement update
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct);

        % Measurement mean and var
        measMean1(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean1(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar1(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k)); % Measurement update variance

        % Expected pred mean/var using KF
        predMeanEst = F*measMean1(:,k);
        predVarEst = F*measVar1(:,:,k)*F' + Q;

        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; % Back to filtering space
        maxF = max(gridBoundWantCorners,[],2); % Min/Max meas corners
        minF = min(gridBoundWantCorners,[],2);
        for ind3 = 1:1:nx % Creation of filtering grid so that it creates wanted predictive grid.
            gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);
            gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));
        end
        measGridNew = combvec(gridDim);
        GridDelta(:,k) = gridStep; % Grid step size


        % INSTEAD of interpolating(finding) the new rectangular grid measurement grid values on the
        % scattered predictive grid from last step, we are interpolating the measurement grid values transformed
        % inversely inv(F)*measGridNew on the old predictive grid transformed also inversly inv(F)*predGrid (which is
        % actually measurement grid  from last step) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew;
        else
            filtGridInterpInvTrsf = inv(F)*(measGridNew);
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:))';

        gridDimOld{1} = gridDim{1}+u(1);
        gridDimOld{2} = gridDim{2}+u(2);
        %------------- Time update-------------------

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size

        % STD time-update
        predDensityProb2cub =  pmfUpdateSTD(measGridNew,measPdf,GridDelta,ffunctSTD,predGrid,nx,k,invQ,predDenDenomW,N);

        predDensityProb = reshape(predDensityProb2cub,N,1); % back to computational space
        predDensityProb = predDensityProb./(sum(predDensityProb)*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predGrid = predGrid+u; % Input

        tocPMF1(k) = toc;

    end


    %% PMF init

    % Initial grid
    [predGrid, GridDelta, gridDimOld] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial weights
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator);% Adding probabilities to points

    % Auxiliary variables
    halfGrid = ceil(N/2);
    for k = 1:1:kf
        %% Efficient FFT
        tic
        %----------------------- Measurement update---------------------
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct);

        % Measurement mean and var
        measMean2(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean2(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar2(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k));

        % Expected pred mean/var
        predMeanEst = F*measMean2(:,k);
        predVarEst = F*measVar2(:,:,k)*F' + Q;

        % Setup the meas grid on which to interpolate on the basis of wanted
        % pred grid sigma
        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; % Back to filtering space
        maxF = max(gridBoundWantCorners,[],2); % Min/Max meas corners
        minF = min(gridBoundWantCorners,[],2);
        for ind3 = 1:1:nx % Creation of filtering grid so that it creates wanted predictive grid.
            gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);
            gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));
        end
        measGridNew = combvec(gridDim);
        GridDelta(:,k) = gridStep; % Grid step size

        % INSTEAD of interpolating(finding) the new rectangular grid measurement grid values on the
        % scattered predictive grid from last step, we are interpolating the measurement grid values transformed
        % inversely inv(F)*measGridNew on the old predictive grid transformed also inversly inv(F)*predGrid (which is
        % actually measurement grid  from last step) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew;
        else
            filtGridInterpInvTrsf = inv(F)*(measGridNew);
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:))';

        gridDimOld{1} = gridDim{1}+u(1); % Input
        gridDimOld{2} = gridDim{2}+u(2);

        %------------- Time update-------------------

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size


        % ULTRA FAST PMF
        filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
        filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa,Npa); % Into physical space

        pom = (predGrid(:,halfGrid)'-(predGrid)');
        TPMrow = ((exp(sum(-0.5*pom*invQ.*pom,2)))/predDenDenomW)';% Middle row of the TPM matrix
        TPMrowCubPom = reshape(TPMrow,Npa,Npa); % Into physical space

        lfftfun = @(l) 2^nextpow2(l); % Next power of two for FFT
        dims = 1:1:nx;

        for dim=dims % FFT over all dimensions
            % compute the FFT length
            l = lfftfun(Npa+Npa-1);
            TPMrowCubPom = fft(TPMrowCubPom,l,dim); % FFT
        end

        predDensityProb2cub = convnfft(filtDenDOTprodDeltasCub, TPMrowCubPom, Npa); % FFT Convolution

        predDensityProb = reshape(predDensityProb2cub,N,1); % back to computational space
        predDensityProb = predDensityProb./(sum(predDensityProb)*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predGrid = predGrid+u; % Input

        tocPMF2(k) = toc;
    end

    %% PMF init

    % Initial grid
    [predGrid, GridDelta, gridDimOld] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial weights
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator);% Adding probabilities to points

    % Auxiliary variables
    for k = 1:1:kf
        %% Efficient SFT 0.01
        tic
        %----------------------- Measurement update---------------------
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct);

        % Measurement mean and var
        measMean3(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean3(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar3(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k));

        % Expected pred mean/var
        predMeanEst = F*measMean3(:,k);
        predVarEst = F*measVar3(:,:,k)*F' + Q;

        % Setup the meas grid on which to interpolate on the basis of wanted
        % pred grid sigma
        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; % Back to filtering space
        maxF = max(gridBoundWantCorners,[],2); % Min/Max meas corners
        minF = min(gridBoundWantCorners,[],2);
        for ind3 = 1:1:nx % Creation of filtering grid so that it creates wanted predictive grid.
            gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);
            gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));
        end
        measGridNew = combvec(gridDim);
        GridDelta(:,k) = gridStep; % Grid step size

        % INSTEAD of interpolating(finding) the new rectangular grid measurement grid values on the
        % scattered predictive grid from last step, we are interpolating the measurement grid values transformed
        % inversely inv(F)*measGridNew on the old predictive grid transformed also inversly inv(F)*predGrid (which is
        % actually measurement grid  from last step) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew;
        else
            filtGridInterpInvTrsf = inv(F)*(measGridNew);
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:))';

        gridDimOld{1} = gridDim{1}+u(1); % Input
        gridDimOld{2} = gridDim{2}+u(2);

        %------------- Time update-------------------

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size


        % ULTRA FAST PMF
        filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
        filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa,Npa); % Into physical space

        % Continuous prediction using fast sine transform
        dt = 0.01; % Numerical method time step
        a = (Qc(1,1)*dt)/(2*GridDelta(1,k)^2); % Finite difference diffusion matrix diagonal values
        b = 1 - (Qc(1,1)*dt)/(2*GridDelta(1,k)^2) - (Qc(2,2)*dt)/(2*GridDelta(2,k)^2) - dt*trace(A);
        c = Qc(2,2)*dt/(2*GridDelta(2,k)^2);

        count = (1:1:Npa); % indexes
        lambdaJ = b + 2*a*cos((count*pi)/(Npa+1)) + 2*c*cos((count'*pi)/(Npa+1)); % Eigenvalues of diffusion matrix

        asd = (dtt2D(filtDenDOTprodDeltasCub,5)'/(Npa+1)); % SFT calculation of prediction
        predDensityProb2cub = dtt2D(lambdaJ.^(1/dt).*asd,5)';

        predDensityProb = reshape(predDensityProb2cub,N,1); % back to computational space
        predDensityProb = predDensityProb./(sum(predDensityProb)*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predGrid = predGrid+u; % Input


        tocPMF3(k) = toc;
    end

    %% PMF init

    % Initial grid
    [predGrid, GridDelta, gridDimOld] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial weights
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator);% Adding probabilities to points

    % Auxiliary variables
    for k = 1:1:kf
        %% Efficient SFT 0.001
        tic
        %----------------------- Measurement update---------------------
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct);

        % Measurement mean and var
        measMean4(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean4(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar4(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k));

        % Expected pred mean/var
        predMeanEst = F*measMean4(:,k);
        predVarEst = F*measVar4(:,:,k)*F' + Q;

        % Setup the meas grid on which to interpolate on the basis of wanted
        % pred grid sigma
        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; % Back to filtering space
        maxF = max(gridBoundWantCorners,[],2); % Min/Max meas corners
        minF = min(gridBoundWantCorners,[],2);
        for ind3 = 1:1:nx % Creation of filtering grid so that it creates wanted predictive grid.
            gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);
            gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));
        end
        measGridNew = combvec(gridDim);
        GridDelta(:,k) = gridStep; % Grid step size

        % INSTEAD of interpolating(finding) the new rectangular grid measurement grid values on the
        % scattered predictive grid from last step, we are interpolating the measurement grid values transformed
        % inversely inv(F)*measGridNew on the old predictive grid transformed also inversly inv(F)*predGrid (which is
        % actually measurement grid  from last step) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew;
        else
            filtGridInterpInvTrsf = inv(F)*(measGridNew);
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:))';

        gridDimOld{1} = gridDim{1}+u(1);
        gridDimOld{2} = gridDim{2}+u(2);

        %------------- Time update-------------------

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size


        % ULTRA FAST PMF
        filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
        filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa,Npa); % Into physical space

        dt = 0.001; % Numerical method time step
        a = (Qc(1,1)*dt)/(2*GridDelta(1,k)^2); % Finite difference diffusion matrix diagonal values
        b = 1 - (Qc(1,1)*dt)/(2*GridDelta(1,k)^2) - (Qc(2,2)*dt)/(2*GridDelta(2,k)^2) - dt*trace(A);
        c = Qc(2,2)*dt/(2*GridDelta(2,k)^2);

        count = (1:1:Npa); %pom
        lambdaJ = b + 2*a*cos((count*pi)/(Npa+1)) + 2*c*cos((count'*pi)/(Npa+1)); % Eigenvalues of diffusion matrix

        asd = (dtt2D(filtDenDOTprodDeltasCub,5)'/(Npa+1)); % SFT calculation of prediction
        predDensityProb2cub = dtt2D(lambdaJ.^(1/dt).*asd,5)';

        predDensityProb = reshape(predDensityProb2cub,N,1); % back to vector
        predDensityProb = predDensityProb./(sum(predDensityProb)*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predGrid = predGrid+u;


        tocPMF4(k) = toc;
    end


    %% PMF init

    % Initial grid
    [predGrid, GridDelta, gridDimOld] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    measGridOld = predGrid;

    % Initial weights
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator);% Adding probabilities to points

    % Auxiliary variables
    for k = 1:1:kf
        %% Efficient SFT 0.001
        tic
        %----------------------- Measurement update---------------------
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct);

        % Measurement mean and var
        measMean5(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean5(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar5(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k));

        % Expected pred mean/var
        predMeanEst = F*measMean5(:,k);
        predVarEst = F*measVar5(:,:,k)*F' + Q;

        % Setup the meas grid on which to interpolate on the basis of wanted
        % pred grid sigma
        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; % Back to filtering space
        maxF = max(gridBoundWantCorners,[],2); % Min/Max meas corners
        minF = min(gridBoundWantCorners,[],2);
        for ind3 = 1:1:nx % Creation of filtering grid so that it creates wanted predictive grid.
            gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);
            gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));
        end
        measGridNew = combvec(gridDim);
        GridDelta(:,k) = gridStep; % Grid step size

        % INSTEAD of interpolating(finding) the new rectangular grid measurement grid values on the
        % scattered predictive grid from last step, we are interpolating the measurement grid values transformed
        % inversely inv(F)*measGridNew on the old predictive grid transformed also inversly inv(F)*predGrid (which is
        % actually measurement grid  from last step) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew;
        else
            filtGridInterpInvTrsf = inv(F)*(measGridNew); %#ok<*MINV> 
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:))';

        gridDimOld{1} = gridDim{1}+u(1);
        gridDimOld{2} = gridDim{2}+u(2);

        %------------- Time update-------------------

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        measGridOld = measGridNew; % Measurement grid for next for cycle
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size


        % ULTRA FAST PMF
        filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
        filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa,Npa); % Into physical space

        dt = 0.005; % Numerical method time step
        a = (Qc(1,1)*dt)/(2*GridDelta(1,k)^2); % Finite difference diffusion matrix diagonal values
        b = 1 - (Qc(1,1)*dt)/(2*GridDelta(1,k)^2) - (Qc(2,2)*dt)/(2*GridDelta(2,k)^2) - dt*trace(A);
        c = Qc(2,2)*dt/(2*GridDelta(2,k)^2);

        count = (1:1:Npa); %pom
        lambdaJ = b + 2*a*cos((count*pi)/(Npa+1)) + 2*c*cos((count'*pi)/(Npa+1)); % Eigenvalues of diffusion matrix

        asd = (dtt2D(filtDenDOTprodDeltasCub,5)'/(Npa+1)); % SFT calculation of prediction
        predDensityProb2cub = dtt2D(lambdaJ.^(1/dt).*asd,5)';

        predDensityProb = reshape(predDensityProb2cub,N,1); % back to vector
        predDensityProb = predDensityProb./(sum(predDensityProb)*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predGrid = predGrid+u;


        tocPMF5(k) = toc;
    end

    for k = 1:1:kf
        %% Partcle filter - Bootstrap
        tic
        % Measurement step
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); %Prediction density grid through measurement EQ
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom); % Measurement weights
        w = w/sum(w);
        xEst(:,k) = ksiPrior*w; % mean filtering estimate
        varEst(:,k) = ksiPrior.^2*w - xEst(:,k).^2; % diagonal of filtering covariance matrix

        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPart);
        I = binarySearch(cumW, randomN, 'first');
        %I = arrayfun(@(x) find(cumW>x,1,"first"),randomN); %alternative

        % Time-update step
        ksi = ksiPrior(:,I);
        ksiPrior = F*ksi + u + mvnrnd(zeros(1,nx),Q,noPart)';

        tocPF(k) = toc;
    end

    % Evaluation
    rmsePMF(:,mc) = sqrt(mean((x-measMean1).^2,2)); %#ok<*SAGROW> 
    astdPMF11(mc) = mean(sqrt(measVar1(1,1,:)));
    astdPMF22(mc) = mean(sqrt(measVar1(2,2,:)));
    tocPMFavg(:,mc) = mean(tocPMF1);

    rmsePMF2(:,mc) = sqrt(mean((x-measMean2).^2,2));
    astdPMF211(mc) = mean(sqrt(measVar2(1,1,:)));
    astdPMF222(mc) = mean(sqrt(measVar2(2,2,:)));
    tocPMFavg2(:,mc) = mean(tocPMF2);

    rmsePMF3(:,mc) = sqrt(mean((x-measMean3).^2,2));
    astdPMF311(mc) = mean(sqrt(measVar3(1,1,:)));
    astdPMF322(mc) = mean(sqrt(measVar3(2,2,:)));
    tocPMFavg3(:,mc) = mean(tocPMF3);

    rmsePMF4(:,mc) = sqrt(mean((x-measMean4).^2,2));
    astdPMF411(mc) = mean(sqrt(measVar4(1,1,:)));
    astdPMF422(mc) = mean(sqrt(measVar4(2,2,:)));
    tocPMFavg4(:,mc) = mean(tocPMF4);

    rmsePMF5(:,mc) = sqrt(mean((x-measMean5).^2,2));
    astdPMF511(mc) = mean(sqrt(measVar5(1,1,:)));
    astdPMF522(mc) = mean(sqrt(measVar5(2,2,:)));
    tocPMFavg5(:,mc) = mean(tocPMF5);



    rmsePF(:,mc) = sqrt(mean((x-xEst).^2,2));
    astdPF11(mc) = mean(sqrt(varEst(1,:)));
    astdPF22(mc) = mean(sqrt(varEst(2,:)));
    tocPFavg(:,mc) = mean(tocPF);
end

%% Evaluation

rmsePMFout = mean(rmsePMF,2);
tocPMFavgOut = mean(tocPMFavg,2);

rmsePMFout2 = mean(rmsePMF2,2);
tocPMFavgOut2 = mean(tocPMFavg2,2);

rmsePMFout3 = mean(rmsePMF3,2);
tocPMFavgOut3 = mean(tocPMFavg3,2);

rmsePMFout4 = mean(rmsePMF4,2);
tocPMFavgOut4 = mean(tocPMFavg4,2);

rmsePMFout5 = mean(rmsePMF5,2);
tocPMFavgOut5 = mean(tocPMFavg5,2);


rmsePFout = mean(rmsePF,2);
tocPFavgOut = mean(tocPFavg,2);



T2 = table([ rmsePMFout(1) rmsePMFout2(1) rmsePMFout3(1) rmsePMFout5(1) rmsePMFout4(1)  rmsePFout(1)]',...
    [ rmsePMFout(2) rmsePMFout2(2) rmsePMFout3(2) rmsePMFout5(2) rmsePMFout4(2)  rmsePFout(2)]',...
    [ mean(astdPMF11)  mean(astdPMF211)  mean(astdPMF311) mean(astdPMF511)  mean(astdPMF411)  mean(astdPF11)]',...
    [ mean(astdPMF22) mean(astdPMF222) mean(astdPMF322) mean(astdPMF522) mean(astdPMF422) mean(astdPF22)]',...
    [ mean(tocPMFavgOut) mean(tocPMFavgOut2) mean(tocPMFavgOut3) mean(tocPMFavgOut5) mean(tocPMFavgOut4)  mean(tocPFavgOut)]',...
    'VariableNames',{'RMSE x1','RMSE x2','ASTD 1,1','ASTD 2,2','TIME'},'RowName',...
    {'PMF','eFFT','eFST 0.01','eFST 0.005','eFST 0.001','PF bootstrap'}) %#ok<NOPTS>