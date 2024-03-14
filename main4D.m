%%
%   Efficient Point-Mass Filter for 4D.
%   author: pesslovany@gmail.com
%   Papers reference will be added after publication

%% Parameters and system simulation
clc
clear
close all

load("mapTAN.mat")

method = 'fft'; % std- efficient pmf using convolution(can be faster, depends on the number of points, dimensions, and hardware)
% fft - efficient pmf using FFT

format shortG

% System parameters
nx = 4; % state dimension
nz = 1; % measurement dimension
kf = 10; % final time

Ts = 1; % time step
q = 10; % noise parameter
turn = deg2rad(30); % turn rate

Q = q*[(2*(turn*Ts-sin(turn*Ts))/turn^3) (1-cos(turn*Ts))/turn^2 0 ((turn*Ts-sin(turn*Ts))/turn^2);
    (1-cos(turn*Ts))/turn^2 Ts -((turn*Ts-sin(turn*Ts))/turn^2) 0;
    0 -((turn*Ts-sin(turn*Ts))/turn^2) (2*(turn*Ts-sin(turn*Ts))/turn^3) (1-cos(turn*Ts))/turn^2;
    ((turn*Ts-sin(turn*Ts))/turn^2) 0 (1-cos(turn*Ts))/turn^2 Ts]; % system noise
Q = Q([1 3 2 4],[1 3 2 4]);% change to state [p_x p_y v_x v_y] - my habbit

invQ = inv(Q);
R = 1; % measurement noise covariance for both modes
invR = inv(R);

% PMF parameters
Npa = 21; % number of points per axis
N = Npa^nx; % number of points - total
sFactor = 4; % scaling factor (number of sigmas covered by the grid)

meanV = [0 20]; % Mean values of components of meas noise
wV = [0.5 0.5]; % weights

V.pdf = gmdistribution(meanV',R,wV); % Meas noise pdf


MC = 10; % monte carlo simulations
for mc = 1:1:MC
    mc
    % Initial condition - Gaussian
    meanX0 = [mean(map_m.x,"all"); mean(map_m.y,"all"); 50; 50];% initial cond mean value
    varX0 = [90 0 0 0;
        0 160 0 0;
        0 0 5 0;
        0 0 0 5]; % initial cond variance

    % PF parameters
    noPart = 1200000; % number of particles for PF
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF


    x = zeros(nx,kf); % state
    x(:,1) = mvnrnd(meanX0,varX0); % initial estimate


    % Dynamics - known turn rate
    F = [1 sin(turn*Ts)/turn 0 (cos(turn*Ts)-1)/turn;
        0 cos(turn*Ts)/turn 0 -sin(turn*Ts);
        0 (1-cos(turn*Ts))/turn 1 sin(turn*Ts)/turn;
        0 sin(turn*Ts) 0 cos(turn*Ts)];
    F = F([1 3 2 4],[1 3 2 4]);


    ffunct = @(x,w,k) F*x + w; % state equation
    hfunct = @(x,v,k) interp2(map_m.x,map_m.y,map_m.z,x(1,:),x(2,:)) + v; % measurement equation

    % System simulation
    w = mvnrnd(zeros(nx,1),Q,kf-1)'; % system noise generation

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
    measGridOld = predGrid;

    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)

    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index


    for k = 1:1:kf
        tic
        %% Measurement update
        [measPdf] = pmfMeasMix(predGrid,nz,k,z(:,k),V,predDensityProb,GridDelta(:,k),hfunct); % Measurement update

        % Measurement mean and var
        measMean(:,k) = predGrid*measPdf*prod(GridDelta(:,k)); % Measurement update mean
        covariance = zeros(nx);
        chip_ = (predGrid-measMean(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar(:,:,k) = chip_w*chip_' * prod(GridDelta(:,k)); % Measurement update variance


        % Expected pred mean/var using kf
        predMeanEst = F*measMean(:,k);
        predVarEst = F*measVar(:,:,k)*F' + Q;

        % Set-up the meas grid on which to interpolate on the basis of wanted
        % pred grid sigma
        [eigVect,eigVal] = eig(predVarEst); % eigenvalue and eigenvectors, for setting up the grid
        eigVal = diag(eigVal);

        gridBoundWant = sqrt(eigVal)*sFactor; % Wanted boundaries of pred grid
        gridBoundWantCorners = boxvertex(nx,gridBoundWant); % Corners for 4D cube normal cube
        gridBoundWantCorners = (gridBoundWantCorners'*eigVect)' + predMeanEst; % Wanted corner of predictive grid
        gridBoundWantCorners = inv(F)*gridBoundWantCorners; %#ok<*MINV> % Back to filtering space
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
        % actually measurement grid  from last step - gridDimOld) - therefore we are able to use griddedInterpolation which is
        % volumes faster than scatteredInterpolation. In esence we are
        % calculating the interpolation in other space.

        % Interpolation
        Fint = griddedInterpolant(gridDimOld,reshape(measPdf,Npa,Npa,Npa,Npa),"linear","nearest");
        if k == 1
            filtGridInterpInvTrsf = measGridNew; % first time step
        else
            filtGridInterpInvTrsf = inv(F)*measGridNew; % other time steps
        end
        measPdf = Fint(filtGridInterpInvTrsf(1,:),filtGridInterpInvTrsf(2,:),filtGridInterpInvTrsf(3,:),filtGridInterpInvTrsf(4,:))';

        gridDimOld = gridDim;

        %% Time update

        % Pred Grid
        predGrid = F*measGridNew; % Predictive grid
        GridDelta(:,k+1) = F*GridDelta(:,k); % Grid step size

        % Precalculated TPM
        pom = (predGrid(:,halfGrid)'-(predGrid)');
        TPMrow = ((exp(sum(-0.5*pom*invQ.*pom,2)))/predDenDenomW)';% Middle row of transition matrix
        TPMrowCub = reshape(TPMrow,Npa,Npa,Npa,Npa); % Transition matrix middle row in physical space


        % ULTRA FAST PMF
        filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
        filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa,Npa,Npa,Npa); % reshape to physical space
        switch method
            case 'std' % efficient pmf using convolution (can be faster, depends on the number of points, dimensions, and hardware)
                predDensityProb2cub = convn(filtDenDOTprodDeltasCub,TPMrowCub,"same");
            case 'fft' % efficient pmf using FFT
                lfftfun = @(l) 2^nextpow2(l); % Power of two for easier FFT
                dims = 1:1:nx;

                for dim=dims % Over dimensions
                    % compute the FFT length
                    l = lfftfun(Npa+Npa-1);
                    TPMrowCub = fft(TPMrowCub,l,dim); % FFT of transition density matrix middle row in dim dimenson
                end
                predDensityProb2cub = convnfft(filtDenDOTprodDeltasCub, TPMrowCub, Npa); % convolution
        end
        predDensityProb2cub = predDensityProb2cub./(sum(predDensityProb2cub,"all")*prod(GridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)
        predDensityProb = reshape(predDensityProb2cub,N,1); % Reshape back to computation space

        tocPMF(k) = toc; % Time evaluation


        %% Partcle filter - to compare
        tic

        %Measurement - Update
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); % Prediction density grid through measurement EQ
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom); % Weights
        w = w/sum(w); % Normalization of weights
        xEst(:,k) = ksiPrior*w; % mean filtering estimate
        varEst(:,k) = ksiPrior.^2*w - xEst(:,k).^2; % diagonal of filtering covariance matrix

        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPart);
        I = binarySearch(cumW, randomN, 'first');

        % Time - Update
        ksi = ksiPrior(:,I);
        ksiPrior = F*ksi + mvnrnd(zeros(1,nx),Q,noPart)';

        tocPF(k) = toc;

    end


    % Evaluation
    rmsePMF(:,mc) = sqrt(mean((x-measMean).^2,2)); %#ok<*SAGROW> 
    astdPMF11(mc) = sqrt(mean(measVar(1,1,:)));
    astdPMF22(mc) = sqrt(mean(measVar(2,2,:)));
    astdPMF33(mc) = sqrt(mean(measVar(3,3,:)));
    astdPMF44(mc) = sqrt(mean(measVar(4,4,:)));


    rmsePF(:,mc) = sqrt(mean((x-xEst).^2,2));
    astdPF11(mc) = sqrt(mean(varEst(1,:)));
    astdPF22(mc) = sqrt(mean(varEst(2,:)));
    astdPF33(mc) = sqrt(mean(varEst(3,:)));
    astdPF44(mc) = sqrt(mean(varEst(4,:)));

    tocPMFavg(:,mc) = mean(tocPMF);
    tocPFavg(:,mc) = mean(tocPF);
end

%% Evaluation
rmsePMFout = mean(rmsePMF,2);
rmsePFout = mean(rmsePF,2);
tocPMFavgOut = mean(tocPMFavg,2);
tocPFavgOut = mean(tocPFavg,2);


T2 = table([ rmsePMFout(1) rmsePFout(1)]',...
    [ rmsePMFout(2) rmsePFout(2)]',...
    [ rmsePMFout(3) rmsePFout(3)]',...
    [ rmsePMFout(4) rmsePFout(4)]',...
    [ mean(astdPMF11) mean(astdPF11)]',...
    [ mean(astdPMF22) mean(astdPF22)]',...
    [ mean(astdPMF33) mean(astdPF33)]',...
    [ mean(astdPMF44) mean(astdPF44)]',...
    [ mean(tocPMFavgOut) mean(tocPFavgOut)]',...
    'VariableNames',{'RMSE x1','RMSE x2','RMSE x3','RMSE x4','ASTD 1','ASTD 2','ASTD 3','ASTD 4','TIME'},'RowName',...
    {'PMF','PF bootstrap'}) %#ok<NOPTS>
