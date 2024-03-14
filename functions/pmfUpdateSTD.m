function predDensityProb = pmfUpdateSTD(measGrid,measPdf,predGridDelta,ffunct,predGrid,nx,k,invQ,predDenDenomW,N)
%standard time-update for PMF

    filtDenDOTprodDeltas = (measPdf*prod(predGridDelta(:,k))); % measurement PDF * measurement PDF step size

    gridNext = ffunct(measGrid,zeros(nx,1),k); % Old grid through dynamics
    for ind2 = 1:1:N %Over number of state of prediction grid
        pom = (predGrid(:,ind2)'-(gridNext)'); % distance
        predDensityProb(ind2,1) = ((exp(sum(-0.5*pom*invQ.*pom,2)))/predDenDenomW)'*filtDenDOTprodDeltas; %One point of predictive density
    end
    predDensityProb = predDensityProb./(sum(predDensityProb)*prod(predGridDelta(:,k+1)))'; % Normalizaton (theoretically not needed)

    
end

