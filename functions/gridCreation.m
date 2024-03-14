function [predGrid, predGridDelta, gridDim] = gridCreation(xp_aux,Pp_aux,sFactor,nx,Npa)
% gridCreation - creates grid with boundaries aligned with state space axis

    pom = sqrt(diag(Pp_aux))*sFactor; % Grid size
    maxF = (pom + xp_aux); 
    minF = (-pom + xp_aux); 

    for ind3 = 1:1:nx %Creation of propagated grid
        gridDim{ind3,1} = linspace(minF(ind3),maxF(ind3),Npa);%New grid axis
        gridStep(ind3,1) = abs(gridDim{ind3,1}(1)-gridDim{ind3,1}(2));%Grid step
    end
    predGrid = combvec(gridDim); %New grid
    predGridDelta = gridStep; %New Grid step size


end

