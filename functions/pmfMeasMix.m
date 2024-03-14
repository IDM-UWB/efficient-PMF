function [measPdf] = pmfMeasMix(predGrid,nz,k,z,V,predDensityProb,predGridDelta,hfunct)
%PMFMEAS - measurement update

    predThrMeasEq = hfunct(predGrid,zeros(nz,1),k); %Prediction density grid through measurement EQ
    pom = z'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
    filterDensityNoNorm =  pdf(V.pdf,pom).*predDensityProb; % Filtration density unnormalized
   
    measPdf = filterDensityNoNorm/sum(prod(predGridDelta)*filterDensityNoNorm); %Normalization
end

