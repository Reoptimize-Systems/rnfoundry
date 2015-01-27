function Results


    % Always loop 10 simulation runs to smooth out "noise"
    nIters = 10;
    ProbDims = 2.^(0:18);
    
    Ave = zeros(length(ProbDims),5);
    for i = 1:length(ProbDims)
         Y = InterpTests(ProbDims(i),nIters);
         
         Ave(i,:) = mean(Y);
         
    end
    
    semilogx(ProbDims,Ave(:,1),'-ro',ProbDims,Ave(:,2),'-go',ProbDims,Ave(:,3),'-bo',...
        ProbDims,Ave(:,4),'-ko',ProbDims,Ave(:,5),'-mo','LineWidth',2);
    legend('1-D (ppuval)','1-D (ppmval)','2-D (ppmval)','3-D (ppmval)','4-D (ppmval)');
    xlabel('Problem dimension');
    ylabel('Performance gain');
end