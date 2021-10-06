function addplot(varargin)

    hold on;
    plot(gca, varargin{:});
    hold off;
    
end