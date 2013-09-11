function h = plotpolymodel(polymodel, varargin)
% plotpolymodel: plots a 1D or 2D polymodel within the data range to which
% it is fitted
%
% Syntax
% 
% plotpolymodel(polymodel);
% plotpolymodel(..., 'Parameter', 'Value', ...)
%     
% Input
%
% polymodel - A 1D or 2D polymodel as produced by polyfitn
%
% If 1D, the polymodel will be evaluated at 1000 points within the range of
% data to which it was originally fitted. To specify the number of data
% points, the Parameter-Value pair denoted by 'Points' can be used.
% 
% If 2D, the polymodel will be evaluated on a 25 x 25 linearly spaced grid
% in the range of data to which it is fitted. To specify the size of the
% grid set the p-v pair denoted by 'GridSize', this will be a scalar value
% which denotes the number of points per dimension of the grid, e.g.
%
% plotpolymodel(polymodel, 'GridSize', 50)
%
% will plot the polymodel over a 50 x 50 grid. 
%
% By default a contour plot is produced for a 2D polymodel. Other plot
% types are available through the PlotType2D p-v pair. The PlotType2D is a
% string with following possible values:
%
%   'contour' - (default) contour plot
%   'scatter' - 3D scatter plot with colour bar
%   'surf' - surface plot with colour bar
%   'surfc' - combined surface and contour plot with colour bar
%
% Output
%
% h is the handle to the figure which is produced
%


    Inputs.PlotType2D = 'contour';
    Inputs.GridSize = 25;
    Inputs.Points = 1000;
    Inputs.PlotOptions = {};
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if ischar(Inputs.PlotOptions)
        Inputs.PlotOptions = {Inputs.PlotOptions};
    end
    
    if size(polymodel.DataRange, 2) > 2
        error('Only 1D and 2D polymodels can be visualised')
    end
    
    h = figure;
    
    if size(polymodel.DataRange, 2) == 1
        
        % 1D polymodel
        xrange = linspace(polymodel.DataRange(1,1), polymodel.DataRange(2,1), Inputs.Points);
        
        y = polyvaln(polymodel, xrange');
        
        plot(xrange, y, Inputs.PlotOptions{:});
        
        xlabel(polymodel.VarNames{1});
        ylabel('Poly Value');
        
    else 
        % 2D polymodel
        xrange = linspace(polymodel.DataRange(1,1), polymodel.DataRange(2,1), Inputs.GridSize);
        yrange = linspace(polymodel.DataRange(1,2), polymodel.DataRange(2,2), Inputs.GridSize);
        [X,Y] = meshgrid(xrange,yrange);
        
        Z = reshape(polyvaln(polymodel, [X(:), Y(:)]), size(X));
        
        switch Inputs.PlotType2D
            
            case 'contour'
                
                contour(X, Y, Z, Inputs.PlotOptions{:});
                colorbar;
                
            case 'scatter'
                
                if isempty(Inputs.PlotOptions)
                    scatter3(X(:), Y(:), Z(:), 40, Z(:), '+');
                else
                    scatter3(X(:), Y(:), Z(:), 40, Z(:), Inputs.PlotOptions{:});
                end
                colorbar;
                
            case 'surf'
                
                surf(X,Y,Z, Inputs.PlotOptions{:})
                colorbar;
                
            case 'surfc'
                
                surfc(X,Y,Z, Inputs.PlotOptions{:})
                colorbar;
                
            otherwise
                error('Unsupported 2D polymodel plot type');
        end
        
        xlabel(polymodel.VarNames{1});
        ylabel(polymodel.VarNames{2});
        zlabel('Poly Value');
        
    end
    
    

end