function znodemat = unifyznodes_FM(zguidecon, zsupp, zframe, tolerance)
% unifyznodes_FM: creates a categorised list of nodes in the z axis for the
% creation of a finite element mesh of a flat profile linear generator

    zsupp = [zsupp, ones(numel(zsupp), 1),  zeros(numel(zsupp), 2)];
    
    zguidecon = [zguidecon, zeros(numel(zguidecon), 1),  ones(numel(zguidecon), 1), zeros(numel(zguidecon), 1)];
    
    zframe = [zframe, zeros(numel(zframe), 1),  zeros(numel(zframe), 1), ones(numel(zframe), 1)];
    
    znodemat = [zsupp; zguidecon; zframe];
    
    znodemat = sortrows(znodemat, 1);
    
    overlapfound = false;
    firstrun = true;
    
    while overlapfound || firstrun

        firstrun = false;
        
        for i = 2:size(znodemat, 1)

            if iswithin(znodemat(i), znodemat(i-1), tolerance)
                
                % Merge the z positions we have found which overlap
                overlapfound = true;
                
                if i > 2
                    startmat = znodemat(1:i-2,:);
                else
                    startmat = [];
                end
                
                mergednode = [znodemat(i,1), znodemat(i,2:end) + znodemat(i-1,2:end)];
                
                if size(znodemat, 1) > i
                    endmat = znodemat(i+1:end,:);
                else
                    endmat = [];
                end
                
                znodemat = [startmat; mergednode; endmat];
                
                break;

            else
                overlapfound = false;
            end

        end

    end
    
end