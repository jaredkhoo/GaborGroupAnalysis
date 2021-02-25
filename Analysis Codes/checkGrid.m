% checkGrid
%
% Takes the range of x and y values used to generate a meshgrid, the grid
% resulting from applying a function to the meshgrid, and a target value.
% Returns complete = true if all values on the edges of the grid are
% greater than the target. Otherwise, incrementally adjusts the range of x 
% or y values and returns complete = false
function [slope_range,int_range,complete] = checkGrid(chi_grid,slope_range,int_range,target)
    complete = true;
        
    % Reduce the slope minimum if any value in column 1 is < target
    if(min(chi_grid(:,1)) < target)
        slope_range(1) = slope_range(1) - 1;
        complete = false;
    end
    
    % Increase the slope maximum if any value in column 100 is < target
    if(min(chi_grid(:,100)) < target)
        slope_range(2) = slope_range(2) + 1;
        complete = false;
    end
    
    % Reduce the intercept minimum if any value in row 1 is < target
    if(min(chi_grid(1,:)) < target)
        int_range(1) = int_range(1) - 5;
        complete = false;
    end

    % Increase the intercept maximum if any value in row 100 is < target
    if(min(chi_grid(100,:)) < target)
        int_range(2) = int_range(2) + 5;
        complete = false;
    end
end