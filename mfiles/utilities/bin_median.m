function [Aavg] = bin_median(A,dt)
% bin_median gives average of data over every dt points
%
% [Aavg] = bin_median(A,dt)
%
% If A is a vector, bin_median will simply calculate median over dt points..
%
% For instance, if you have 10 minute data and you want hourly medians, 
% set dt = 6. This code then finds the median over the six points in every hour.
%
% if A is a structure. A.time MUST be one of the fields in A.
% dt is the number of points you want to average (median) into each time step of the
%       new vector.
% All variabies in A with the same lenth as time will be bin averaged (median).



if isstruct(A)      % case where A is a structure
    
    vars = fields(A);

    N = length(A.time);
    NN = floor(N/dt);
    leftover = mod(N,dt);

    for vv = 1:length(vars)
        if length(size(A.(vars{vv}))) > 2
            disp('Matrix is too big! It needs to be 1D or 2D.')
            return
        end
        if length(A.(vars{vv})) == length(A.time)
            if size(A.(vars{vv}),1) == 1 | size(A.(vars{vv}),2) == 1
                Aavg.(vars{vv}) = nanmedian(reshape(A.(vars{vv})(1:N-leftover),dt,NN),1);
            % cases where A is a 2-D matrix, not an array    
            elseif size(A.(vars{vv}),1) == length(A.time)
                for ii = 1:size(A.(vars{vv}),2)
                    Aavg.(vars{vv})(:,ii) = nanmedian(reshape(A.(vars{vv})(1:N-leftover,ii),dt,NN),1);
                end
            elseif size(A.(vars{vv}),2) == length(A.time)
                for ii = 1:size(A.(vars{vv}),1)
                    Aavg.(vars{vv})(ii,:) = nanmedian(reshape(A.(vars{vv})(ii,1:N-leftover),dt,NN),1);
                end
            end
        end
    end
    
else            % case where A is a vector
    
    N = length(A);
    NN = floor(N/dt);
    leftover = mod(N,dt);
    
    [ii,jj] = size(A);
    
    if ii == 1 | jj == 1
        Aavg = nanmedian(reshape(A(1:N-leftover),dt,NN),1);
    else
        disp('Error in bin_median. A must be a 1xN or Nx1 vector')
    end

end

