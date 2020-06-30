function [Amo] = monthly_median(A)
% monthly_median loops through every month in a time series,
% and finds the median of each month's data.
%
% [Amo] = monthly_median(A)
%
% A is your input structure. It MUST have a field A.time
% Amo is the median in each month for all of the fields in A
%
% Note: If you data has dt < 1 day, this code may take a very long time to
% run. It is MUCH quicker to first use bin_meidan to convert high
% frequency data to daily averages BEFORE running monthly_average.

 
count = 0;

year1 = str2num(datestr(A.time(1),'yyyy'));
year2 = str2num(datestr(A.time(end),'yyyy'));

vars = fields(A);

for yy = year1:year2
%     disp(['year = ' num2str(yy)])
    for mm = 1:12
%         disp(['year = ' num2str(yy) ', month = ' num2str(mm)])
        clear ind timetest
        ind = str2num(datestr(A.time,'mm')) == mm & ...
            str2num(datestr(A.time,'yyyy')) == yy;

        timetest = nanmean(A.time(ind));
        if ~isnan(timetest) % ignores part of first year before record starts
            count = count + 1;            
            Amo.time(count) = timetest;

            for vv = 1:length(vars)
                % check that variable isn't the time vector
                if ~strcmp(vars{vv},'time')
                    % check that variable has the same length as A.time
                    if length(A.(vars{vv})) == length(A.time)   
                        % one dimensional case
                        if size(A.(vars{vv}),1) == 1 | size(A.(vars{vv}),2) == 1
                            % check that there are at least 10 good days in the month
                            if sum(~isnan(A.(vars{vv})(ind))) >= 10
                                Amo.(vars{vv})(count) = ...
                                    nanmedian(A.(vars{vv})(ind));
                            else
                                Amo.(vars{vv})(count) = NaN;
                            end
                        % two dimensional case 1    
                        elseif size(A.(vars{vv}),1) == length(A.time)
                            for ii = 1:size(A.(vars{vv}),2)
                                % check that there are at least 10 good days in the month
                                if sum(~isnan(A.(vars{vv})(ind,ii))) >= 10
                                    Amo.(vars{vv})(count,ii) = ...
                                        nanmedian(A.(vars{vv})(ind,ii));
                                else
                                    Amo.(vars{vv})(count,ii) = NaN;
                                end
                            end
                        % two dimensional case 2    
                        elseif size(A.(vars{vv}),2) == length(A.time)
                            for ii = 1:size(A.(vars{vv}),1)
                                % check that there are at least 10 good days in the month
                                if sum(~isnan(A.(vars{vv})(ii,ind))) >= 10
                                    Amo.(vars{vv})(ii,count) = ...
                                        nanmedian(A.(vars{vv})(ii,ind));
                                else
                                    Amo.(vars{vv})(ii,count) = NaN;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

end

