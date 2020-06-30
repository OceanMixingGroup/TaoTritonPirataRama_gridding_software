function [ time ] = taotime
% Create the same time vector for all TaoTritonPirataRama analysis

% typical chipod times
starttime = datenum(2005,09,23,04,05,00);

% endtime = datenum(2017,12,31,23,55,00);
% endtime = datenum(2019,07,30,23,55,00);
% endtime = datenum(2020,01,01,00,00,00);
endtime = datenum(2020,06,01,00,00,00);

% dates for Bill's analysis that goes back to Jan 1, 1990
% starttime = datenum(1990,01,01,00,00,00);
% endtime = datenum(2020,05,01,00,00,00);


jump = datenum(0,0,0,0,10,0);
time = starttime:jump:endtime;


end

