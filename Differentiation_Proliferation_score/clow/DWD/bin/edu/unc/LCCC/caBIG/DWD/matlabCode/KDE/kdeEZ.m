function kdeEZ()

% Everett Zhou 3/31/2005

% This function does:
%   Read the original data from a text file into array.
%   Pass one column each time to kdeSMEZ to get xgrid and kde and merge
%   together
%   Save the mrege array into a text file.

[dataAll] =ReadKDEInput;
[n d]=size (dataAll);
xgridkde = [];

for i = 1:d;
    data =dataAll(:,[i]);
    [temp_kde,temp_xgrid] = kdeSMEZ(data);
    xgridkde = [xgridkde temp_xgrid  temp_kde];
end ;

WriteKDEOutput (xgridkde);