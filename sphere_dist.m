function d1 = sphere_dist(lat1, lon1,lat2, lon2)
% format: d1 =sphere_dist(latlon1,latlon2)
% Distance:
% d1: distance based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
% --Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of destination point [lat lon]
%
% --Outputs:
%   d1: distance calculated by Haversine formula

% First version: 14 April 2017
%--------------------------------------------------------------------------

deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sind((deltaLat)/2)^2 + cosd(lat1)*cosd(lat2) * sind(deltaLon/2)^2;
d1=2*atan2(sqrt(a),sqrt(1-a));    %Haversine distance


end