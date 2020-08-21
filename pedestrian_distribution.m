function [res] = pedestrian_distribution (N)

% pedestrian distribution along the sidewalk
%inter-pedestrian distance
Total = ((419+236)*2)-((413+230)*2);

res = Total/N;

end