function [num_boxes] = Box_count(lower_bounds, pos, speed, tau_array, box_size)
% Find the number of boxes of given size required to cover attractor (Idea by Dr. Ramsden)
% 
% Inputs:
% lower_bounds = array of minimum values of x, dx/dt and tau values respectively
% pos = array of solution values in x
% speed = array of solution values in dx/dt
% tau_array = array of values of tau
% box_size = size of boxes to cover attractor with
%
% Outputs:
% num_boxes = number of boxes of size box_size neeeded to cover given attractor
%
% This idea by Dr. Ramsden is truly genius. Originally, the implementation was of O(N^2) and was
% essentially sphere counting rather than box counting, and if attempted, box couting would have
% been a very messy job in this 3D case. However, as we had a finite set of points to work with
% rather than a continuous 2 or 3D image, we could work around that. The idea now was to translate
% the entire attractor so that all coordinates in each dimension are positive (using knowledge of
% its minimum) and then divide the translated coordinates by the box size, to box them essentially,
% and finally floor them, to force all the points in a "box" into its lowermost corner (in each
% dimension). Now, to find the number of boxes covering the attractor we find the number of boxes
% containing one or more points, which can be done by finding the number of unique points in the new
% data sets we have! Without knowing exactly how Matlabs "unique" function works, we cannot
% determine the time complexity of this new implementation, but at worst case we can say it would be
% O(N^2) but with MUCH smaller coeffcients due to inbuilt compiled C functions being used. 

    % Finding lower bounds in each dimension
    pos_min = lower_bounds(1);
    spd_min = lower_bounds(2);
    tau_min = lower_bounds(3);

    % Boxing the attractor after translating it into positive quadrant (or octant in 3D I guess?)
    pos_boxed_points = floor((pos - pos_min) / box_size);
    spd_boxed_points = floor((speed - spd_min) / box_size);
    tau_boxed_points = floor((tau_array - tau_min) / box_size);

    % Finding number of unique points left and hence number of boxes required
    boxed_coords = transpose([pos_boxed_points; spd_boxed_points; tau_boxed_points]);
    unique_points = unique(boxed_coords, 'rows');
    num_boxes = length(unique_points);

end