% Add additional inputs after the given ones if you want to
% Example:
% your_input = 1;
% eskf_map_handle = @(sensor) eskf1(sensor, your_input);
%
% We will only call eskf_map_handle in the test function.
% Note that this will only create a function handle, but not run the function

eskf_map_handle = @(sensor) eskf_map(sensor);
