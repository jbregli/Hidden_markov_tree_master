function [ signal ] = generate_edge(noise, size, edge_pos)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    noise = false;
end
if nargin < 2
    size = 1080;
end
if nargin < 3
    edge_pos = floor(4 * size / 10);
end


% create an empty image (all zeros)
signal = zeros(size,1);

% to add the square, make the top left quarter white
%# by setting the pixel values to true (i.e. 1)

signal(edge_pos:end) = 1;

% Add noise:
if noise == true
    signal = signal + (-0.1 + (0.1 + 0.1).*rand(size,1));
end

end
