function [ image ] = generate_circle(empty, noise, size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    empty = true;
end
if nargin < 2
    noise = false;
end
if nargin < 3
    size = [640,640];
end

% Parameters:
radius = floor(size(1) / 3);
center = floor(size ./2);

% Circle with radius centered at 40,40 in an 80x80image
[xMat,yMat] = meshgrid(1:size(1),1:size(2));
distFromCenter = sqrt((xMat-center(1)).^2 + (yMat-center(2)).^2);

image = distFromCenter<=radius;

% To get just the border
if empty == true
    image = image & ~bwmorph(image,'erode',1);
end

% Add noise:
if noise == true
    image = image + (-0.1 + (0.1 + 0.1).*rand(size));
end


end

% To get just the border
% B = circleMat & ~bwmorph(circleMat,'erode',1);
% figure, imshow(B)
