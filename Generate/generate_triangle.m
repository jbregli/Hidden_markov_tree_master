function [ image ] = generate_triangle(empty, noise, translate, rotate, size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 1
        empty = true;
    end
    if nargin < 2
        noise = false;
    end
    if nargin < 3
        translate = false;
    end
    if nargin < 4
        rotate = false;
    end
    if nargin < 5
        size = [640,640];
    end

    % create an empty image (all zeros)
    image = zeros(size);
    
    x_size = size(1);
    y_size = size(2);
        
    tx_size = floor(x_size ./ 3);
    ty_size = floor(y_size ./ 3);

    x_pos= floor(x_size/2);
    y_pos= floor(y_size/2);
    
    xCoords = [x_pos-tx_size x_pos+tx_size x_pos];
    yCoords = [y_pos y_pos y_pos+ty_size];
    mask = poly2mask(xCoords, yCoords, size(1), size(2));
    image(mask) = 1; % or whatever value you want.

    
    % To get just the border
    if empty == true
        image = image & ~bwmorph(image,'erode',1);
    end

    % Add noise:
    if noise == true
        image = image + (-0.05 + (0.05 + 0.05).*rand(size));
    end

    % Translate:
    if translate == true
        xmax = max(floor(size(1)/20), 5);
        ymax = max(floor(size(2)/20), 5);
        xTrans = randi([-xmax,xmax]);
        yTrans = randi([-ymax, ymax]);
        image =  circshift(image,[xTrans, yTrans]);
    end

    % Rotate:
    if rotate == true
        degRot = randi([-50,50]);
        image =  imrotate(image, degRot);
    end
end
