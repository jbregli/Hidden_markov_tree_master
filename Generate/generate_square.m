function [ image ] = generate_square(empty, noise, translate, rotate, size)
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

    mask = floor(size./2);

    fill_g = floor((size-mask) ./ 2);
    fill_d = size - fill_g - mask;

    % to add the square, make the top left quarter white
    %# by setting the pixel values to true (i.e. 1)

    image(fill_g(1)+1:(mask(1)+fill_d(1)),...
            fill_g(2)+1:(mask(2)+fill_d(2))) = 1;

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
