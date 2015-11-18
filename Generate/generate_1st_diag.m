function [ image ] = generate_1st_diag(empty, size, translate)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 1
        empty = true;
    end
    if nargin < 2
        size = [640,640];
    end
    if nargin < 3
        translate = false;
    end
    if empty == true
        image = tril(ones(size),-1) + triu(ones(size),1);
    else
        image = tril(ones(size),-1);
    end

    % Translate:
    if translate == true
        xmax = max(floor(size(1)/10), 5);
        ymax = max(floor(size(2)/10), 5);
        xTrans = randi([-xmax,xmax]);
        yTrans = randi([-ymax, ymax]);
        image =  circshift(image,[xTrans, yTrans]);
    end


end
