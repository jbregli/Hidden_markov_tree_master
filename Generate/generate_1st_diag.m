function [ image ] = generate_1st_diag(empty, size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    empty = true;
end
if nargin < 2
    size = [640,640];
end

if empty == true
    image = tril(ones(size),-1) + triu(ones(size),1);
else
    image = tril(ones(size),-1);
end

end
