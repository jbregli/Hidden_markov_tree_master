function [] = plot_states(transform, image, layer, scale )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    figure
    subplot(2,1,1)
    imagesc(transform{image}{layer}.signal{scale})
    subplot(2,1,2)
    imagesc(transform{image}{layer}.state{scale})


end

