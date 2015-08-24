function [] = plot_states(set_S, image, layer, scale )
% plot_states: PLOT ST COEFFICIENT AND STATES AT A GIVEN LAYER AND SCALE

    figure
    subplot(2,1,1)
    imagesc(set_S{image}{layer}.signal{scale})
    subplot(2,1,2)
    imagesc(set_S{image}{layer}.state{scale})
end

