% plot_vetor_potential_TORUS_SLOTLESS

% plots the vector potential polynomials fitted in the gap between the yoke
% and magnets in the slotless torus machine

scatter3(reshape(design.X(:,:,1), numel(design.X(:,:,1)), []), ...
    reshape(design.Y(:,:,1), numel(design.X(:,:,1)), []), ...
    polyvaln(design.APoly(1), [reshape(design.X(:,:,1), numel(design.X(:,:,1)), []), reshape(design.Y(:,:,1), numel(design.X(:,:,1)), [])]));

hold all

scatter3(reshape(design.X(:,:,2), numel(design.X(:,:,1)), []), ...
    reshape(design.Y(:,:,2), numel(design.X(:,:,1)), []), ...
    polyvaln(design.APoly(2), [reshape(design.X(:,:,2), numel(design.X(:,:,1)), []), reshape(design.Y(:,:,2), numel(design.X(:,:,1)), [])]));

hold off

legend('lhs coil part', 'rhs coil part')