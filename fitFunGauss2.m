function [err] = fitFunGauss2(q,x,y)

% q(1) = alpha, q(2) = mu, q(3) = sigma
lambda = q(1) .* exp(-(((x - q(2)).^2) / (2.*q(3).^2)));
err = -sum(log(poisspdf(y,lambda)));   %log likelihood