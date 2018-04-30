function Z = randlangevinoutliers(n, N, kappa, p)
% Z is an n-by-n-by-N matrix such that each n-by-n matrix is distributed
% with probability p according to a Langevin distribution of mean identity
% and of concentration kappa (isotropic), and with probability 1-p is just
% picked uniformly at random on SO(n).
%
% Nicolas Boumal, Oct. 19, 2011.

    Z = zeros(n, n, N);
    outliers = rand(N, 1) > p;
    noutliers = nnz(outliers);
    Z(:, :,  outliers) = randrot(n, noutliers);
    Z(:, :, ~outliers) = randlangevin(kappa(~outliers), n, N-noutliers);

end
