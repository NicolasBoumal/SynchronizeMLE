function [C msecrb] = synchrocrb(problem)
% Returns the Cramér-Rao lower bound matrix for the synchronization problem.
%
% See the associate paper:
% http://imaiai.oxfordjournals.org/content/early/2013/09/23/imaiai.iat006.abstract
% N. Boumal, A. Singer, P.-A. Absil and V. D. Blondel,
% Cramér-Rao bounds for synchronization of rotations,
% in Information and Inference, a journal of the IMA, 2013.
%
% C is of size NxN and is such that C(i, i) is a lower bound on the
% expected squared distance between rotation i and the estimator of
% rotation i. This is valid for anchored synchronization. For anchor-free
% synchronization or 1-anchor synchronization,
% C(i, i) + C(j, j) - 2*C(i, j) is approximately a lower bound on the
% expected squared distance between the relative rotation Ri*Rj^T and the
% estimator of that relative rotation.
%
% For SO(2), there are no curvature terms in the CRB.
% For SO(3), the curvature terms are included in the CRB.
% For SO(n), n >= 4, the curvature terms are neglected in the CRB.
%
% msecrb is the CRB on the mean squared error (MSE), normalized by the
% number of rotations to estimate. If there are no anchors, the number of
% "free" rotations is N-1 and not N, because of the global rotation
% ambiguity.
%
% See also: synchromse synchrolaplacian
%
% Nicolas Boumal, UCLouvain, Jan. 16, 2013.

    n = problem.n;
    N = problem.N;
    A = problem.A;
    
    % Construct the synchronization Laplacian, where edges are weighted by
    % information quality.
    L = synchrolaplacian(problem);
    
    % Force rows and columns corresponding to anchors to zero.
    LA = L;
    LA(:, A) = 0;
    LA(A, :) = 0;
    
    % The CRB is essentially the pseudoinverse of LA.
    % There may be better ways of computing this.
    iLA = pinv(full(LA));
    
    % For n = 2, there is no curvature. For n = 3, we have an explicit
    % expression. For n >= 4, we have no expression (but the curvature
    % terms in the CRB should be negligible).
    switch n
        case 2
            C = iLA;
        case 3
            diLA = diag(diag(iLA));
            C = 9 * ( iLA - .25*(diLA*iLA + iLA*diLA) );
        otherwise
            d = nchoosek(n, 2);
            C = (d^2) * iLA;
    end
    
    % CRB on the MSE
    msecrb = trace(C) / (N-max(1, numel(A)));

end
