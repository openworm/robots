%Given a number of curves specified by tangent angles, get the eigencurves

% Inputs: curves: an nxm matrix of curves.  n = # of curves, m = number of
% points per curve

% Outputs: the eigenvectors of the covariance matrix of the curves

%In theory, this only has to be run once for a particular family of curves, then stored.

function [vecs, vals] = familyEigenVectorsFromAngles(curves)
	%Next, compute the covariance matrix.
	covMtx = cov(curves);
	
	%Now find the eigenvectors, and we're done!
	[vecs,vals] = eig(covMtx);
    vals = diag(vals);
end