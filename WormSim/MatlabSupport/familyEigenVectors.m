%Given a number of curves, compute the eigencurves

% Inputs: curves: a 2*nxm matrix of curves.  curves are made up of two adjacent rows.  Row i is the x(s) for curve i/2, Row i+1 is the y(s) for that curve

% Outputs: the eigenvectors

%In theory, this only has to be run once for a particular family of curves, then stored.

function [vecs, vals] = familyEigenVectors(curves)
	
	sizeCurves = size(curves);
	nCurves = sizeCurves(1)/2;
	lCurves = sizeCurves(2);
	
	%First, get the theta vector for each curve
	thetas = zeros(nCurves,lCurves-1); %due to numerical differentiation, there is one less angle than there is number of points
	
	%TODO: this is strictly speaking not very general.  It assumes we are using a simple first differences approach.  The entire library should be made much more
	%robust to the method of differentiation
	
	for i = 1:nCurves
		thetas(i,:) = curveAngles(curves(2*i-1:2*i,:));
	end
	
	%Next, compute the covariance matrix.
	covMtx = cov(thetas);
	
	%Now find the eigenvectors, and we're done!
	[vecs,vals] = eig(covMtx);
    vals = diag(vals);
end