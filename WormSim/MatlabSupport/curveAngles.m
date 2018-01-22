%Function to compute the angles of a curve

%TODO: better documentation

function theta = curveAngles(curve)

	%get tangent vectors
	tangents = curveTangents(curve);
	
	%get tangent angles
	theta = tangentAngles(tangents);

end