
double jeffreyRiemannDivergence (RealVector mu1, RealMatrix S1, RealVector mu2, RealMatrix S2) {
	RealVector u = mu1 - mu2;
	RealMatrix iS1 = inv (S1);
	RealMatrix iS2 = inv (S2);

	RealMatrix psi = iS1 + iS2;
	
	v = geigen(S1, S2)$values
	dR = sqrt( tr( log(v)^2 ))
	
	double div = sqrt( 0.5 * t(u) %*% psi %*% u) + dR
	return (div);
}

