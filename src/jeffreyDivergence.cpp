
double hellingerDistance (RealVector mu1, RealMatrix S1, RealVector mu2, RealMatrix S2) {
	RealVector u = mu1 - mu2;
	RealMatrix iS1 = inv (S1);
	RealMatrix iS2 = inv (S2);
	RealMatrix G = 0.5 * S1 + 0.5 * S2;
	RealMatrix iG = inv (G);

	
	double det1, sign;
	log_det(det1, sign, S1);

	double det2, sign;
	log_det(det2, sign, S2);
	
	
	double rho = 1.0/sqrt (det(G)) * powf (det1, 0.25) * powf(det2, 0.25) * exp (-0.125 * t(u)*iG*u)
	return (sqrt (1 - rho));
}




# jeffreys=symmetric kl-divergence, but it is not a metric.
jeffreyDivergence <- function (mu1, S1, mu2, S2)
{
	d = dim(S1)[1] 
	u = mu1 - mu2
	iS1 = inv(S1)
	iS2 = inv(S2)
	psi = iS1 + iS2
	div = 1/2 * t(u) %*% psi %*% u + 1/2 * tr(iS1 %*% S2 + iS2 %*% S1) - d
	
	return (div)
}
