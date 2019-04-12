Gamma1 = lambda^2 * sin(omegas) ./ omegas ./ (omegas.^2 - lambda^2) + (cos(omegas) + lambda^2 * sin(omegas) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
Gamma2 = lambda^2 * (1 - cos(omegas)) ./ omegas ./ (omegas.^2 - lambda^2) + (sin(omegas) + lambda^2 * (1 - cos(omegas)) ./ (omegas .* (omegas.^2 - lambda^2))) ./ (2 *(omegas.^2 - lambda^2) ./ (lambda^2 * (2 * H * tan(theta) - 1)) - 1);
xc = 0:1e-6:1;
idx = 2;
-omegas(idx)^2 .* Gamma1(idx) + lambda^2 * (sin(omegas(idx))/omegas(idx)+Gamma1(idx))+lambda^2*(H*tan(theta)-1/2) * (cos(omegas(idx))+Gamma1(idx))
-omegas(idx)^2 .* Gamma2(idx) + lambda^2 * ((1-cos(omegas(idx)))/omegas(idx)+Gamma2(idx))+lambda^2*(H*tan(theta)-1/2) * (sin(omegas(idx))+Gamma2(idx))