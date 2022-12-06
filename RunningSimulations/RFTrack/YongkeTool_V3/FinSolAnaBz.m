
% Z (vector), R, l in [m]
% Bc, Bz (vector) in [T]

function Bz = FinSolAnaBz (Z, R, l, Bc)
  u0_N_I = 2.0 * Bc / ( (l/2.0-0)/l/hypot(R,l/2.0-0) + (l/2.0+0)/l/hypot(R,l/2.0+0) );
  Bz = [];
  Bz = u0_N_I/2.0 * ( (l/2.0-Z)/l./hypot(R,l/2.0-Z) + (l/2.0+Z)/l./hypot(R,l/2.0+Z) ); % T
endfunction

