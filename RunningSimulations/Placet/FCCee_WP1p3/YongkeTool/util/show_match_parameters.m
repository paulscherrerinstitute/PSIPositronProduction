 
LM = [];
for i = 1:5
  A = dlmread(["Beamlines/M" int2str(i) ".tcl"]);
  P = [];
  for j = 1:rows(A)
    if (mod(j,2)==1)
      P(j) = A(j,3);
    else
      P(j) = A(j,5);
    endif
  endfor
  P
endfor
