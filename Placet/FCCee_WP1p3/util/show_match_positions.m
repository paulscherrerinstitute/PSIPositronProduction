 
LM = [];
for i = 1:5
  A = dlmread(["Beamlines/M" int2str(i) ".tcl"]);
  LM(i) = sum(A(:,3));
endfor

LS = [];
LS(1) = 8*(0.4+0.15*2+1.5) + 0.4; 
LS(2) = 18*(0.4+0.09*2+1.5) + 0.4; 
LS(3) = 14*(0.4*2+0.075*4+1.5*2) + 0.4; 
LS(4) = 7*(0.4*2+0.05*6+1.5*4) + 0.4; 
LS(5) = 6*(0.4*2+0.1*8+1.5*6) + 0.4; 

ZMi(1) = 0;
ZMo(1) = ZMi(1) + LM(1);

ZMi(2) = ZMo(1) + LS(1);
ZMo(2) = ZMi(2) + LM(2);

ZMi(3) = ZMo(2) + LS(2);
ZMo(3) = ZMi(3) + LM(3);

ZMi(4) = ZMo(3) + LS(3);
ZMo(4) = ZMi(4) + LM(4);

ZMi(5) = ZMo(4) + LS(4);
ZMo(5) = ZMi(5) + LM(5);

ZMi
ZMo

for i = 1:5
  printf("INFO:: Match section %i located at z = [%f, %f] m\n", i, ZMi(i), ZMo(i));
endfor
