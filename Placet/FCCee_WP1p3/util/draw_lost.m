 
load ("output/Loss.dat");

plot(A_LOSS(:,1), A_LOSS(:,2));
ylim([0,1e-2]);

lost = sum(A_LOSS(:,2));

printf("INFO:: Particles lost: %.1f%%\n", lost*100);

text(20,0.008, sprintf("Total particles lost: %.1f%%\n", lost*100));

%% show matching sections
source("util/show_match_positions.m");
for i = 1:5
  rectangle ("Position", [ZMi(i), 0.005, ZMo(i)-ZMi(i), 0.002], "EdgeColor", "r");
endfor

pause();
