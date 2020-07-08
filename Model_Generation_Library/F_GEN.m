function F=F_GEN(L1,L2,the1,the2)

F1 = exp(-1i*pi*(0:L1-1)'*(the1))/sqrt(L1);
F2 = exp(-1i*pi*(0:L2-1)'*(the2))/sqrt(L2);
F=kron(F2,F1);

end