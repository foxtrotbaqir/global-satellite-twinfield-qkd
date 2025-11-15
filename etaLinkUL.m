function eta=etaLinkUL(L)
lam=1550e-9; Dtx=1.0; Drx=0.3; alpha=0.003; sigma=3e-7;
wL=(lam/(pi*Dtx))*L;
eta=(Drx/(2*wL))^2*exp(-alpha*(L/1e3))*exp(-2*(sigma*L)^2/wL^2)*0.8*0.3;
end
