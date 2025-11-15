function eta=etaLinkISL(L)
lam=1550e-9; Dtx=0.3; Drx=0.3; sigma=1e-7;
wL=(lam/(pi*Dtx))*L;
eta=(Drx/(2*wL))^2*exp(-2*(sigma*L)^2/wL^2)*0.8*0.8;
end
