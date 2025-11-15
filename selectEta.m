function eta=selectEta(objTx,objRx,L)
if isa(objTx,'matlabshared.satellitescenario.GroundStation') && ~isa(objRx,'matlabshared.satellitescenario.GroundStation')
    eta=etaLinkUL(L);
elseif ~isa(objTx,'matlabshared.satellitescenario.GroundStation') && isa(objRx,'matlabshared.satellitescenario.GroundStation')
    eta=etaLinkDL(L);
else
    eta=etaLinkISL(L);
end
end