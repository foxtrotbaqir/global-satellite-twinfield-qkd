function skr=computePathSKR_asym(nodes,t0,satObjs,gs1,gs2,satNames)
nh=numel(nodes)-2; skv=zeros(1,nh);
for p=1:nh
    A=nodes{p};B=nodes{p+1};C=nodes{p+2};
    objA=getObj(A,satObjs,gs1,gs2,satNames);
    objB=getObj(B,satObjs,gs1,gs2,satNames);
    objC=getObj(C,satObjs,gs1,gs2,satNames);
    pA=getECEF(objA,t0); pB=getECEF(objB,t0); pC=getECEF(objC,t0);
    [L1,~]=rangeangle(pB,pA); [L2,~]=rangeangle(pB,pC);
    eta1=selectEta(objA,objB,L1); eta2=selectEta(objC,objB,L2);
    skv(p)=tfqkd_skr_asym(eta1,eta2,t0,gs1);
end
skr=min(skv);
end