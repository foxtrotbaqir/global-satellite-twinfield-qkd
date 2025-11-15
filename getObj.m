function obj=getObj(name,satObjs,gs1,gs2,satNames)
if strcmp(name,'Munich'), obj=gs1;
elseif strcmp(name,'NewYork'), obj=gs2;
else, idx=find(strcmp(satNames,name),1); obj=satObjs(idx);
end
end