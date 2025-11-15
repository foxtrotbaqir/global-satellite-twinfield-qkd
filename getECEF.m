
function pos=getECEF(obj,t)
if isa(obj,'matlabshared.satellitescenario.GroundStation')
    pos=myLla2ecef([obj.Latitude,obj.Longitude,0]);
else
    pos=states(obj,t,'CoordinateFrame','ecef');
end
end