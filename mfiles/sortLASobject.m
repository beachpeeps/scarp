function [T, x, y, z, amp, sortedInd] = sortLASobject(T,c)
% sorts the x,y,z,time,and amplitude based on given variable (here T)

[T, sortedInd] = sort(T);

% get values
x = c.x;
y = c.y;
z = c.z;
intensity = get_intensity(c);

x = x( sortedInd);
y = y( sortedInd);
z = z( sortedInd);
amp = intensity( sortedInd);