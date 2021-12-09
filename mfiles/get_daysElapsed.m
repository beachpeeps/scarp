function daysElapsed = get_daysElapsed(a)
[m,n] = size(a);
if m>n
    a = a';
end
aL = length(a);
b = 1:aL;
ind = find(a==1);

if sum(ind) == 0
    daysElapsed = b;
    return
end

inds = [ind(1:end-1); ind(2:end)];
indL = length(inds);

if a(end) == 0
    inds(:,indL+1) = [ind(end);aL];
end

if a(1) == 0
    inds = cat(2,[1;inds(1,1)],inds);
end
indL = length(inds);

aa = [];

for i=1:indL
    aa = [aa (inds(1,i):inds(2,i)-1)-inds(1,i)];
end
if a(end) ==0
    aa = [aa aa(end)+1];
elseif a(end) == 1
    aa = [aa 0];
end





daysElapsed = aa;
