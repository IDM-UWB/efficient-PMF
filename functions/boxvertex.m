function boxVertex = boxvertex(n,bound)


% Generate all vertices of box [-1,1]^n
bound = flipud(bound);
for i = 1:n
    if i==1
        boxVertex = [-bound(i) bound(i)];
    else
        nbox = size(boxVertex,2);
        boxVertex = [-bound(i)*ones(1,nbox) bound(i)*ones(1,nbox); boxVertex boxVertex];
    end
end


end