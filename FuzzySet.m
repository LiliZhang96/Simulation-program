function phi=fuzzySet(h,center)
for i=1:length(h)
    for j=1:length(center)
        A(i,j)=exp(-0.5*(h(i)-center(j))'*(h(i)-center(j)));
    end
end
M=0;
for j=1:length(center)
    P(j)=1;
    for i=1:length(h)
        P(j)=P(j)*A(i,j);
    end
    M=M+P(j);
end
phi=[];
for j=1:length(center)
    phi=[phi;P(j)/M];
end

