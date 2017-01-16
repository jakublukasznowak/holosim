function ap = apod1 (x,x0,r0,dr)

r=abs(x-x0);

ap=cos(pi/2/dr*(r-(r0-0.5*dr))).^2;
ap(r<r0-0.5*dr)=1;
ap(r>r0+0.5*dr)=0;

end