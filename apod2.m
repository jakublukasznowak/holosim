function ap = apod2 (x,y,x0,y0,r0,dr)

r=sqrt((x-x0).^2+(y-y0).^2);

ap=cos(pi/2/dr*(r-(r0-0.5*dr))).^2;
ap(r<r0-0.5*dr)=1;
ap(r>r0+0.5*dr)=0;

end