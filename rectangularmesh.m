function [nodecoo, nodenum]=rectangularmesh(nx, ny,lx,ly)
for i= 1:ny
    for j=1:nx
        nodenum(j+((i-1)*nx),:)=[j j+1 j+nx+2 j+nx+1]+(((i-1)*nx)+i-1)*ones(1,4);
    end
end
for i=1:ny+1
    for j=1:nx+1
        nodecoo(j+((i-1)*(nx+1)),:)=[(j-1)*lx/nx (i-1)*ly/ny];
    end
end