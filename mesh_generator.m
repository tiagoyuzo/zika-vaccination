clear
clc

% loading boundary coordinates
load rndata2.dat
lla = rndata2;
lla(:,3) = 0; % z components
p=lla2ecef(lla); % transforming coordinate system

% projecting to a plane
u=p(2,:)-p(1,:);
v=p(3,:)-p(2,:);
normal=cross(u,v);
uv=null(normal(:).');
p = p*uv; % projected coordinates

% reescaling
mpx = min(p(:,1));
mpy = min(p(:,2));
a=p(:,1)-mpx;
b=p(:,2)-mpy;
amax = max(a);
bmax = max(b);
% new coordinates
c=[a/amax b/amax];

% identifying cities

load rncities.dat
lla = rncities;
lla(:,3) = 0;
p=lla2ecef(lla);

% projecting to a plane
p = p*uv; % projected coordinates

% reescaling
a=p(:,1)-mpx;
b=p(:,2)-mpy;
% new coordinates
cc=[a/amax b/amax];
m = length(cc);

% generating mesh

n = length(c);
cont = [(1:(n-1))' (2:n)'; n 1];
gd = [2;n;c(cont(:,1),1);c(cont(:,1),2)];
dl = decsg(gd);
% dlmwrite('meshdl.dat',dl,'delimiter','\t')
% return

model=createpde;
geometryFromEdges(model,dl);
% generateMesh(model);
generateMesh(model,'Hmax',0.05)


coord=model.Mesh.Nodes';
conect=model.Mesh.Elements;
nn=length(coord);
nel = length(conect);

figure
hold on
pdeplot(model); axis equal
plot(cc(:,1),cc(:,2),'ok','linewidth',2)

% check cities index
cnames={'Apodi'
    'Areia Branca'
    'Caico'
    'Canguaretama'
    'Ceara Mirim'
    'Joao Camara'
    'Mossoro'
    'Natal'
    'Parnamirim'
    'Rio do fogo'
    'Santa Cruz'
    'Santo Antonio'
    'Sao Goncalo do Amarante'};
for k=1:14
text(cc(k,1),cc(k,2),cnames(k))%num2str(k))
end

% identifying cities

for i=1:nn
    for j=1:m
        dist(i,j)=(coord(i,1)-cc(j,1))^2+(coord(i,2)-cc(j,2))^2;
    end
end

for i=1:m 
    [a1(i),b1(i)]=min(dist(:,i));
end
indcity=b1';

% for i=1:m
%     plot(coord(b1(i),1),coord(b1(i),2),'o','linewidth',2);
%     pause
% end

% finding elements index for cities
elem = zeros(14,7);

[i,j]=find(conect==b1(1)); i1=[i j];
elem(1,1:length(j)) = j;
[i,j]=find(conect==b1(2)); i2=[i j];
elem(2,1:length(j)) = j;
[i,j]=find(conect==b1(3)); i3=[i j];
elem(3,1:length(j)) = j;
[i,j]=find(conect==b1(4)); i4=[i j];
elem(4,1:length(j)) = j;
[i,j]=find(conect==b1(5)); i5=[i j];
elem(5,1:length(j)) = j;
[i,j]=find(conect==b1(6)); i6=[i j];
elem(6,1:length(j)) = j;
[i,j]=find(conect==b1(7)); i7=[i j];
elem(7,1:length(j)) = j;
[i,j]=find(conect==b1(8)); i8=[i j];
elem(8,1:length(j)) = j;
[i,j]=find(conect==b1(9)); i9=[i j];
elem(9,1:length(j)) = j;
[i,j]=find(conect==b1(10)); i10=[i j];
elem(10,1:length(j)) = j;
[i,j]=find(conect==b1(11)); i11=[i j];
elem(11,1:length(j)) = j;
[i,j]=find(conect==b1(12)); i12=[i j];
elem(12,1:length(j)) = j;
[i,j]=find(conect==b1(13)); i13=[i j];
elem(13,1:length(j)) = j;
[i,j]=find(conect==b1(14)); i14=[i j];
elem(14,1:length(j)) = j;

% k=14;
% for j=1:length(i1)
% close(figure(1))
% figure(1)
% fig1 = figure(1);set(fig1, 'Position', [100 100 1200 900])
% hold on
% pdeplot(model); axis equal
% plot(cc(k,1),cc(k,2),'x','linewidth',2)
% for i=1:3
%     x=coord(indcity(k),1);
%     y=coord(indcity(k),2);
%     plot(x,y,'*','linewidth',2)
%     x=coord(conect(i,i14(j,2)),1);
%     y=coord(conect(i,i14(j,2)),2);
%     plot(x,y,'o','linewidth',2)
% end
% j
% pause
% end
