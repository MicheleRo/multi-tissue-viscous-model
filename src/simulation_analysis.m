clear
clc
value=[];
aa =0;
bb =300*10^(-6);
cc=4500*10^(-6);
filename='X';
 
B=load(strcat(filename,'allparam.mat'), 'N','g1','g2','beta1','beta2','Tmax','x_centers','y_centers', 'mu','aireM','e','injpsm','injnt','nmax');
C=load(strcat(filename,'variables.mat'),'v1_xav','v2_xav','v1_yc','v2_yc','F1_xav','F2_xav','F1_yc','F2_yc','N1','N2','elongtip');
D=load(strcat(filename,'overtime.mat'),'overtime');
E=load(strcat(filename,'TIME'));
delta_x=B.x_centers(2)-B.x_centers(1);
delta_y=B.y_centers(2)-B.y_centers(1);
[X,Y] = meshgrid(B.x_centers,B.y_centers);
em=1;
v1normx=[];
v1normy=[];
v2normx=[];
v2normy=[];


%% plot pressure profile at final time
P1(:)=B.e.*(C.N1(:,end)+C.N2(:,end))./...
    (B.nmax-(C.N1(:,end)+C.N2(:,end)));

figure()
surf(X,Y,reshape(P1(:,end),[B.N,B.N]));
shading interp
colorbar
hold off
view(0,90)

%% plot density profile at final time
figure()
surf(X,Y,reshape(C.N2(:,end),[B.N,B.N]));
hold on
surf(X,Y,reshape(C.N1(:,end),[B.N,B.N]));
shading interp
colorbar
hold off
view(0,90)

%% normalize the velocity vectors for better visual
for i=1:size(C.v1_yc,1)
    nom=sqrt(C.v1_yc(i,end)^2+C.v1_xav(i,end)^2);
    v1normy(i,1)=C.v1_yc(i,end)/nom;
    v1normx(i,1)=C.v1_xav(i,end)/nom;
end
for i=1:size(C.v1_yc,1)
    if (C.N1(i,end-em)<=10^(12))
        v1normx(i,1)=0;
        v1normy(i,1)=0;
    end
end

for i=1:size(C.v2_yc,1)
    nom=sqrt(C.v2_yc(i,end-em)^2+C.v2_xav(i,end-em)^2);
    v2normy(i,1)=C.v2_yc(i,end-em)/nom;
    v2normx(i,1)=C.v2_xav(i,end-em)/nom;
end
for i=1:size(C.v2_yc,1)
    if (C.N2(i,end-em)<=10^(12))
        v2normy(i,end-em)=0;
        v2normx(i,end-em)=0;
    end
end
%% restrain velocities only in area occupied by tissue
for i=1:size(C.v2_yc,1)
    if (C.N2(i,end-em)<=10^(12))
        C.v2_yc(i,end-em)=0;
        C.v2_xav(i,end-em)=0;
    end
end
for i=1:size(C.v1_yc,1)
    if (C.N1(i,end-em)<=10^(12))
        C.v1_yc(i,end-em)=0;
        C.v1_xav(i,end-em)=0;
    end
end

%% plot velocity vectors
figure()
quiver(X,Y,reshape(C.v1_yc(:,end-em),[B.N, B.N]),reshape(C.v1_xav(:,end-em),[B.N, B.N]),'r')
%use the curl function to visualize the vortices in the tissue
axis([aa cc aa bb])
figure()
quiver(X,Y,reshape(C.v2_yc(:,end-em),[B.N, B.N]),reshape(C.v2_xav(:,end-em),[B.N, B.N]),'g')
axis([aa cc aa bb])
%use the div function to plot the divergence (expansion/compression) in the tissue


%% compute NT & PSM width

last_one = find(C.elongtip == 0, 1);
last=floor((C.elongtip-B.N)/B.N)+1;
mid_nt_width=find(C.N2(floor(last(end)/2)*B.N+1:B.N*(floor(last(end)/2)+1),end)>=10^(10));
mid_nt_width=numel(mid_nt_width)*delta_y*10^6;
mid_psm_width=find(C.N1(floor(last(end)/2)*B.N+1:B.N*(floor(last(end)/2)+1),end)>=10^(10));
mid_psm_width=floor(numel(mid_psm_width)/2 )*delta_y*10^6;
if last(end)==0
    ll=find(last == 0, 1)-1;
mid_psm_width=find(C.N1(floor((last(ll))/2)*B.N+1:B.N*(floor((last(ll))/2)+1),last_one-1)>=10^(10));
mid_psm_width=floor(numel(mid_psm_width)/2 )*delta_y*10^6;
end

%% compute elongation rate

if numel(last_one)==0
elongation_rate=60*(delta_x*10^6/(120*6))*(floor((C.elongtip(end)-C.elongtip(60))/B.N)+1);
else
elongation_rate=60*(delta_x*10^6/(120*6))*(floor((C.elongtip(last_one-1)-C.elongtip(60))/B.N)+1);
end


%% compute sliding 
b=5;
a =0;
bb =300*10^(-6);
c=4500*10^(-6);

%compute sliding by averaging on 3 points taken in the middle of the tissue
%(formla in supplementary materials and methods)

pointsinpsm=[3 4 5];
pointsinnt=[10 11 12];

[X,Y] = meshgrid(B.x_centers,B.y_centers);
last_relative_case = find(C.elongtip == 0, 1);
last_relative=floor( (C.elongtip-B.N)/B.N);
firstZeroIndex=floor( (C.elongtip(last_relative_case-1)-B.N)/B.N);
num_valuesonAP=5;

%compute on equally spaced points along the antero-posterior axis
if numel(firstZeroIndex)==0
    time_values=size(C.elongtip,1)-10:size(C.elongtip,1);
    indices= round(linspace(1, floor(last_relative(end)*0.9), num_valuesonAP));
else
    time_values=min(size(C.elongtip,1),last_relative_case)-10:min(size(C.elongtip,1),last_relative_case);
    indices= round(linspace(1, floor(firstZeroIndex*0.9), num_valuesonAP));
end

%average velocity computation
vypsmtraj=zeros(numel(indices),numel(time_values));
vxpsmtraj=zeros(numel(indices),numel(time_values));
vynttraj=zeros(numel(indices),numel(time_values));
vxnttraj=zeros(numel(indices),numel(time_values));

for koo=1:numel(time_values)
    k=time_values(koo);
    ix=0;
    for g=1:numel(indices)
        i=indices(g);
        snpsm=0;
        ix=ix+1;
        for j=pointsinpsm(1,1):pointsinpsm(1,size(pointsinpsm,2))
            vypsmtraj(ix,koo)=vypsmtraj(ix,koo)+C.N1(j+(i-1)*B.N,k)*C.v1_yc(j+(i-1)*B.N,k);
            vxpsmtraj(ix,koo)=vxpsmtraj(ix,koo)+C.N1(j+(i-1)*B.N,k)*C.v1_xav(j+(i-1)*B.N,k);
            snpsm=snpsm+C.N1(j+(i-1)*B.N,k);

        end
        vypsmtraj(ix,koo)=vypsmtraj(ix,koo)/snpsm;
        vxpsmtraj(ix,koo)=vxpsmtraj(ix,koo)/snpsm;
    end
end
for koo=1:numel(time_values)
    k=time_values(koo);

    ix=0;
    for g=1:numel(indices)
        i=indices(g);
        snnt=0;
        ix=ix+1;
        for j=pointsinnt(1,1):pointsinnt(1,size(pointsinnt,2))
            vynttraj(ix,koo)=vynttraj(ix,koo)+C.N2(j+(i-1)*B.N,k)*C.v2_yc(j+(i-1)*B.N,k);
            vxnttraj(ix,koo)=vxnttraj(ix,koo)+C.N2(j+(i-1)*B.N,k)*C.v2_xav(j+(i-1)*B.N,k);
            snnt=snnt+C.N2(j+(i-1)*B.N,k);
        end
        vynttraj(ix,koo)=vynttraj(ix,koo)/snnt;
        vxnttraj(ix,koo)=vxnttraj(ix,koo)/snnt;
    end
end
vpsm=mean(vypsmtraj,2);
vnt=mean(vynttraj,2);

%sliding computation from average tissue velocities
sliding=vnt-vpsm;

figure()
plot(3600*10^(6)*(sliding),'b','LineWidth',2)
