function [N1_sol, N2_sol]=multi_tissue_brinkman(N,filename,Tmax,g1,g2,injpsm,injnt,beta1,beta2,mu,e,thresh_pz,alpha,nmax_growth, nmax)
tStart = cputime;
%% in this code we consider the 2D of the same depth as the slice of 10um
%% Numerical scheme: Finite volume n staggered grid

%% parameters
s_iter=1;
time=0;% time in seconds
DT=zeros(Tmax/360+1,1);

%% domain definition
a =0; %start of domain
c =300*10^(-6); %width of domain 300 um
b=4500*10^(-6); %length of domain 4500 um

x_gp = linspace(a,b,N+1);
y_gp = linspace(a,c,N+1);

delta_x = x_gp(2) - x_gp(1);
delta_y = y_gp(2) - y_gp(1);

x_centers = a+delta_x/2:delta_x:b;
y_centers = a+delta_y/2:delta_y:c;

[X,Y] = meshgrid(x_centers,y_centers);

%% initializing vectors to save later
v1_yav=zeros(N*N,1);
v2_yav=zeros(N*N,1);
v1_xav=zeros(N*N,Tmax/360);
v2_xav=zeros(N*N,Tmax/360);
v1_yc=zeros(N*N,Tmax/360);
v2_yc=zeros(N*N,Tmax/360);
F1_xav=zeros(N*N,Tmax/360);
F2_xav=zeros(N*N,Tmax/360);
F1_yav=zeros(N*N,1);
F2_yav=zeros(N*N,1);
F1_yc=zeros(N*N,Tmax/360);
F2_yc=zeros(N*N,Tmax/360);
P1_sol=zeros(N*N,1);
gradp1_x=zeros(N*(N+1),1);
gradp1_y=zeros((N+1)*N,1);
F1_x=zeros(N*(N+1),1);
F2_x=zeros(N*(N+1),1);
F1_y=zeros((N+1)*N,1);
F2_y=zeros((N+1)*N,1);
A_v1x=zeros(N*(N+1),N*(N+1));
A_v2x=zeros(N*(N+1),N*(N+1));
A_v1y=zeros(N*(N+1),N*(N+1));
A_v2y=zeros(N*(N+1),N*(N+1));
N0=zeros(N*N);
overtime=zeros(Tmax/360+1,1);
elongtip=zeros(Tmax/360,1);

%% initializing some variables to save compuation time
aireM=delta_x*delta_y;
cpx=mu*delta_x;
cpy=mu*delta_y;
cax1=beta1/(mu*(delta_x)^2);
cay1=beta1/(mu*(delta_y)^2);
cax2=beta2/(mu*(delta_x)^2);
cay2=beta2/(mu*(delta_y)^2);

%% initial condition on tissue densities 

%NT
for i=1:N
    for j=1:N
        if((j>=floor(N/3)+2 && j<=floor(2*N/3)-1) && X(j,i)<=b/3)
            N0(i,j)=0.017536408/10^(-12);
        else
            N0(i,j)=0;
        end
    end
end

N2=zeros(N*N,Tmax/360+1);
for i=1:N
    for j=1:N
        N2(j+(i-1)*N,1)=N0(i,j);
    end
end
N2_sol(:,1)=N2(:,1);

%PSM

N0=zeros(N*N);
minipsm=0.011035977 /10^(-12);
maxpsm=0.013127577/10^(-12);
coef=(minipsm-maxpsm)/(b/3-b/8)^2;
for i=1:N
    for j=1:N
        if ((j<=floor(N/3)+1 || j>=floor(2*N/3)) && X(j,i)<=b/8)
            N0(i,j)=maxpsm;
        elseif ((j<=floor(N/3)+1 || j>=floor(2*N/3)) && X(j,i)>=b/8 && X(j,i)<=b/3)
            N0(i,j)=coef*(X(j,i)-b/8)^2+maxpsm;
        else
            N0(i,j)=0;
        end
    end
end

N1=zeros(N*N,Tmax/360+1);
for i=1:N
    for j=1:N
        N1(j+(i-1)*N,1)=N0(i,j);
    end
end
N1_sol(:,1)=N1(:,1);

%% PZ tip definition
pz=sum(X(1,:)<=b/3);
k2=13+(pz-1)*N;
xp_nt=k2+N;

%% saving simulation parameters
save(strcat(filename,'allparam.mat'),'maxpsm','minipsm','N','g1','g2','beta1','beta2','Tmax','x_centers','y_centers', 'mu','aireM','e','injpsm','injnt','thresh_pz','nmax','nmax_growth');

%% Velocity (in 2D) matrix for NT and PSM, Dirichlet homogeneous BC for the velocity
for k=1:N*(N+1)
    if (mod(k,N+1)==1)
        A_v1x(k,k)=1;
        A_v2x(k,k)=1;
    elseif (mod(k,N+1)==0)
        A_v1x(k,k)=1;
        A_v2x(k,k)=1;

    elseif (k>=2 && k<=N)
        A_v1x(k,k)=5*(cax1)+1;
        A_v1x(k,k+1)=-(cax1);
        A_v1x(k,k-1)=-(cax1);
        A_v1x(k,k+N+1)=-(cax1);

        A_v2x(k,k)=5*(cax2)+1;
        A_v2x(k,k+1)=-(cax2);
        A_v2x(k,k-1)=-(cax2);
        A_v2x(k,k+N+1)=-(cax2);
    elseif (k>N*N && k<N*(N+1))
        A_v1x(k,k)=5*(cax1)+1;
        A_v1x(k,k+1)=-(cax1);
        A_v1x(k,k-1)=-(cax1);
        A_v1x(k,k-(N+1))=-(cax1);

        A_v2x(k,k)=5*(cax2)+1;
        A_v2x(k,k+1)=-(cax2);
        A_v2x(k,k-1)=-(cax2);
        A_v2x(k,k-(N+1))=-(cax2);
    else
        A_v1x(k,k)=4*(cax1)+1;
        A_v1x(k,k+1)=-(cax1);
        A_v1x(k,k-1)=-(cax1);
        A_v1x(k,k+N+1)=-(cax1);
        A_v1x(k,k-(N+1))=-(cax1);

        A_v2x(k,k)=4*(cax2)+1;
        A_v2x(k,k+1)=-(cax2);
        A_v2x(k,k-1)=-(cax2);
        A_v2x(k,k+N+1)=-(cax2);
        A_v2x(k,k-(N+1))=-(cax2);
    end
end

for k=1:N*(N+1)
    if (mod(k,N+1)==1)
        A_v1y(k,k)=1;
        A_v2y(k,k)=1;
    elseif (mod(k,N+1)==0)
        A_v1y(k,k)=1;
        A_v2y(k,k)=1;
    elseif (k>=2 && k<=N)

        A_v1y(k,k)=5*(cay1)+1;
        A_v1y(k,k+1)=-(cay1);
        A_v1y(k,k-1)=-(cay1);
        A_v1y(k,k+N+1)=-(cay1);

        A_v2y(k,k)=5*(cay2)+1;
        A_v2y(k,k+1)=-(cay2);
        A_v2y(k,k-1)=-(cay2);
        A_v2y(k,k+N+1)=-(cay2);
    elseif (k>N*N && k<N*(N+1))
        A_v1y(k,k)=5*(cay1)+1;
        A_v1y(k,k+1)=-(cay1);
        A_v1y(k,k-1)=-(cay1);
        A_v1y(k,k-(N+1))=-(cay1);

        A_v2y(k,k)=5*(cay2)+1;
        A_v2y(k,k+1)=-(cay2);
        A_v2y(k,k-1)=-(cay2);
        A_v2y(k,k-(N+1))=-(cay2);
        %mailles interieures
    else
        A_v1y(k,k)=4*(cay1)+1;
        A_v1y(k,k+1)=-(cay1);
        A_v1y(k,k-1)=-(cay1);
        A_v1y(k,k+N+1)=-(cay1);
        A_v1y(k,k-(N+1))=-(cay1);

        A_v2y(k,k)=4*(cay2)+1;
        A_v2y(k,k+1)=-(cay2);
        A_v2y(k,k-1)=-(cay2);
        A_v2y(k,k+N+1)=-(cay2);
        A_v2y(k,k-(N+1))=-(cay2);
    end
end
%% sparse matrix & LU decomposition for better performance
%A_v2y=sparse(A_v2y);
%A_v1y=sparse(A_v1y);
%A_v2x=sparse(A_v2x);
%A_v1x=sparse(A_v1x);

%[L2y,U2y,P2y]=lu(A_v2y);
%[L2x,U2x,P2x]=lu(A_v2x);
%[L1y,U1y,P1y]=lu(A_v1y);
%[L1x,U1x,P1x]=lu(A_v1x);

%% predefine some variables to use later: width of NT & PSM; injection rate with right units
largpsm=size(1:floor(N/3)+1,2);
largnt=size(floor(N/3)+2:floor(2*N/3)-1,2);
ant=injnt./(aireM.*largnt);
apsm=injpsm./(aireM.*largpsm);

%% time loop
while(time<=Tmax)

    %% pressure
    P1_sol(:,1)=e.*(N1_sol(:,1)+N2_sol(:,1))./(nmax-(N1_sol(:,1)+N2_sol(:,1)));

    %% logistic growth with strong non linearity (alpha)
    G1(:,1)=max(g1.*(1-(N1_sol(:,1)./nmax_growth).^alpha),0);
    G2(:,1)=max(g2.*(1-(N2_sol(:,1)./nmax_growth).^alpha),0);

    %% pressure gradient in x and y
    for i=1:N
        for j=2:N
            gradp1_x(j+(i-1)*(N+1),1)=(P1_sol(j+(i-1)*N,1)-P1_sol(j-1+(i-1)*N,1))/cpx;
        end
    end

    for i=2:N
        for j=1:N
            gradp1_y(i+(j-1)*(N+1),1)=(P1_sol(j+(i-1)*N,1)-P1_sol(j+(i-2)*N,1))/cpy;
        end
    end

    %% compute right hand side of elliptic equations of the velocities with BC into account (homogeneous Dirichlet)
b1x=zeros(N*(N+1),1);
b2x=zeros(N*(N+1),1);
b1y=zeros(N*(N+1),1);
b2y=zeros(N*(N+1),1);

for k=1:N*(N+1)
    if (mod(k,N+1)==1)
        b1x(k)=0;
        b2x(k)=0;
    elseif mod(k,N+1)==0
        b1x(k)=0;
        b2x(k)=0;
    else
        b1x(k)=-gradp1_x(k);
        b2x(k)=-gradp1_x(k);
    end
end

for k=1:N*(N+1)
    if (mod(k,N+1)==1)
        b1y(k)=0;
        b2y(k)=0;
    elseif mod(k,N+1)==0
        b1y(k)=0;
        b2y(k)=0;
    else
        b1y(k)=-gradp1_y(k);
        b2y(k)=-gradp1_y(k);
    end
end

v1x=A_v1x\b1x;
v1y=A_v1y\b1y;
v2x=A_v2x\b2x;
v2y=A_v2y\b2y;

    %% compute the fluxes in x and y for PSM & NT
    for i=1:N
        for j=2:N
            index=j+(i-1)*(N+1);
            F1_x(index,1)=max(v1x(index,1),0)*N1_sol(j-1+(i-1)*N,1)-max(-v1x(index,1),0)*...
                N1_sol(j+(i-1)*N,1);
            F2_x(index,1)=max(v2x(index,1),0)*N2_sol(j-1+(i-1)*N,1)-max(-v2x(index,1),0)*...
                N2_sol(j+(i-1)*N,1);
        end
    end

    for j=1:N
        for i=2:N
            index=i+(j-1)*(N+1);
            F1_y(index,1)=max(v1y(index,1),0)*N1_sol(j+(i-2)*N,1)-max(-v1y(index,1),0)*...
                N1_sol(j+(i-1)*N,1);
            F2_y(index,1)=max(v2y(index,1),0)*N2_sol(j+(i-2)*N,1)-max(-v2y(index,1),0)*...
                N2_sol(j+(i-1)*N,1);
        end
    end
   
    %% compute beginning of PZ
    psm_entry_left2_right2=zeros(N*N,1);
    nt_entry2=zeros(N*N,1);
    test=(N2_sol(9+(pz-1)*N:13+(pz-1)*N,1)>=thresh_pz);
    if (sum(test)==5)
        k2=13+(pz-1)*N;
        xp_nt=k2+N;
        pz=pz+1;
    end
    if pz>=N
        break
    end
    %% compute zone of cell injection in each tissue
    psm_entry_left2_right2(xp_nt- largnt -largpsm+1:xp_nt- largnt)=1;
    psm_entry_left2_right2(xp_nt+1:xp_nt+ largpsm)=1;
    nt_entry2(xp_nt- largnt+1:xp_nt)=1;

    %% compute CFL conditions (positive densities and smaller than nmax)
    lambdax=4*max([max(max(-v1x(:,1),0)), max(max(-v2x(:,1),0)),...
        max(max(v1x(:,1),0)),max(max(v2x(:,1),0))]);
    lambday=4*max([max(max(-v1y(:,1),0)),max(max(-v2y(:,1),0)),...
        max(max(v1y(:,1),0)),max(max(v2y(:,1),0))]);
    delta_t=min([delta_x^2/lambdax, delta_y^2/lambday, 0.8*(nmax./max([max(N1_sol(:,1)),max(N2_sol(:,1))])-1)/(max(max(G1(:),G2(:)))+(1/delta_x)*lambdax+(1/delta_y)*lambday)]);
    
    %% predefine some variables to use later multiple times
    rat=delta_t/delta_x;
 
    psm_entry_left2_right2=delta_t.*psm_entry_left2_right2.*apsm;
    nt_entry2=nt_entry2.*delta_t.*ant;

    %% save variables each 6 minutes
    if (time+delta_t>=s_iter*360)
        DT(s_iter+1)=delta_t;
        for i=1:N
            for j=1:N
                v1_xav(j+(i-1)*N,s_iter)=(v1x(j+1+(i-1)*(N+1),1)+v1x(j+(i-1)*(N+1),1))/2;
                v2_xav(j+(i-1)*N,s_iter)=(v2x(j+1+(i-1)*(N+1),1)+v2x(j+(i-1)*(N+1),1))/2;
            end
        end

        for j=1:N
            for i=1:N
                v1_yav(i+(j-1)*N,1)=(v1y(i+1+(j-1)*(N+1),1)+v1y(i+(j-1)*(N+1),1))/2;
                v2_yav(i+(j-1)*N,1)=(v2y(i+1+(j-1)*(N+1),1)+v2y(i+(j-1)*(N+1),1))/2;
            end
        end

        for i=1:N
            for j=1:N
                v1_yc(j+(i-1)*N,s_iter)=v1_yav(i+(j-1)*N,1);
                v2_yc(j+(i-1)*N,s_iter)=v2_yav(i+(j-1)*N,1);
            end
        end

        for i=1:N
            for j=1:N
                F1_xav(j+(i-1)*N,s_iter)=(F1_x(j+1+(i-1)*(N+1),1)+F1_x(j+(i-1)*(N+1),1))/2;
                F2_xav(j+(i-1)*N,s_iter)=(F2_x(j+1+(i-1)*(N+1),1)+F2_x(j+(i-1)*(N+1),1))/2;
            end
        end

        for j=1:N
            for i=1:N
                F1_yav(i+(j-1)*N,1)=(F1_y(i+1+(j-1)*(N+1),1)+F1_y(i+(j-1)*(N+1),1))/2;
                F2_yav(i+(j-1)*N,1)=(F2_y(i+1+(j-1)*(N+1),1)+F2_y(i+(j-1)*(N+1),1))/2;
            end
        end

        for i=1:N
            for j=1:N
                F1_yc(j+(i-1)*N,s_iter)=F1_yav(i+(j-1)*N,1);
                F2_yc(j+(i-1)*N,s_iter)=F2_yav(i+(j-1)*N,1);
            end
        end
        N1(:,s_iter+1)=N1_sol(:,1);
        N2(:,s_iter+1)=N2_sol(:,1);
        overtime(s_iter,1)=time;
        titi=strcat(filename,'overtime');
        save(titi,'overtime');
        tim=strcat(filename,'TIME');
        save(tim,'DT')
        o=strcat(filename,'variables');
        elongtip(s_iter)=xp_nt;
        save(o,'v1_xav','v2_xav','v1_yc','v2_yc','F1_xav','F2_xav','F1_yc','F2_yc','N1','N2','elongtip');
        s_iter=s_iter+1;
    end
    %% compute densities
    time=time+delta_t;
    F1_x=rat.*F1_x;
    F1_y=rat.*F1_y;
    F2_x=rat.*F2_x;
    F2_y=rat.*F2_y;

    for i=1:N
        for j=1:N
            index = j + (i - 1) * N;
            N1_sol(index,2)= (1+G1(index,1)*delta_t)*N1_sol(index,1)-...
                (F1_x(j+1+(i-1)*(N+1),1)-F1_x(j+(i-1)*(N+1),1))-(F1_y(i+1+(j-1)*(N+1),1)-F1_y(i+(j-1)*(N+1),1))+psm_entry_left2_right2(index);


            N2_sol(index,2)= (1+G2(index,1)*delta_t)*N2_sol(index,1)-(F2_x(j+1+(i-1)*(N+1),1)-F2_x(j+(i-1)*(N+1),1))-(F2_y(i+1+(j-1)*(N+1),1)-...
                F2_y(i+(j-1)*(N+1),1))+nt_entry2(index);
        end
    end
%% update densities
    N1_sol(:,1)=N1_sol(:,2);
    N2_sol(:,1)=N2_sol(:,2);
end
tEnd = cputime - tStart


