%% to have the tissue colours (green for NT and red for PSM) download first the multigradient.m from Mathworks
f='multigradient([1 1 1; 0 0.7 0],''pts'', [8 15])';
f2='multigradient([1 1 1; 1 0 0],''pts'', [8 15])';

filename='X'; %enter the filename you gave when running the main simulation
a =0;
b =300*10^(-6);
c=4500*10^(-6);
A=load(strcat(filename,'TIME'));
t=A(1,1);
it=1;
B=load(strcat(filename,'allparam.mat'),'maxpsm','minipsm','N','g1','g2','beta1','beta2','Tmax','x_centers','y_centers', 'mu','aireM','e','injpsm','injnt','thresh_pz','nmax','nmax_growth');
C=load(strcat(filename,'variables.mat'),'v1_xav','v2_xav','v1_yc','v2_yc','F1_xav','F2_xav','F1_yc','F2_yc','N1','N2','elongtip');
D=load(strcat(filename,'overtime.mat'),'overtime');
[X,Y] = meshgrid(B.x_centers,B.y_centers);
delta_x=B.x_centers(2)-B.x_centers(1);
delta_y=B.y_centers(2)-B.y_centers(1);
aireM=delta_x*delta_y;

%% create video, give video name
press=VideoWriter('movie_densities_final.avi');
press.Quality=100;
open(press);

%% plot densities 
while(it<=length(A.DT)-1)
    t=A.DT(it);
    figure(2)
    l1=surf(X,Y,reshape(C.N2(:,it),[B.N,B.N]));
    h1=colormap(eval(f));
    freezeColors
    hold on
    h2=colormap(gca, eval(f2));
    colorbar
    shading interp
    hold on
    l2=surf(X,Y,reshape(C.N1(:,it),[B.N,B.N]));
    hold on
    colorbar
    caxis([0; nmax])
    shading interp
    view(90,90)
    title(['Time= ',num2str(D.overtime(it)/(3600)),' hours'])
    currFrame1=getframe(gcf);
    writeVideo(press,currFrame1)
    clf
    hold on
    it=it+1;
end
close(press);