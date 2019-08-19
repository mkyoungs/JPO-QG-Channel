global tau H1 H2 u2meandagger nx u1meandagger kappa psimeandagger psibar1 x yp...
    q1totm psi2meandagger psibar2 q2totm uq1 vq1 uq2 vq2 teapechange fr1...
    r H f0 ny v2mean v1mean


tic
tau0s =  [0 0.0125 0.025 0.05 0.1 0.2 0.3 0.4 0.6 0.8 1.5];
ots = [-4 -2 0 2];
top0s = 1;
kappa = 1.3;

r = 0.03;
counterind = 20;

%tau0s = [0 0.0001 0.0002 0.0004 0.0008 0.0016 0.0032 0.0064];
%ots = -0.01;
%tau0s = 0;

if top0s == 0
    dbname = 'flat';
elseif top0s == 1
    dbname = 'ridge';
elseif top0s == 2
    dbname = 'ridge2';
else
    dbname = '';
end

load(dbname)
%ots = [-4 -2 -1 -0.5 -0.25 -0.1 -0.05 -0.025 -0.01 -0.005 0];

%kappa = 20;
for j = 4:5%length(tau0s)
    for k = 2%1:length(ots)
        close all
        
        tau0=tau0s(j) ;
        H0 = ots(k); % overturning in Sv
        
        %tau0 = 0.0001;
        
        tau=tau0*1000*3600^2*24^2*cos(pi*yp/width); % in km/day units
        %         tau=tau0*1000*3600^2*24^2*cos(pi*yp/width*2);
        %         tau(1:ny/4) = 0;
        %         tau(3*ny/4:end) = 0;% in km/day units
        tau = tau/1000^4/H1;
        
        top0 = top0s*1000;
        %         filename = '~/Desktop/ch2/results/test1.dat';
        %titlepart = ' test';
        %
        %filename = sprintf('~/Desktop/ch2/results/tau%gtop%got%gr%g.dat',tau0,top0s,H0,r);
        filename = sprintf('~/Desktop/ch2/results/tau%gtop%got%g.dat',tau0,top0s,H0);
        titlepart = sprintf(' tau_0 = %g ridge = %g m ot = %g Sv',tau0,top0,H0);
        
        factor1 = 1000/3600/24*4*1000;
        H_0 = -H0*pi/width/factor1;
        H = H_0*sin(pi*yp/width);
        
        
        try
            fid=fopen(filename,'r');
            
            
            topography = f0/H2*top0s*ones(size(yp))*exp(-(x-800).^2/200^2);
            topo = 1000*top0s*ones(size(yp))*exp(-(x-800).^2/200^2);
            %mycolor;
            
            
            fr1=1.0/rd/rd/(1.0+del);fr2=del*fr1;
            
            
            ts=[];
            q1m=[];
            q2m=[];
            q1tot = [];
            q2tot = [];
            u1m=[];
            u2m=[];
            qy1m=[];
            qy2m=[];
            sm=[];
            smx=[];
            phim=[];
            psim=[];
            psidagger=[];
            psi2dagger =[];
            v2dagger=[];
            v1dagger=[];
            u2dagger=[];
            u1dagger=[];
            eddyv=[];
            eddyv2=[];
            hm=[];
            psibt = [];
            psibc = [];
            vpqp1 = [];
            vpqp2 = [];
            psibar1mean = [];
            psibar2mean = [];
            psi1totmean = [];
            psi2totmean = [];
            ape = [];
            u1sq = [];
            u2sq = [];
            v1sq = [];
            v2sq = [];
            
            
            counter = 0;
            psimeandagger = [];
            psi2meandagger = [];
            hmmean = [];
            eddyvmean = [];
            eddyv2mean = [];
            u1mean = [];
            u2mean = [];
            u1meandagger = [];
            u2meandagger = [];
            v1mean = [];
            v2mean = [];
            phimean = [];
            psimean = [];
            bctrans = [];
            bttrans = [];
            pv1flux = [];
            pv2flux = [];
            q1mean = [];
            q2mean = [];
            q1sq = [];
            q2sq = [];
            q1totm = [];
            q2totm = [];
            teapechange = [];
            psibar1m = [];
            psibar2m = [];
            uq1 = [];
            vq1 = [];
            uq2 = [];
            vq2 = [];
            reynolds1 = [];
            reynolds2 = [];
            h2 = [];
            uv = [];
            uh = [];
            vh = [];
            totalape =[];
                                psi1sq = [];
                    psi2sq = [];
                    psi1psi2 = [];
                    meanpsi1tot = [];
                    meanpsi2tot = [];
            
            
            while(1)
                [t,nn]=fread(fid,1,'float');
                counter = counter+1;
                
                if nn ~= 1
                    break;
                end
                psi1=fread(fid,[nx,nyp],'float')';
                psi2=fread(fid,[nx,nyp],'float')';
                q1=fread(fid,[nx,nyp],'float')';
                q2=fread(fid,[nx,nyp],'float')';
                u1=fread(fid,[nx,nyp],'float')';
                u2=fread(fid,[nx,nyp],'float')';
                v1=fread(fid,[nx,nyp],'float')';
                v2=fread(fid,[nx,nyp],'float')';
                ubar1=fread(fid,[1,nyp],'float')';
                ubar2=fread(fid,[1,nyp],'float')';
                qbary1=fread(fid,[1,nyp],'float');
                qbary2=fread(fid,[1,nyp],'float');
                phi=fread(fid,[1,nyp],'float')';
                qbar1=cumsum(qbary1)*dx-qbary1(1)*dx/2-qbary1(nyp)/2*dx;
                qbar2=cumsum(qbary2)*dx-qbary2(1)*dx/2-qbary2(nyp)/2*dx;
                qtot1=q1+qbar1'*ones(1,nx);
                qtot2=q2+qbar2'*ones(1,nx)+topography;
                psibar1=-cumsum(ubar1)*dx+ubar1(1)*dx/2+ubar1(nyp)/2*dx;
                psibar2=-cumsum(ubar2)*dx+ubar2(1)*dx/2+ubar2(nyp)/2*dx;
                psitot1=psibar1*ones(1,nx)+psi1;
                psitot2=psibar2*ones(1,nx)+psi2;
                
                %  purge_tmp_files;
                %  title(sprintf('t = %6.2f',t));
                %  contourg(x,yp,qtot1);
                dv=[qtot2;zeros(5,nx);qtot1];
                %  dv=psitot1-psitot2;
                %  dv=[psitot2;zeros(5,nx);psitot1];
                %  dv=[qtot1];
                %  df(t,dv,min(min(dv)),max(max(dv)),4);
                %  df(t,dv,min(min(dv)),max(max(dv)),1);
                imagesc(dv);axis('xy','equal');title(num2str(t));drawnow
                %pause;
                ts=[ts,t];
                hm=[hm,(psibar1-psibar2)*fr2/f0*H2];
                eddyv=[eddyv,mean(v1.*(psi1-psi2)*fr1/f0,2)];
                eddyv2=[eddyv2,mean(v2.*(psi2-psi1)*fr2/f0,2)];
                %q1m=[q1m,[min(min(qtot1));max(max(qtot1))]];
                %q2m=[q2m,[min(min(qtot2)),max(max(qtot2))]];
                u1m=[u1m,ubar1];
                u2m=[u2m,ubar2];
                qy1m=[qy1m,qbary1'];
                qy2m=[qy2m,qbary2'];
                q1m = cat(3,q1m,q1);
                q2m = cat(3,q2m,q2);
                q1tot = cat(3,q1tot,qtot1);
                q2tot = cat(3,q2tot,qtot2);
                phim = [phim,phi];
                psim = cat(3,psim,(psitot1*fr2+psitot2*fr1)/(fr1+fr2));
                %   psibt = cat(3,psibt,(fr2*psi1+fr1*psi2)/(fr1+fr2));
                %   psibc = cat(3,psibc,1/2*(psi1-psi2));
                psidagger = cat(3,psidagger,psi1);
                psi2dagger= cat(3,psi2dagger,psi2);
                v1dagger = cat(3,v1dagger,v1);
                v2dagger = cat(3,v2dagger,v2);
                u1dagger = cat(3,u1dagger,u1+ubar1*ones(1,nx));
                u2dagger = cat(3,u2dagger,u2+ubar2*ones(1,nx));
                psibar1mean = [psibar1mean, psibar1];
                psibar2mean = [psibar2mean, psibar2];
                psi1totmean = cat(3,psi1totmean, psitot1);
                psi2totmean = cat(3,psi2totmean, psitot2);

                ape = cat(3,ape,(psitot1-psitot2).^2);
                
                [t,max(ubar1),max(max(q1))];
                
                if counter == counterind
                    counter = 0;
                    psimeandagger = cat(3,psimeandagger,mean(psidagger,3));
                    psi2meandagger = cat(3,psi2meandagger,mean(psi2dagger,3));
                    hmmean = [hmmean, mean(hm,2)];
                    eddyvmean = [eddyvmean, mean(eddyv,2)];
                    eddyv2mean = [eddyv2mean, mean(eddyv2,2)];
                    u1mean = [u1mean, mean(u1m,2)];
                    u2mean = [u2mean, mean(u2m,2)];
                    phimean = [phimean, mean(phim,2)];
                    psimean = cat(3,psimean,mean(psim,3));
                    v1mean = cat(3,v1mean,mean(v1dagger,3));
                    v2mean = cat(3,v2mean,mean(v2dagger,3));
                    u1meandagger = cat(3,u1meandagger,mean(u1dagger,3));
                    u2meandagger = cat(3,u2meandagger,mean(u2dagger,3));
                    u1sq = cat(3,u1sq,mean(u1dagger.^2,3));
                    u2sq = cat(3,u2sq,mean(u2dagger.^2,3));
                    v1sq = cat(3,v1sq,mean(v1dagger.^2,3));
                    v2sq = cat(3,v2sq,mean(v2dagger.^2,3));
                    
                    bctrans = [bctrans,(sum(mean(u1m,2)) - sum(mean(u2m,2)))*1000/3600/24*4*dx]; % right now its u1-u2
                    bttrans = [bttrans,(fr2*sum(mean(u1m,2))+fr1*sum(mean(u2m,2)))/(fr1+fr2)*1000/3600/24*4*dx];
                    pv1flux = cat(3,pv1flux,mean(v1dagger.*diffxy(yp,q1m,1),3)+mean(u1dagger.*diffxy(x,q1m,2),3));
                    pv2flux = cat(3,pv2flux,mean(v2dagger.*diffxy(yp,q2m,1),3)+mean(u2dagger.*diffxy(x,q2m,2),3));
                    q1mean = cat(3,q1mean,mean(q1m,3));
                    q2mean = cat(3,q2mean,mean(q2m,3));
                    q1sq = cat(3,q1sq,mean(q1m.^2,3));
                    q2sq = cat(3,q2sq,mean(q2m.^2,3));
                    q1totm = cat(3,q1totm,mean(q1tot,3));
                    q2totm = cat(3,q2totm,mean(q2tot,3));
                    teapechange = cat(3,teapechange, mean(Jacobian(x,yp,psi1totmean,psi2totmean),3));
                    psibar1m = [psibar1m, mean(psibar1mean,2)];
                    psibar2m = [psibar2m, mean(psibar2mean,2)];
                    uq1 = cat(3,uq1,mean((u1dagger).*q1tot,3));
                    vq1 = cat(3,vq1,mean((v1dagger).*q1tot,3));
                    uq2 = cat(3,uq2,mean((u2dagger).*q2tot,3));
                    vq2 = cat(3,vq2,mean((v2dagger).*q2tot,3));
                    reynolds1 = [reynolds1, mean(mean(v1.*diffxy(yp,u1),2),3)];
                    reynolds2 = [reynolds2, mean(mean(v2.*diffxy(yp,u2),2),3)];
                    uv= cat(3,uv,mean(H1*(u1dagger.^2+v1dagger.^2) + H2*(u2dagger.^2+v2dagger.^2),3));
                    h2 = cat(3,h2,mean((psidagger-psi2dagger).^2,3));
                    uh = cat(3,uh,mean(u1dagger.*(psidagger-psi2dagger),3));
                    vh = cat(3,vh,mean(v1dagger.*(psidagger-psi2dagger),3));
                    totalape = cat(3,totalape,mean(ape,3));
                    psi1sq = cat(3,psi1sq,mean(psi1totmean.^2,3));
                    psi2sq = cat(3,psi2sq,mean(psi2totmean.^2,3));
                    psi1psi2 = cat(3,psi1psi2,mean(psi1totmean.*psi2totmean,3));
                    meanpsi1tot = cat(3,meanpsi1tot,mean(psi1totmean,3));
                    meanpsi2tot = cat(3,meanpsi2tot,mean(psi2totmean,3));

                    q1m=[];
                    q2m=[];
                    u1m=[];
                    u2m=[];
                    qy1m=[];
                    q2tot=[];
                    q1tot=[];
                    qy2m=[];
                    sm=[];
                    smx=[];
                    phim=[];
                    psim=[];
                    psidagger=[];
                    psi2dagger =[];
                    v2dagger=[];
                    v1dagger=[];
                    u2dagger=[];
                    u1dagger=[];
                    eddyv=[];
                    eddyv2=[];
                    hm=[];
                    psibar1mean = [];
                    psibar2mean = [];
                    psi1totmean = [];
                    psi2totmean = [];
                    ape = [];
                    
                end
                
                
                
            end
            
            clear q1m q2m u1m u2m qy1m qy2m phim psim psidagger psi2dagger
            clear v2dagger v1dagger u2dagger u1dagger eddyv eddyv2 
            clear psibar1mean psibar2mean
            disp(ts(end))
            %pclose(fid);
            fclose(fid);
            
            factor = 1000/3600/24;
            
            lim = min(300,floor(length(ts)/counterind));
            
            psimeandagger = mean(psimeandagger(:,:,end-lim+1:end),3);
            psi2meandagger =  mean(psi2meandagger(:,:,end-lim+1:end),3);
            hmmean = mean(hmmean(:,end-lim+1:end),2);
            eddyvmean = mean(eddyvmean(:,end-lim+1:end),2);
            eddyv2mean = mean(eddyv2mean(:,end-lim+1:end),2);
            u1mean = mean(u1mean(:,end-lim+1:end),2);
            u2mean = mean(u2mean(:,end-lim+1:end),2);
            u1meandagger = mean(u1meandagger(:,:,end-lim+1:end),3);
            u2meandagger = mean(u2meandagger(:,:,end-lim+1:end),3);
            phimean =  mean(phimean(:,end-lim+1:end),2);
            psimean =  mean(psimean(:,:,end-lim+1:end),3);
            v1mean =  mean(v1mean(:,:,end-lim+1:end),3);
            v2mean =  mean(v2mean(:,:,end-lim+1:end),3);
            q1mean = mean(q1mean(:,:,end-lim+1:end),3);
            q2mean =  mean(q2mean(:,:,end-lim+1:end),3);
            pv1flux = mean(pv1flux(:,:,end-lim+1:end),3) - v1mean.*diffxy(yp,q1mean,1) ...
                + u1mean.*diffxy(x,q1mean,2);
            pv2flux = mean(pv2flux(:,:,end-lim+1:end),3)  - v2mean.*diffxy(yp,q2mean,1) ...
                + u2mean.*diffxy(x,q2mean,2);
            q1sq = mean(q1sq(:,:,end-lim+1:end),3)-q1mean.^2;
            q2sq = mean(q2sq(:,:,end-lim+1:end),3)-q2mean.^2;
            q1totm = mean(q1totm(:,:,end-lim+1:end),3);
            q2totm = mean(q2totm(:,:,end-lim+1:end),3);
            teapechange = mean(teapechange(:,:,end-lim+1:end),3);
            psibar1m =  mean(psibar1m(:,end-lim+1:end),2);
            psibar2m =  mean(psibar2m(:,end-lim+1:end),2);
            uq1 = mean(uq1(:,:,end-lim+1:end),3);
            vq1 = mean(vq1(:,:,end-lim+1:end),3);
            uq2 = mean(uq2(:,:,end-lim+1:end),3);
            vq2 = mean(vq2(:,:,end-lim+1:end),3);
            reynolds1 = mean(reynolds1(:,end-lim+1:end),2);
            reynolds2 = mean(reynolds2(:,end-lim+1:end),2);
            uv= mean(uv(:,:,end-lim+1:end),3);
            h2 = mean(h2(:,:,end-lim+1:end),3);
            uh = mean(uh(:,:,end-lim+1:end),3);
            vh = mean(vh(:,:,end-lim+1:end),3);
            totalape = mean(totalape(:,:,end-lim+1:end),3);
            
            eke(k,j) = mean(mean(H1*(mean(u1sq,3)-u1meandagger.^2 + mean(v1sq,3)-v1mean.^2)+...
                    H2*(mean(u2sq,3)-u2meandagger.^2 + mean(v2sq,3)-v2mean.^2)));
            epe(k,j) = mean(mean(fr1*H1*(mean(psi1sq,3)+mean(psi2sq,3) -2*mean(psi1psi2,3)-...
                        mean(meanpsi1tot,3).^2 - mean(meanpsi2tot,3).^2 +...
                        2*(mean(meanpsi1tot,3).*mean(meanpsi2tot,3)))));
            
            %ratio(k,j) = EnergyBudget;
            u1meanave(k,j) = mean(u1mean)
            u2meanave(k,j) = mean(u2mean)
            meanape(k,j) = mean(mean(totalape));
            %vhmean(:,j) = nanmean(vh,1);
            diffpvflux(k,j) = dx*trapz(mean(vq1-vq2,2))*10^6/3600/24/3600/24;
   
            %%
            
            plotch
            %mombudget
            %chvort
            
            
        catch errorname
            % set things to NaNs
            disp('no file')
            baroclinictransport(k,j) = NaN;
            barocliniccenter(k,j) = NaN;
            barotropictransport(k,j) = NaN;
            residualottot(k,j) = NaN;
%             formdragcoeff(k,j) = NaN;
            dragcoeff(k,j) = NaN;
            sehf(k,j) = NaN;
            tehf(k,j) = NaN;
            alpha(k,j) = NaN;
            kappamean(k,j) = NaN;
            ratio(k,j) = NaN;
%             vstarcoeff(k,j) = NaN;
%             teddycoeff(k,j) = NaN;
%             seddycoeff(k,j) = NaN;
            psimean1coeff(k,j) =  NaN;
            psimean2coeff(k,j) = NaN;
%             formbccoeff(k,j) = NaN;
%             formbtcoeff(k,j) = NaN;
%             formbtgyre(k,j) = NaN;
%             formbtjet(k,j) = NaN;
            seddygyre(k,j) = NaN;
            seddyjet(k,j) = NaN;
            psibtmax(k,j) = NaN;
            psibtmin(k,j) = NaN;
            %             totalenergy(k,j) = NaN;
            %             meanwind(k,j) = NaN;
            %             meanapete(k,j) = NaN;
            %             meanapese(k,j) = NaN;
            %             meanheatflux(k,j) = NaN;
            %             meandamping(k,j) = NaN;
            %             meankappa(k,j) = NaN;
        end
        
        
    end
end
toc
% h7 = figure;
% contourf(x,yp,v2meandagger*1000/3600/24)
% ti = strcat('Lower Layer v- Velocity' , titlepart);
% xlabel('X (km)')
% ylabel('Y (km)')
% colorbar
% title(ti)
% SaveFigureThesis(h7,ti)

% %% compute enstrophy budget
%
% % % compute enstrophy
% % figure; contourf(Z1)
% %
% % % mean advection should be zero
% %
% %
% % % eddy advection
% % figure; contourf(DZ1u+DZ1v)
% eddyadvec = DZ1u+DZ1v;
%
%
% % mean pv
% Qyvpqp1 = qy1m.*vpqp1;
% Qyvpqp2 = qy2m.*vpqp2;
%
%
% % kappa damping
% kddyqp1 = kappa*diffxy(yp,Z1,1,2);
% kddyqp2 = kappa*diffxy(yp,Z2,1,2);
% % figure; contourf(kddyqp1)
%
% % bottom drag damping
%
%
%
% figure; plot(yp,mean(Qyvpqp1(:,end-lim+1:end),2),yp,mean(eddyadvec(:,end-lim+1:end),2),...
%     yp,mean(kddyqp1(:,end-lim+1:end),2))
% legend('Mean PV term','Eddy Advection','Kappa','location','best')


            %%


%%


%

save(dbname,'baroclinictransport','barotropictransport','residualottot', 'formdragcoeff', ...
    'dragcoeff', 'sehf', 'tehf', 'tau0s', 'ots', 'top0s', 'barocliniccenter',...
    'vstarcoeff','teddycoeff','seddycoeff','psimean1coeff','psimean2coeff',...
    'formbtcoeff','formbccoeff','seddygyre','seddyjet','psibtmax','psibtmin','meanape','vhmean')

