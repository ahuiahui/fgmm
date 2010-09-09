
clf
nstates = 4;
dim = 2;
dat = [];
for i=1:nstates
  mu = ones(1,dim);
  mu(1,mod(i,dim)+1) = i;
  dat_ = randn(1000,dim)*.1 + repmat(mu,1000,1);
  dat = [dat; dat_];
end

%plot(dat(:,1),dat(:,2),'k+')
hold on

disp('starting k-means')
[Priors, Mu , Sigma] = EM_init_kmeans(dat',nstates)

gmm = struct('Prior',Priors,'Mu',Mu,'Sigma',Sigma,'data',dat', ...
	     'vars',ones(dim,1),'means',zeros(dim,1));

drawGMM(gmm,1:2,[1 0 0]);

disp('end k-mean')
Mu = Mu';
size(dat)
[Priors, Mu , Sigma] = gmm_em(dat',nstates,Priors,Mu,Sigma)

gmm = struct('Prior',Priors,'Mu',Mu','Sigma',Sigma,'data',dat', ...
 	     'vars',ones(dim,1),'means',zeros(dim,1));
drawGMM(gmm);

figure()

[Priors, Mu , Sigma] = gmm_em(dat',nstates)

gmm = struct('Prior',Priors,'Mu',Mu','Sigma',Sigma,'data',dat', ...
  	     'vars',ones(dim,1),'means',zeros(dim,1));
drawGMM(gmm);
