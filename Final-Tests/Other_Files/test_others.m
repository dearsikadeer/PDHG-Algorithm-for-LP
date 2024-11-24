% test non-random problems

test_imgs
test_hard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_imgs(img)
% test on the inpainting problems
%
% change the following path to suit your case
if nargin < 1, img = 2; end
eval(['load data' int2str(img) '_inpaint']);

% display information
fsnr = @(Sig,Noi) 20*log10(norm(Sig(:))/norm(Noi(:)));
SNR1 = fsnr(Img1,Img1-Img0);
fprintf('\n--------------------------------------------\n');
fprintf('  Image inpainting. Initial SNR = %g\n',SNR1)
fprintf('--------------------------------------------\n\n');

xstr(1,:) = 'Recovery by linprog';
xstr(2,:) = 'Recovery by yz-pdhg';
xstr(3,:) = 'Recovery by my-pdhg';
%close all

showoff = 0;
nSolvers = 1;
if exist('yz_pdhg','file'), nSolvers = nSolvers + 1; end
if exist('my_pdhg','file'), nSolvers = nSolvers + 1; end

for k = 1:nSolvers % different solvers
    
    % recover image
    tic; 
    fprintf([xstr(k,:) ' ...\n']);
    Img2 = Inpaint(Img1,Omega,k);
    SNR2 = fsnr(Img2,Img2-Img0);
    fprintf('Recovered SNR = %g\n',SNR2);
    toc
    strL = 'Left: Original. \t';
    strM = 'Middle: Contaminated (SNR: %.2f). \t';
    strR = 'Right: Recovered (SNR: %.2f)';
    str = sprintf([strL strM strR],SNR1,SNR2);
    
    % display images
    if ~showoff, fprintf('\n'), continue; end
    warning('off','images:initSize:adjustingMag');
    figure(k); imshow([Img0 Img1 Img2],[])
    h1 = title(str); h2 = xlabel(xstr(k,:));
    set(h1,'fontsize',12,'FontWeight','normal')
    set(h2,'fontsize',14,'FontWeight','normal')
    drawnow; shg
    
end

end

%% call LP solver to do recovery
function X = Inpaint(X1,Omega,k)
% construct LP data
m = numel(Omega);
n = size(X1,1);
n2 = n^2;

x = X1(:);
y = x(Omega);

e = ones(n,1); I = speye(n);
D = spdiags([e -e],0:1,n-1,n);
D = [kron(I,D); kron(D,I)];
S = speye(n2);
S = S(Omega,:);

% solve LP
% min u + v, st. Dx - u + v = 0, Sx = y, x,u,v >= 0
nD = size(D,1);
A = [D -speye(nD) speye(nD);
    S   sparse(m,2*nD)];
b = [zeros(nD,1); y];
c = [zeros(n2,1); ones(2*nD,1)];
fprintf('m = %i. n = %i.\n',size(A))
tol = 1e-2; maxit = 5000; prt = 0;

switch k
    case 1
        options = optimoptions(@linprog,'display','none',...
            'OptimalityTolerance',tol,...
            'Algorithm','interior-point');
        [x,~,~,output] = ...
            linprog(c,[],[],A,b,zeros(n2+2*nD,1),[],options);
        iter = output.iterations;
        fprintf('number of iterations: %i\n',iter)
    case 2
        if exist('yz_pdhg','file')
            [x,~,iter] = yz_pdhg(A,b,c,tol,maxit,prt);
            fprintf('number of iterations: %i\n',iter)
        end
    case 3
        if exist('my_pdhg','file')
            [x,~,iter] = my_pdhg(A,b,c,tol,maxit,prt);
            fprintf('number of iterations: %i\n',iter)
        end
end

X = reshape(x(1:n2),n,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_hard
% test 2 netlib problems

LPs = {'israel_min'};
solvers = {'yz_pdhg','my_pdhg'};
run = 1:2;

for i = 1:1
    %% load problem
    eval(['load ' LPs{i}])
    [m,n] = size(A);
    
    fprintf('\n------------------------------------\n');
    fprintf(['   ' LPs{i} ': (m,n) = (%i,%i)\n'],m,n)
    fprintf('------------------------------------\n\n');
    %% run linprog
    fprintf('  --- Matlab linprog ---\n');
    LP_options = optimoptions(@linprog,...
        'Display','off',...
        'OptimalityTolerance',1e-8,...
        'Algorithm','interior-point');
    t0 = tic;
    [x,~,~,output,lambda] = linprog...
        (c,[],[],A,b,zeros(n,1),[],LP_options);
    toc(t0)
    y = -lambda.eqlin;
    iter = output.iterations;
    rp = norm(A*x-b)/norm(b);
    rd = norm(min(0,c-A'*y))/norm(c);
    rc = abs(c'*x-b'*y)/max(1e-8,abs(b'*y));
    fprintf('P_res: %14.8e\n',rp);
    fprintf('D_res: %14.8e\n',rd);
    fprintf('pdGap: %14.8e\n',rc);
    fprintf('P_obj: %14.8e\n',c'*x);
    fprintf('Number of iter: %i\n\n',iter)
    %xs = x; ys = y;
    %% run FOM
    for j = run
        if j < 1 || ~exist(solvers{j},'file'), continue, end
        disp(['  --- ' solvers{j} ' ---'])
        t0 = tic;
        tol = 1e-6; maxit = 5e5;  %%%%% call solver %%%%%%
        [x,y,iter] = feval(solvers{j},A,b,c,tol,maxit,0); 
        toc(t0);
        rp = norm(A*x-b)/norm(b);
        rd = norm(min(0,c-A'*y))/norm(c);
        rc = abs(c'*x-b'*y)/max(1e-8,abs(b'*y));
        fprintf('P_res: %14.8e\n',rp);
        fprintf('D_res: %14.8e\n',rd);
        fprintf('pdGap: %14.8e\n',rc);
        fprintf('P_obj: %14.8e\n',c'*x);
        fprintf('Number of iter: %i\n\n',iter)
    end
end
end