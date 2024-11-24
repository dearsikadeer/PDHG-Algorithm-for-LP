%      Tests for the final project 2024
% -------------------------------------------------
% (1) Copy your code my_pdhg.m into this directory.
% (2) Run this test script with all the codes here.
% (3) Submit the generated file: MyOutput.txt.
% -------------------------------------------------

addpath ./Other_Files/
dfile = 'MyOutput.txt';
if exist(dfile, 'file'); delete(dfile); end
Run = [0 1 0 1]; prt = 0; 
rep = repmat('*',1,20);

%% starting
diary(dfile), diary on
t00 = tic;

%% small sizes
ptype = ' p = 1:2 tol = 1e-14 ';
fprintf(['\n***** ' ptype rep '\n'])
for p = 1:2
    test_rand(p,Run,prt,1e-14);
end

%% mid-large sizes
ptype = ' p = 5:5:25 tol = 1e-6 ';
fprintf(['\n***** ' ptype rep '\n'])
for p = [5 10 15 20 25]
    test_rand(p,Run,prt,1e-6);
end

%% test others
ptype = ' test on image and netlib problem ';
fprintf(['\n***** ' ptype rep '\n'])
test_others

%% finishing
toc(t00)
disp(datetime)
diary off