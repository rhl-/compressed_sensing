# Date created: 12/2/2014
# Author: Victor Minden
# Title: ClusterDriverMatlab.py
# Purpose: Construct many different matlab scripts and submit them to a Stanford cluster using qsub
#	       for David Donoho's compressed sensing project.
# Comment: Don't forget to put Matlab files in necessary path!

import os
import subprocess


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



Ns = [32, 64]
Measurements = ['Entry','Perm','RSPerm','CSPerm','RGPerm','CGPerm','RDirac','CDirac','RGauss','CGauss']
#Mats = ['RPSD','RSYM','HPSD','HERM']
Mats = ['HPSD']


matlabTemplate = '''
cd ~/CSProject/{3}
CVX_PATH = '~/Matlab/cvx';
CS_PATH = '~/Matlab/CS';
pathstr = sprintf('addpath %s -begin', CS_PATH);
eval(pathstr)
pathstr = sprintf('addpath %s -begin',CVX_PATH);
eval(pathstr)

cvx_setup;
eval(pathstr);
%cvx_quiet(true);

nMC = 8;
n = {0};
if n == 128
	nMC = 1;
end
if n == 64
	nMC = 2;
end
if n == 32
	nMC = 10;
end

Meas = '{1}';
Mat = '{2}';


mydate = date;
RandStream.setGlobalStream(RandStream('mt19937ar','seed', sum(100*clock)));
defaultStream = RandStream.getGlobalStream;
savedState = defaultStream.State;
fname = sprintf('randState{0}.mat');
save(fname, 'mydate', 'savedState');
r = [1:5, randi([6,n/4],1,2), randi([n/4+1,n/2],1,2), randi([n/2+1,n],1,2)];

%upper bound
u = n^2;
if strcmp(Meas,'Entry')
	% There are only this many entries
	u = n^2/2 + n/2;
end
if strcmp(Meas(2:end),'Dirac') && strcmp(Mat(1),'R')
	% There are only this many Diracs that are symmetric
	u = n^2/2 + n/2;
end
u = n^2/2 + n/2;
tuples = [];
for i_r = 1:length(r)
	rho = r(i_r)/n;
	l = rho;%round(2*rho - rho^2);
	l = round(rho - rho^2/2);
	for i_m = 1:10
		m = randi([n^2*l,u],1);
		tuples = [tuples; ...
				   m, r(i_r) ];
	end
end


filename = 'exp_n_{0}_Meas_{1}_Mat_{2}_suid_vminden_summaries.txt';

outputFile = fopen(filename,'w+');
result = run_experiment(outputFile,n,nMC,tuples,Meas,Mat);
fclose(outputFile);

filename = 'exp_n_{0}_Meas_{1}_Mat_{2}_suid_vminden_results.txt';
FID = fopen(filename,'w+');
for i_list = 1:length(tuples)
	fprintf(FID,'result(%i,:)= [%i %i %i %i ', i_list,n,tuples(i_list,2),tuples(i_list,1),nMC);
    fprintf(FID,'%6.3f %6.3f]\\n',...
       result(i_list,5),...
       result(i_list,6));
end
''' 


qsubfileTemplate = '''#!/bin/bash

#$ -N run{0}{1}{2}
#$ -o job{0}{1}{2}.out
#$ -e job{0}{1}{2}.error
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=156:00:00,h_vmem=20G

module load MATLAB-R2014a
matlab -nodesktop -singleCompThread < run{0}.m


'''

for A_type in Measurements:
	for X_type in Mats:
		if A_type[0] == 'C' and X_type[0] == 'R':
			print "Skipping combination (%s, %s) because it is not valid" % (X_type, A_type)
		else:
			# Make a subdirectory to put the output of this combination
			trial_dir = 'NewXtype%sAtype%s' % (X_type, A_type)
			try:
				os.mkdir(trial_dir)
			except:
				print 'Directory exists...'
			# For each of the matrix sizes, generate the appropriate Matlab script
			for n in Ns:
				if n == 128 and X_type[0] != 'HPSD':
					pass
				with open('%s/run%d.m' % (trial_dir, n), 'w') as runfile:
					with open('%s/run%d.sh' % (trial_dir, n), 'w') as qsubfile:
						# Create Matlab script
						matlabScript = matlabTemplate.format(n, A_type, X_type, trial_dir)
						runfile.write(matlabScript)
						qsubScript = qsubfileTemplate.format(n, A_type, X_type)
						qsubfile.write(qsubScript)
						#fi
				# filename = os.path.expanduser('~/CSProject/%s/run%d.sh' % (trial_dir,n))
				# print 'Printing: qsub ' + filename
				# #os.system('qsub ' + filename)
				with cd(os.getcwd() + '/' + trial_dir):
					subprocess.call(['qsub','run%d.sh' % n])
