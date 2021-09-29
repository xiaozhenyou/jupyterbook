function [dynmat,dynmatsd]=dynamic(inputmat, window,step, motioncovar)
% assume inputmat is parcel x time (TRs)matrx or a txt file contains the matrix
%  slide window is 1 TR or anything
%motioncovar is the 1xtime vector such as FD vector can also be a txt file,
%or fmriprep output
%sub-DCNLRDOC002_task-SocStroop1_desc-confounds_regressors.tsv
%contains the values
% example input
inputmat=rand(400,180);
motioncovar=rand(1,180);
%motioncovar='sub-DCNLRDOC004_task-rest_desc-confounds_regressors.tsv';
window=10;
step=1;
% Xiaozhen You 2021 September 10
if ischar(inputmat)
    inputmat=load(inputmat);
end
if ischar(motioncovar)
    if strcmp(motioncovar(end-2:end),'tsv')
        motioncovar=readtable(motioncovar,'FileType','text');
        motioncovar=motioncovar.framewise_displacement;
        motioncovar(1)=0;
        
    else
        motioncovar=load(motioncovar);
    end
end
totalwindownumber=round((size(inputmat,2)-window)/step);
parcelnumber=size(inputmat,1);

dynmat=zeros(parcelnumber,parcelnumber,totalwindownumber);
dynmatsd=zeros(parcelnumber,parcelnumber);
motionmat=zeros(totalwindownumber,1);
for ti=1:totalwindownumber
    
    selmat=inputmat(:,(ti-1)*step+1:(ti-1)*step+window);
    dynmat(:,:,ti)=corrcoef(selmat');
    
    motionmat(ti)=mean(motioncovar((ti-1)*step+1:(ti-1)*step+window));
end

% regress out motion car
for pari=1:parcelnumber
    for parj=1:parcelnumber
        dynmat(pari,parj,:)=regress(squeeze(dynmat(pari,parj,:)),motionmat);
        dynmatsd(pari,parj)=std(dynmat(pari,parj,:));
    end
end
figure;imagesc(dynmatsd);

end