function out=svm_cross_validation(labels,features,options)


svm_kernel=getoptions(options,'svm_kernel','rbf');

switch lower(svm_kernel)
    case 'poly'
    polynomial=1;
    otherwise
    polynomial=0;
end

nfold=getoptions(options,'cv_nfold',1);
alpha=getoptions(options,'cv_alpha',0.2);

if nfold==0
%best_c = getoptions(options,'svm_margin',2.82842712475);
%best_g=getoptions(options,'svm_gamma',0.00728932024638);
best_c = getoptions(options,'svm_margin',5.00);
best_g=getoptions(options,'svm_gamma',0.0018);
tempo=sprintf('-g %g -c %g -h 0', best_g, best_c);
out.model=tempo;
out.best_c=best_c;
out.best_g=best_g;
return;
end

L=length(labels);
LL=round(L*alpha);

steps=getoptions(options,'cv_steps',5);
depth=getoptions(options,'cv_depth',2);


if polynomial
cst=8/(steps-1);
c_grid=2.^[9:15]
g_grid=[1:1];
finegrid=2.^[-1/2:1/4:1/2];
else
cst=4/(steps-1);
c_grid=2.^[0:cst:4];
%c_grid=2.^[0:cst:8];
gst=4/(steps-1);
g_grid=2.^[-11:gst:-7];
finegrid=2.^[-1/2:1/4:1/2];
end


for n=1:nfold
%split randomly the data into validation and training
	fprintf('--cv fold %d \n', n)
	rien=rand(1,L);
	[ri,plant]=sort(rien);
	lab_tmp=labels(plant);
	feat_tmp=features(plant,:);
        perf=zeros(1,length(g_grid)*length(c_grid));
	parfor ir=1:length(c_grid)*length(g_grid)
	%for ir=1:length(c_grid)*length(g_grid)
	%for ic=1:length(c_grid)
	%		for ig=1:length(g_grid)
                        ic=1+floor((ir-1)/length(g_grid));
                        ig=1+mod(ir-1,length(g_grid));
                          if polynomial
				tempo=sprintf('-t 1 -d %g -c %g -h 0 -q', g_grid(ig), c_grid(ic));
                          else
				tempo=sprintf('-g %g -c %g -h 0 -q', g_grid(ig), c_grid(ic));
                          end
				svm_model = svmtrain(lab_tmp(LL+1:end),feat_tmp(LL+1:end,:),tempo);
				[o1,o2,o3]=svmpredict(lab_tmp(1:LL),feat_tmp(1:LL,:),svm_model);
				%perf(length(g_grid)*(ic-1)+ig)=o2(1);
				perf(ir)=o2(1);
                                fprintf('obtained with parameters ic:%g ig:%g \n',c_grid(ic),g_grid(ig)) 
	%		end
	end
	[optperf,pos]=max(perf);
	opt_c=1+floor((pos-1)/length(g_grid));
	opt_g=1+mod(pos-1,length(g_grid));
	chosen_c(n)=c_grid(opt_c);
	chosen_g(n)=g_grid(opt_g);
        if polynomial 
        fineperf=zeros(1,length(finegrid));
        else
        fineperf=zeros(1,length(finegrid)*length(finegrid));
        end  
	%fine grid optimization
	if(depth>1)
		%find optimium on coarse grid
		fprintf('--cv fine tune %d \n', n)
		fine_c_grid=c_grid(opt_c)*finegrid;
                if polynomial 
		fine_g_grid=g_grid(opt_g);
                else
		fine_g_grid=g_grid(opt_g)*finegrid;
                end
		parfor ir=1:length(fine_c_grid)*length(fine_g_grid)
		%for ir=1:length(fine_c_grid)*length(fine_g_grid)
		%for ic=1:length(fine_c_grid)
		%	for ig=1:length(fine_g_grid)
                        ic=1+floor((ir-1)/length(fine_g_grid));
                        ig=1+mod(ir-1,length(fine_g_grid));	
                                if polynomial
				tempo=sprintf('-t 1 -d %g -c %g -h 0 -q', fine_g_grid(ig), fine_c_grid(ic));
                                else
				tempo=sprintf('-g %g -c %g -h 0 -q', fine_g_grid(ig), fine_c_grid(ic));
                                end
				svm_model = svmtrain(lab_tmp(LL+1:end),feat_tmp(LL+1:end,:),tempo);
				[o1,o2,o3]=svmpredict(lab_tmp(1:LL),feat_tmp(1:LL,:),svm_model);
				fineperf(ir)=o2(1);
		%	end
		end
		[optperf,pos]=max(fineperf);
		opt_c=1+floor((pos-1)/length(fine_g_grid));
		opt_g=1+mod(pos-1,length(fine_g_grid));	
		chosen_c(n)=fine_c_grid(opt_c);
		chosen_g(n)=fine_g_grid(opt_g);
	end

end

best_c=mean(chosen_c);
best_g=mean(chosen_g);

if polynomial
tempo=sprintf('-t 1 -d %g -c %g -h 0', best_g, best_c);
else
tempo=sprintf('-g %g -c %g -h 0', best_g, best_c);
end

out.model=tempo;
out.best_c=best_c;
out.best_g=best_g;


