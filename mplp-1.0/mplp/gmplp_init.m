function [gmplp_state] = gmplp_init(regions,intersects,lambda,orig_local,del,orig_adj,orig_lambda) %; %,lp_pair_bel)

gmplp_state.lambda = lambda;
gmplp_state.orig_local = orig_local;
gmplp_state.orig_lambda = orig_lambda;
gmplp_state.del = del;
gmplp_state.orig_adj = orig_adj;

fprintf('Entering GMPLP init\n');
tic;
% Sort regions (we need this for the intersections)
nV = length(orig_local);
for li=1:length(orig_local)
    var_size(li) = length(orig_local{li});
end

region_id = sparse(length(regions),nV);
for ri=1:length(regions)
    regions{ri} = sort(regions{ri});
    region_id(ri,regions{ri}) = 1;
    region_subsets{ri}=[];
end

intersect_id = sparse(length(intersects),nV);
for si=1:length(intersects)
    intersects{si} = sort(intersects{si});    
    intersect_id(si,intersects{si}) = 1; 
    intersect_size(si) = prod(var_size(intersects{si}));               
    intersect_with{si} = [];
end
intersect_lens = sum(intersect_id,2);

region_adj = region_id*region_id';
region_adj = region_adj-diag(diag(region_adj));


%nVals = var_size(li);
bel=0;


% For each intersection, specify the regions it intersects with
cluster_subsets = {};
all_bel_change = [];
all_int = [];
all_obj=[];

fprintf('First part of preparation took %g\n',toc);
tic;
[rlist_ii,rlist_jj,vv] = find(region_adj);
mask_for_max = 0;
mask = 0;
all_leninters = intersect_id*region_id';
% Go over regions, and find their intersections with the given sets
for ri = 1:length(regions)
    leninters = all_leninters(:,ri);;
    inters_withri = find(leninters==intersect_lens);
    for i=1:length(inters_withri)
        si = inters_withri(i);
        intersect_with{si}(end+1) = ri;
        region_subsets{ri}(end+1) = si;
        tmptmp=[];
        for ii=1:length(intersects{si})
            jj = find(intersects{si}(ii)==regions{ri});
            tmptmp(end+1) = jj;
        end
        inds_in_region{ri,length(region_subsets{ri})} = tmptmp;
        
%        mask{ri,length(region_subsets{ri})} = get_marginal_mask_multi1(var_size(regions{ri}),length(regions{ri}),inds_in_region{ri,length(region_subsets{ri})});        
%        tt = mask{ri,length(region_subsets{ri})} ;
%        tt(find(tt==0))=NaN;
%        mask_for_max{ri,length(region_subsets{ri})} = tt;        
    end
end



for si=1:length(intersects)
    intersect_with{si} = unique(intersect_with{si});
end
fprintf('Second part of preparation took %g\n',toc);
tic;
  
if ~exist('stop_on_int','var')
    stop_on_int  =1;
end

%rand('seed',0);
% Initialize messages
for ri=1:length(regions)
    % Go over all its intersections
    for sj=1:length(region_subsets{ri})
        curr_msgs{ri,sj} = 0*(rand(intersect_size(region_subsets{ri}(sj)),1)-0.5);
    end
    if length(regions{ri})==1
        curr_msgs{ri,1} = orig_local{regions{ri}(1)}(:);
        sum_into_subsets{region_subsets{ri}(1)} = orig_local{regions{ri}(1)}(:);
    end
end

% Get the sum of messages going into the intersection sets
sum_into_subsets = cell(length(intersects),1);
for si=1:length(intersects)
    tmp = 0;
    for iwi = 1:length(intersect_with{si})
        nregion = intersect_with{si}(iwi);
        index_in_region = find(region_subsets{nregion}==si);
        tmp = tmp + curr_msgs{intersect_with{si}(iwi),index_in_region};
    end
    sum_into_subsets{si} = tmp;
end

gmplp_state.sum_into_subsets = sum_into_subsets;
gmplp_state.intersect_with = intersect_with;
gmplp_state.region_subsets = region_subsets;
gmplp_state.regions = regions;
gmplp_state.regions_abstract = cell(size(regions));
gmplp_state.intersects = intersects;
gmplp_state.curr_msgs = curr_msgs;
gmplp_state.mask = mask;
gmplp_state.mask_for_max = mask_for_max;
gmplp_state.region_id = region_id;
gmplp_state.intersect_id = intersect_id;
gmplp_state.var_size = var_size;
gmplp_state.intersect_size = intersect_size;
gmplp_state.inds_in_region = inds_in_region;

fprintf('Last part of preparation took %g\n', toc);