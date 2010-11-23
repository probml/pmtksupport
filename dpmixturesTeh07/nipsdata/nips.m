% load up nips data
load nips12raw_str602
% author list
alist1 = [
1855     %Viola_P
506      %Freeman_W
...
850      %Kearns_M
1735     %Sutton_R
1770     %Tesauro_G
1780     %Thrun_S
1800     %Touretzky_D
...
553      %Giles_C
961      %LeCun_Y
1637     %Simard_P
175      %Bottou_L
...
1837     %Vapnik_V
1571     %Scholkopf_B
1666     %Smola_A
1618     %Shawe-Taylor_J
...
1766     %Tenenbaum_J
1054     %MacKay_D
1925     %Williams_C
1424     %Rasmussen_C
108      %Bengio_Y
703      %Hinton_G
132      %Bishop_C
546      %Ghahramani_Z
810      %Jordan_M
1549     %Saul_L
783      %Jaakkola_T
...
1391     %Pouget_A
887      %Koch_C
...
2000     %Zemel_R
361      %Dayan_P
1600     %Sejnowski_T
];

alist = alist1;                  % indices of authors

numa = length(alist);            % number of authors
numw = 500;                      % number of word types
totp = length(ptitles);          % total number of papers

apmap = zeros(numa,totp);
for aa = 1:numa
  apmap(aa,apapers{alist(aa)}) = 1;
end
plist = find(sum(apmap,1));      % indices of papers
nump  = length(plist);           % number of papers
apmap = apmap(:,plist);          % entry (i,j)=1 if author i wrote paper j

wpcounts = counts(:,plist);      % counts of words in papers
pcounts = full(sum(wpcounts,1)); % lengths of papers

% pick informative words.
% estimate informativeness by first do tf-idf, then compute entropy.
ww = wpcounts ./ log(1+pcounts(ones(size(wpcounts,1),1),:));
np = length(plist);
h = zeros(1,1000);
for i=1:1000
  pi = 10/np + ww(i,:);
  pi = pi / sum(pi);
  h(i) = -sum(pi.*log(pi));
end
[a i] = sort(h);
wlist = i(1:numw);              % indices of numw most informative words

wpcounts = counts(wlist,plist);  % counts of words in papers
wcounts = full(sum(wpcounts,2)); % counts of word tokens
pcounts = full(sum(wpcounts,1)); % lengths of papers

papers  = ptitles(plist);
authors = anames(alist);
words   = wl(wlist);
