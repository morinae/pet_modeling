function check_theoretical_assumptions(K,Mu,Lam,t) % see Theorem 15 in Paper
    assert(size(K,2)==3,'There must be exactly three types of metabolic parameters'); % three types of metabolic parameters
    p = size(Mu,2); % degree of plasma concentration
    n = size(K,1);  % number of regions
    T = size(t,2);  % number of times
    assert(size(unique(t),2)==T,'Time points must be distinct');            % distinct times
    assert(n>=3,'Need at least three regions');                             % minimal number of regions
    assert(p>=3,'C_P must be at least triexponential');                     % minimal degree of plasma concentration
    assert(prod((K>0),'all')==1,'Metabolic parameters must be positive');   % positive metabolic parameters
    assert(prod((Mu~=0),'all')==1,'C_P exponents must be nonzero');         % nonzero exponents of plasma concentration
    % Sufficient condition for Assumption (A) in Paper
    k3 = K(:,3);            
    k23 = k3+K(:,2);
    [~,w3]=unique(k3,'stable');
    [~,w23]=unique(k23,'stable');
    w = intersect(w3,w23);
    nw = size(w,1);
    assert(nw>=3,'At least three regions must exist where the k_3 and the k_2+k_3 are pairwise distinct, respectively');
    mask = logical(dec2bin(sum(nchoosek(2.^(0:nw-1),3),2))-'0');
    m = size(mask,1);
    df = zeros(1,p);
    for i=1:m
       wi = w(mask(i,:));
       k3i = k3(wi);
       k23i = k23(wi);
       ind1 = prod((Mu+k3i~=0));
       ind2 = (Mu+k23i~=0);
       S=(Lam.*ind2)./(k23i+Mu+(1-ind2));
       df = df+(prod((sum(S')~=0)'+1-ind2)).*ind1;
    end
    df = (df~=0);
    df = prod(df);
    assert(df==1,'Assumption (A) is not fulfilled')         % assumption (A)
    assert((T>=2*(p+3))==1,'Not enough time evaluations');  % minimal number of time evaluations
    disp('Theoretical assumptions are fulfilled')
end