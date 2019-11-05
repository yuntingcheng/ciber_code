function [p,C] = profile_cov_binning(pin,Cin,win,idx)
    % bin idx in prof w/ Covariance C into one bin
    Nin = numel(pin);
    Nout = Nin - numel(idx) + 1;
    
    w = win(idx);
    w = w./sum(w);
    C = Cin;
    
    cbin = zeros([1,Nin]);
    for i=1:numel(idx)
        cbin = cbin + C(idx(i),:).*w(i);
    end
    C(idx(1),:) = cbin;
    C(idx(2:end),:) = [];
    
    cbin = zeros([Nout,1]);
    for i=1:numel(idx)
        cbin = cbin + C(:,idx(i)).*w(i);
    end
    C(:,idx(1)) = cbin;
    C(:,idx(2:end)) = [];
    
    p = pin;
    p(idx(1)) = sum(p(idx).*w);
    p(idx(2:end)) = [];
    
return
