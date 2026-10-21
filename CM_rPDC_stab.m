function CM = CM_rPDC_stab(CM,xy,bad)

%  CM :         structure with the following fields:
%  .p :         order of the MVAR model
%  .ds :        downsampling factor
%  .Fmax :      maximal frequency inculded in the analysis
%  .Fs :        sampling frequency
%  .to_compute: list of parameters to be computed
%     .Coh :      Coherence
%     .PDC :      Partial Directed Coherence
%     .rPDC :     renormalized Partial Directed Coherence
%  .show_plots: whether to show plots (1) or not (0). (default: 1)
%  .stat :      informations about stats
%    .alpha :     alpha value (default: 0.05)
%
%  xy :         data matrix whose dimensions correspond to 
%                  number of signals x number of time points 
%
%  bad :        a vector with the same number of time points as xy,
%               informing for each time point whether it is bad (value 1) 
%               or not (value 0).



% set up analysis parameters
Nsig = size(xy,1);      % number of signals
Nt = size(xy,2);        % number of time points
Nep = 1;                % number of epochs (1 if non-epoched data)

if nargin == 2
    bad = zeros(1,length(xy));
end

CM.freq = (0:CM.p-1)/CM.p*CM.Fs/CM.ds;
if isfield(CM,'Fmax')
    CM.keepF = find(CM.freq <= CM.Fmax);
else
    CM.keepF = 1:length(CM.freq);
end
if ~isfield(CM,'to_compute')
    CM.to_compute.Coh = 1;
    CM.to_compute.PDC = 0;
    CM.to_compute.rPDC = 1;
else
    if ~isfield(CM.to_compute,'Coh')
        CM.to_compute.Coh = 1;
    end
    if ~isfield(CM.to_compute,'PDC')
        CM.to_compute.PDC = 0;
    end
    if ~isfield(CM.to_compute,'rPDC')
        CM.to_compute.rPDC = 1;
    end
end
if ~isfield(CM,'show_plots')
    CM.show_plots = 1;
end
if ~isfield(CM,'stat')
    CM.stat.alpha = 0.05;
end



% low pass filter the signals and resample
if CM.ds ~= 1
    F_h_x   = cosine_filter(Nt,{'low'},1/CM.ds,1/CM.ds/10);
    Fxy = fft(xy,[],2).*repmat(F_h_x,[Nsig 1 Nep]);
    xy = real(ifft(Fxy,[],2));
    xy = xy(:,round(1:CM.ds:Nt),:);
    bad = reshape(bad,[CM.ds Nt/CM.ds]);
    bad = max(bad,[],1);
    Nt = size(xy,2);        % number of time points
end
alk = find(diff([1 bad]) == -1);
lop = find(diff([bad 1]) == 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  simple coherence analysis                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lxy = fix(Nt/CM.p)*CM.p;
if CM.to_compute.Coh
    % Compute coherence
    ovlp = linspace(1,CM.p+1,6);
    ovlp = round(ovlp(1:5));
    tmp = [];
    bad_coh = [];
    for k = 1:length(ovlp)
        Lxytmp = Lxy-mod(Lxy-ovlp(k)+1,CM.p);
        tmp = cat(3,tmp,reshape(xy(:,ovlp(k):Lxytmp,:),Nsig,CM.p,(Lxytmp-ovlp(k)+1)/CM.p,Nep));
        bad_coh = cat(2,bad_coh,max(reshape(bad(1,ovlp(k):Lxytmp),CM.p,(Lxytmp-ovlp(k)+1)/CM.p)));
    end
    tmp = tmp(:,:,find(~bad_coh));
    
    Ftmp = fft(tmp,[],2);
    CSD = zeros(Nsig,Nsig,CM.p);
    for f = 1:CM.p
        CSD(:,:,f) = squeeze(Ftmp(:,f,:))*squeeze(Ftmp(:,f,:))';
    end
    if isfield(CM,'smooth')
        CSD = (CSD(:,:,[end-1:end 1:end-2]) + CSD(:,:,[end 1:end-1]) + CSD + CSD(:,:,[2:end 1]) + CSD(:,:,[3:end 1:2]))/5; %%%%% to remove (smoothing)
    end
    S2 = CM.p;
    CM.nave_coh = size(tmp,3);
    CSD = CSD/size(tmp,3)/S2/(CM.Fs/CM.ds)*2;
    % sqrt(CSD(1,1,1)/CM.p^2/CM.nave_coh)
    % sqrt(CSD(1,1,1)/CM.p^2/CM.nave_coh)/sqrt(CM.Fs/CM.p)
    
    CM.CSD = CSD;
    CM.Coh = abs(CSD).^2;
    for ind = 1:Nsig
        CM.Coh(ind,:,:) = CM.Coh(ind,:,:)./repmat(abs(CSD(ind,ind,:)),[1 Nsig]);
        CM.Coh(:,ind,:) = CM.Coh(:,ind,:)./repmat(abs(CSD(ind,ind,:)),[Nsig 1]);
    end
    
    % display results
    if CM.show_plots
        figure
        scale = max(CM.Coh(CM.Coh < 0.999));
        for ind1 = 1:Nsig
            for ind2 = 1:Nsig
                if ind1 ~= ind2
                    subplot(Nsig,Nsig,Nsig*(ind1-1)+ind2)
                    hold off
                    plot(CM.freq(CM.keepF),squeeze(CM.Coh(ind1,ind2,CM.keepF)))
                    title(['Coh ' num2str(ind1) ' and ' num2str(ind2)])
                    axis([0 CM.freq(CM.keepF(end)) 0 scale])
                end
            end
        end
    end
end

% remove few samples at the end of the signals for convenience...
Lxy = fix(Nt/CM.p)*CM.p;
xy = xy - repmat(mean(xy(:,find(~bad)),2),[1 Nt Nep]);
xy = xy./repmat(std(xy(:,find(~bad)),[],2),[1 Nt Nep]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Prepare PDC/rPDC analysis (fit the MVAR model)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CM.to_compute.PDC | CM.to_compute.rPDC
    pmin = CM.p; pmax = CM.p;

    xy_cell = [];
    L_cell = [];
    for t = 1:length(alk)
        xy_cell{t} = permute(xy(:,alk(t):lop(t),:),[2 1 3]);
        L_cell(t) = size(xy_cell{t},1);
    end
    xy_cell(find(L_cell < CM.p)) = [];


    %%%% to estimate the model order
    % [west, Aest, Cest, CM.SBC, CM.FPE, th] = arfit_stab(xy_cell, 1, 100);
    %%%%
    [west, Aest, Cest, SBC, FPE, th] = arfit_stab(xy_cell, pmin, pmax);
    A = reshape(Aest,[Nsig Nsig CM.p]);

    FA = zeros(Nsig,Nsig,CM.p);
    for f = 0:CM.p-1
        FA(:,:,f+1) = eye(Nsig);
        for k = 1:CM.p
            FA(:,:,f+1) = FA(:,:,f+1) - A(:,:,k)*exp(-f*k*2i*pi/CM.p);
        end
    end

    if isfield(CM,'smooth')
        FA = (FA(:,:,[end-1:end 1:end-2]) + FA(:,:,[end 1:end-1]) + FA + FA(:,:,[2:end 1]) + FA(:,:,[3:end 1:2]))/5; %%%%% to remove (smoothing)
    end
    CM.FA = FA;
end

% FA = repmat(eye(Nsig),[1 1 p]) - fft(A,[],3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        partial directed coh Baccala and Sameshima (2001)        %
%            also described in Schelter (2006,2009)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CM.to_compute.PDC
    CM.PDC = zeros(size(FA));
    for ind = 1:Nsig
        % formula described in Baccala and Sameshima (2001) or Schelter (2006)
        CM.PDC(:,ind,:) = abs(FA(:,ind,:))./repmat(sqrt(sum(abs(FA(:,ind,:)).^2,1)),[Nsig 1]);

        % formula proposed by Schelter (2009)
        % CM.PDC(ind,:,:) = abs(FA(ind,:,:))./repmat(sqrt(sum(abs(FA(ind,:,:)).^2,2)),[1 Nsig]);
    end

    % display results
    if CM.show_plots
        figure
        title('PDC analysis')
        for ind1 = 1:Nsig
            for ind2 = 1:Nsig
                if ind1 ~= ind2
                    subplot(Nsig,Nsig,Nsig*(ind1-1)+ind2)
                    plot(CM.freq(CM.keepF),squeeze(CM.PDC(ind1,ind2,CM.keepF)))
                    title(['PDC ' num2str(ind1) ' <- ' num2str(ind2)])
                    axis([0 CM.freq(CM.keepF(end)) 0 1])
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        renormalized PDC described in Schelter (2006,2009)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CM.to_compute.rPDC
    % parameters used in rPDC estimation
    Sigma = Cest;
    R = zeros(Nsig*CM.p,Nsig*CM.p);
    for k = 1:CM.p
        xy_k = [];
        for t = 1:length(alk)
            xy_k = [xy_k xy(:,(CM.p+alk(t):lop(t))-k,:)];
        end
        for l = 1:CM.p
            xy_l = [];
            for t = 1:length(alk)
                xy_l = [xy_l xy(:,(CM.p+alk(t):lop(t))-l,:)];
            end
            R((-Nsig+1:0)+k*Nsig,(-Nsig+1:0)+l*Nsig) = xy_k(:,:)*xy_l(:,:)'/sum(lop-alk-CM.p+1)/Nep;
        end
    end
    H = pinv(R,1e-5);
    V = zeros(2,2,Nsig,Nsig,CM.p);
    c = cos((1:CM.p)'*(0:CM.p-1)*2*pi/CM.p);
    s = sin((1:CM.p)'*(0:CM.p-1)*2*pi/CM.p);
    CS = cat(2, cat(1,repmat(reshape(c,[1 1 CM.p 1 CM.p]),[1 1 1 CM.p 1]).*repmat(reshape(c,[1 1 1 CM.p CM.p]),[1 1 CM.p 1 1]), ...
                      repmat(reshape(s,[1 1 CM.p 1 CM.p]),[1 1 1 CM.p 1]).*repmat(reshape(c,[1 1 1 CM.p CM.p]),[1 1 CM.p 1 1])), ...
                cat(1,repmat(reshape(c,[1 1 CM.p 1 CM.p]),[1 1 1 CM.p 1]).*repmat(reshape(s,[1 1 1 CM.p CM.p]),[1 1 CM.p 1 1]), ...
                      repmat(reshape(s,[1 1 CM.p 1 CM.p]),[1 1 1 CM.p 1]).*repmat(reshape(s,[1 1 1 CM.p CM.p]),[1 1 CM.p 1 1])));
    for ix = 1:Nsig
        for jx = 1:Nsig
            V(:,:,ix,jx,:) = V(:,:,ix,jx,:) + sum(sum(repmat(reshape(H(jx+(0:CM.p-1)*Nsig,jx+(0:CM.p-1)*Nsig),[1 1 CM.p CM.p]),[2 2 1 1 CM.p]).*CS*Sigma(ix,ix),3),4);
            % V(:,:,i,j,w) equ 10 (Schelter 2009)
        end
    end

    % compute rPDC
    CM.rPDC = zeros(Nsig,Nsig,CM.p);
    scale = 0;
    for f = 0:CM.p-1
        for ix = 1:Nsig
            for jx = 1:Nsig
                X = [real(FA(ix,jx,f+1)) ; imag(FA(ix,jx,f+1))];
                CM.rPDC(ix,jx,f+1) = X'*pinv(V(:,:,ix,jx,f+1),1e-10)*X;
                if ix ~= jx
                    scale = max(scale,max(CM.rPDC(ix,jx,f+1)));
                end
            end
        end
    end

    % compute statistics
    % CM.stat.Th_rPDC = chi2inv(1-CM.stat.alpha,2)/(Nt-CM.p)/Nep;
    CM.nfree = sum(L_cell-CM.p+1);
    CM.stat.Th_rPDC = chi2inv(1-CM.stat.alpha,2);
    CM.rPDC = CM.rPDC*CM.nfree;
    scale = scale*CM.nfree;

    % display results
    if CM.show_plots
        figure
        for ind1 = 1:Nsig
            for ind2 = 1:Nsig
                if ind1 ~= ind2
                    subplot(Nsig,Nsig,Nsig*(ind1-1)+ind2)
                    plot(CM.freq(CM.keepF),squeeze(CM.rPDC(ind1,ind2,CM.keepF)))
                    title([num2str(ind1) ' <- ' num2str(ind2)])
                    hold on
                    plot(CM.freq(CM.keepF([1 end])),squeeze(CM.stat.Th_rPDC*[1 1]),'r')
                    title(['rPDC ' num2str(ind1) ' <- ' num2str(ind2)])
                    axis([0 CM.freq(CM.keepF(end)) 0 max(scale,CM.stat.Th_rPDC)*1.1])
                end
            end
        end
    end
end


