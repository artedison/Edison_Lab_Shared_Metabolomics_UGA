function [xCS,ints,ind,target] = icoshift(xT,xP,inter,n,options,Scal)
% interval Correlation Optimized shifting
%
% [xCS,ints,ind,target] = icoshift(xT,xP,inter[,n[,options[,Scal]]])
%
% Splits a spectral database into "inter" intervals and coshift each vector
% left-right to get the maximum correlation toward a reference or toward an
% average spectrum in that interval. Missing parts on the edges after
% shifting are filled with "closest" value or with "NaNs".
%
% INPUT
% xT (1 × mT)    : target vector.
%                  Use 'average' if you want to use the average spectrum as a reference
%                  Use 'median' if you want to use the median spectrum as a reference
%                  Use 'max' if you want to use for each segment the corresponding actual
%                  spectrum having max features
% xP (nP × mT)   : Matrix of sample vectors to be aligned as a sample-set
%                  towards common target xT
% inter          : definition of alignment mode
%                  'whole'         : it works on the whole spectra (no intervals).
%                  nint            : (numeric) number of many intervals.
%                  'ndata'         : (string) length of regular intervals
%                                    (remainders attached to the last).
%                  [I1s I1e,I2s...]: interval definition. ('I(n)s' interval
%                                    n start, 'I(n)e' interval n end).
%                  (refs:refe)     : shift the whole spectra according to a
%                                    reference signal(s) in the region
%                                    refs:refe (in sampling points)
%                  'refs-refe'     : shift the whole spectra according to a
%                                    reference signal(s) in the region
%                                    refs-refe (in Scal units)
% n (1 × 1)      : (optional)
%                  n = integer n.: maximum shift correction in data
%                                  points/Scal units (cf. options(5))
%                                  in x/rows. It must be >0
%                  n = 'b' (best): the algorithm search for the best n
%                                  for each interval (it can be time consuming!)
%                  n = 'f' (fast): fast search for the best n for each interval (default)
%                  A warning is displayed for each interval if "n" appears too small
% options (1 × 5): (optional)
%                  (1) triggers plots & warnings:
%                      0 : no on-screen output
%                      1 : only warnings (default)
%                      2 : warnings and plots
%                  (2) selects filling mode
%                      0 : using not a number
%                      1 : using previous point (default)
%                  (3) turns on Co-shift preprocessing
%                      0 : no Co-shift preprocessing (default)
%                      1 : Executes a Co-shift step before carrying out iCOshift
%                  (4) max allowed shift for the Co-shift preprocessing (default = equal to n if not specified)
%                      it has to be given in Scal units if option(5)=1
%                  (5) 0 : intervals are given in No. of datapoints  (deafult)
%                      1 : intervals are given in ppm --> use Scal for inter and n
% Scal           : vector of scalars used as axis for plot (optional)
%
% OUTPUT
% xCS  (nP × mT): shift corrected vector or matrix
% ints (nI × 3) : defined intervals (No., starting point, ending point)
% ind  (nP × nI): matrix of indexes reporting how many points each spectrum
%                 has been shifted for each interval
% target (1 x mP): actual target used for the final alignment
%
% Authors:
% Francesco Savorani - Department of Food Science
%                      Quality & Technology - Spectroscopy and Chemometrics group
%                      Faculty of Life Sciences
%                      University of Copenhagen - Denmark
% email: frsa@life.ku.dk - www.models.life.ku.dk
% Giorgio Tomasi -     Department of Basic Science and Environment
%                      Soil and Environmental Chemistry group
%                      Faculty of Life Sciences
%                      University of Copenhagen - Denmark
% email: gto@life.ku.dk - www.igm.life.ku.dk

% 170508 (FrSa) first working code
% 211008 (FrSa) improvements and bugs correction
% 111108 (Frsa) Splitting into regular intervals (number of intervals or wideness in datapoints) implemented
% 141108 (GT)   FFT alignment function implemented
% 171108 (FrSa) options implemented
% 241108 (FrSa) Automatic search for the best or the fastest n for each interval implemented
% 261108 (FrSa) Plots improved
% 021208 (FrSa) 'whole' case added & robustness improved
% 050309 (GT)   Implentation of interpolation modes (NaN); Cosmetics; Graphics
% 240309 (GT)   Fixed bug in handling missing values
% 060709 (FrSa) 'max' target and output 'target' added. Some speed, plot and robustness improvements
% 241109 (GT)   Interval and band definition in units (added options(5))
% 021209 (GT)   Minor debugging for the handling of options(5)
% 151209 (FrSa) Cosmetics and minor debugging for the handling of options(5)
% 151110 (FrSa) Option 'Max' works now also for alignment towards a reference signal
% 310311 (FrSa) Bugfix for the 'whole' case when mP < 101
%% ERRORS CHECK
% Constant
BLOCKSIZE = 2^25; % To avoid out of memory errors when 'whole', the job is divided in blocks of 32MB

% Number of input arguments
if ~nargin
   
   fprintf('\n\n icoshift(target,dataset,alignment mode,max allowed shift,options,scale);')
   fprintf('\n or')
   fprintf('\n icoshift(target,dataset,alignment mode);\n')
   fprintf('\n Examples:')
   fprintf('\n icoshift(''median'',xP,20);\n')
   fprintf('splits in 20 regular intervals using the default settings and align them towards the median spectrum\n')
   fprintf('\n icoshift(''average'',xP,''200'',''40'');\n')
   fprintf('splits in regular intervals made of 200 data points using the default settings\n')
   fprintf('and align them towards the average spectrum allowing a 40 points max shift\n')
   fprintf('\n icoshift(target,xP,(1000:2000),''b'',[2 1 1]);\n')
   fprintf('aligns the whole spectra according to the reference signal in target(1000:2000)\n')
   fprintf('searching for the best max allowed shift for each interval and performing a pre-coshift step.\n')
   fprintf('Plots the achieved results\n')
   fprintf('\n icoshift(''average'',xP,user_ints,''f'',[2 0],ppm_scale);\n')
   fprintf('splits the dataset into the intervals defined by the user (user_ints) and aligns them\n')
   fprintf('towards the average spectrum using the highest allowed shift (fastest method) for each interval\n')
   fprintf('and filling the interval edges with NaNs. Plots the achieved results using the ppm_scale provided\n')
   fprintf('\n For further specifications type ''help icoshift''\n')
   return
   
end
if nargin < 3, error('Not enough input arguments'); end
if nargin < 4 || isempty(n)
   n = 'f'; % Change here the default assigned "n" value
   fprintf('\n Fast search for the best "n" set by default (initial value 50)\n')
end
if nargin < 5, options = []; end
if nargin < 6, Scal = 1:size(xP,2); end

%% Check input
% Scal #GT To handle options(5) properly
if numel(Scal) ~= length(Scal), error('Scal must be a vector'), end
if length(Scal) ~= size(xP,2), error('X and Scal are not of compatible length'), end
dScal   = diff(Scal);
incScal = Scal(1) - Scal(2);
if incScal == 0 || ~all(sign(dScal) == -sign(incScal)), error('Scal must be strictly monotonic'); end
flag_ScalDir = incScal < 0;
flag_DiScal  = any(abs(dScal) > 2 * min(abs(dScal))); % #GT Check for discontinuous Scal vector

% xT and xP
max_flag = 0;
if  strcmpi(xT,'average')
   xT   = MeanwNaN(xP);
   note = 'average';
elseif strcmpi(xT,'median')
   xT = nanmedian(xP);
   note = 'median';
elseif strcmpi(xT,'max')        % MAX
   xT = zeros(1,size(xP,2));    % MAX
   max_flag = 1;                % MAX
   note = 'max';                % MAX
else note = 0;
end
[nT,mT] = size(xT);
[nP,mP] = size(xP);
if (mT ~= mP), error('Target "xT" and sample "xP" must have the same number of columns'); end
if isnumeric(inter) %FRSA bugfix for 'whole' case
if inter > mP
   error('ERROR: number of intervals "inter" must be smaller than number of variables in xP');
end
end

% Options
% options_co   = [1 1 0 0 nargin > 5];
options_co   = [1 1 0 0 0]; % FrSa: default options(5) has to be = 0
if length(options) < length(options_co), options(end + 1:length(options_co)) = options_co(length(options) + 1:end); end
if options(1) < 0 || options(1) > 2,   error('options(1) must be 0, 1 or 2'); end
if options(5), 
   prec = abs(min(unique(dScal))); 
   if flag_DiScal, 
      warning('iCOshift:discontinuous_Scal','Scal vector is not continuous, the defined intervals might not reflect actual ranges'); 
   end
end

% Execute preliminary coshift if so required
flag_coshift = ~strcmpi(inter,'whole') && options(3);
if flag_coshift
   
   if options(4) == 0
      n_co = n;
   else
      n_co = options(4);
   if nargin >= 6 && options(5), n_co = dscal2dpts(n_co,Scal,prec); end % GT To handle scalar input
   end
   if max_flag             % MAX
      xT = MeanwNaN(xP);   % MAX
   end                     % MAX
   xPo = xP;    %FRSA save the original xP for plotting purpose
   [xP,nil,wint] = icoshift(xT,xP,'whole',n_co,[0 1 0]); %#ok<*ASGLU> %FRSA preprocessing COhift doesn't use NaNs
   title_suff    = ' (after coshift)';
   if  strcmpi(note,'average')
      xT = nanmean(xP); %FRSA uses the new CO-shifted average as a target
   end
   if  strcmpi(xT,'median')
      xT = nanmedian(xP); %FRSA uses the new CO-shifted median as a target
   end
   
else
   title_suff = '';
end

% "Inter" % GT: handle scalars automatically
flag_allsegs = false;
whole        = false;
flag2        = false;
if ischar(inter)
   
   if  strcmpi(inter,'whole')
      inter = [1,mP];
      whole = true;
   else
      
      if ~any(inter == '-')
      
         interv = str2double(inter);
         if nargin < 6 || ~options(5), interv = round(interv);
         else interv = dscal2dpts(interv,Scal,prec);
         end
         [inter]      = defints(xP,interv,options(1));
         flag_allsegs = true;
      
      else
         
         interv = regexp(inter,'(-{0,1}\d*\.{0,1}\d*)-(-{0,1}\d*\.{0,1}\d*)','tokens');
         interv = sort(scal2pts(str2double(cat(1,interv{:})),Scal,prec));
         if numel(interv) ~= 2, error('Invalid range for reference signal'), end
         inter = interv(1):interv(2);
         flag2 = true;
         
      end
      
   end
   
elseif isa(inter,'numeric')
   
   if numel(inter) == 1
      
      if fix(inter) == inter
         
         flag_allsegs = true;
         remain       = mod(mP,inter);
         N            = fix(mP/inter);
         % Distributes vars_left_over in the first "vars_left_over" intervals
         startint     = [(1:(N+1):(remain - 1)*(N+1)+1)'; ((remain - 1) * (N + 1) + 1 + 1 + N:N:mP)'];
         % endint=[startint(2:inter)-1; mP]; %ints with adjacent edges
         endint       = [startint(2:inter); mP]; % ints with common edges
         inter        = [startint endint]';
         inter        = inter(:)';
         
      else error('The number of intervals must be an integer');
      end
      % else intdif = diff(inter);
      %     if any(intdif(2:2:length(intdif)) < 0)
      %         uiwait(msgbox('The user-defined intervals are overlapping: is it intentional?','Warning','warn'));
      %     end
      
   else      
      flag2 = isequal(fix(inter),inter) && length(inter) > 1 && ...
         isequal([1,max(inter) - min(inter) + 1],size(inter)) && isequal(unique(diff(inter,1,2)),1);
      if ~flag2 && options(5), 
         inter = scal2pts(inter,Scal,prec);
         if any(inter(1:2:end) > inter(2:2:end)) && ~flag_ScalDir 
            inter = flipud(reshape(inter,2,length(inter)/2));
            inter = inter(:)';
         end
      end % Converts scalars to points 
   end 
   
end

[nint,mint] = size(inter);
wmsg = [];
scfl = isequal(fix(Scal),Scal) & ~options(5);
if ischar(n) && (numel(n) > 1 || ~any(ismember('bf',lower(n))))
   error('"n" must be a scalar/''f''/''b''')
elseif isa(n,'numeric')
   
   if any(n <= 0),  error('ERROR: shift(s) "n" must be larger than zero'); end
   if numel(n) > 1, wmsg = sprintf('"n" must be a scalar/character; first element (i.e. %i) used',round(n)); end
   if scfl && n ~= fix(n),  wmsg = sprintf('"n" must be an integer if Scal is ignorde; first element (i.e. %i) used',round(n));
   elseif options(5), n = dscal2dpts(n,Scal,prec);
   end
   if ~flag2 && any(diff((reshape(inter,2,mint/2)),1,1) < n)
      error('ERROR: shift "n" must be not larger than the size of the smallest interval');
   end
   n = round(n(1));
   if ~isempty(wmsg)
      warning('iCoShift:Input',wmsg)
      fprintf('\n press a key to continue...')
      pause
   end
   
end

%% Missing values
flag     = isnan(cat(1,xT,xP));
frag     = false;
Ref      = @(c) reshape(c,2,length(c)/2)';
vec      = @(A) A(:);
[mi,pmi] = min(inter);
[ma,pma] = max(inter);
if any(vec(flag))
   
   % If there are missing values in the target, check that the same pattern
   % of missing is present in the samples. If not an error is produced.
   if isequal(flag(ones(nP,1),:),flag(2:end,:))
      Select = @any;
   else
      Select = @all;
   end
   if flag2
      
      intern = RemoveNaN([1,pma - pmi + 1],cat(1,xT(:,inter),xP(:,inter)),Select);
      if size(intern,1) ~= 1, error('Reference region contains a pattern of missing that cannot be handled consistently')
      elseif ~isequal(intern,[1,inter(end) - inter(1) + 1]), warning('iCoShift:miss_refreg','The missing values at the boundaries of the reference region will be ignored')
      end
      intern = inter(1) + intern(1) - 1:inter(1) + intern(2) - 1;
      
   else
      [intern,flag_nan] = RemoveNaN(Ref(inter),cat(1,xT,xP),Select);
      intern            = vec(intern')';
   end
   if isempty(intern), error('iCoShift cannot handle this pattern of missing values.'), end
   if length(intern) ~= length(inter) && ~flag2
      
      if whole
         
         if length(intern) > 2
            
            if options(1) == 2
               xPBU = xP;
               xTBU = xT;
            end
            [XSeg,InOr] = ExtractSegments(cat(1,xT,xP),Ref(intern));
            inter       = [InOr(1),InOr(end)];
            InOr        = cat(2,Ref(intern),InOr);
            xP          = XSeg(2:end,:);
            xT          = XSeg(1,:);
            frag        = true;
            
         end
         
      else
         warning('iCoShift:Missing_values','\nTo handle the pattern of missing values, %i segments are created/removed',abs(length(intern) - length(inter))/2)
         inter       = intern;
         [nint,mint] = size(inter);
      end
      
   end
   
end

%% INITIALISE OUTPUT
xCS = xP;

%% ALIGN
[mi,pmi] = min(inter);
[ma,pma] = max(inter);
flag     = length(inter) > 1 && isequal([1,pma - pmi + 1],size(inter)) && isequal(unique(diff(inter,1,2)),1);
if flag % a Reference signal has been chosen
   
   if options(1)
      
      if  strcmpi(n,'b')      %     Option 'best search'
         fprintf('\n Automatic searching for the best "n" for the reference window "refW" enabled \n That can take a longer time \n')
      elseif  strcmpi(n,'f')  %     Option 'fast search'
         fprintf('\n Fast automatic searching for the best "n" for the reference window "refW" enabled \n')
      end
      
   end
   if max_flag                                  % MAX
         [amax,bmax] = max(sum(xP(:,mi:ma),2)); % MAX
         xT(mi:ma) = xP(bmax,mi:ma);            % MAX
   end
   ind                         = NaN(nP,1);
   missind                     = ~all(isnan(xP),2);
   [xCS(missind,:),ind(missind)] = coshifta(xT,xP(missind,:),inter,n,[options(1:3),BLOCKSIZE]);
   ints                        = [1 mi ma];
   
else % splits the dataset in intervals
   
   % Intervals definition
   if mint > 1
      if rem(mint,2), error('Wrong definition of intervals ("inter")'); end
      if ma>mP, error('Intervals ("inter") exceed samples matrix dimension'); end
      allint=[(1:round(mint/2))' inter(1:2:mint)' inter(2:2:mint)'];
   end
   %  Cheks interval overlapping
   sallint = sortrows(allint,2);
   sinter = reshape(sallint(:,2:3)',1,numel(sallint(:,2:3)));
   intdif = diff(sinter);
   if any(intdif(2:2:length(intdif)) < 0)
      uiwait(msgbox('The user-defined intervals are overlapping: is that intentional?','Warning','warn'));
   end
   clear sallint; clear sinter; clear intdif
   
   ints = allint;
   ind  = zeros(nP,size(allint,1));
   % end Interval definition
   
   if options(1)
      
      if  strcmpi(n,'b')
         fprintf('\n Automatic searching for the best "n" for each interval enabled \n That can take a longer time \n')
      elseif  strcmpi(n,'f')
         fprintf('\n Fast automatic searching for the best "n" for each interval enabled \n' )
      end
      
   end
   
   % Local co-shifting
   %     top_int = NaN(1,size(allint,1));
   for i = 1:size(allint,1)
      
      if options(1) ~= 0;
         if whole, fprintf('\n Co-shifting the whole %g samples... \r',nP);
         else fprintf('\n Co-shifting interval no. %g of %g... \r',i,size(allint,1));
         end
      end
      intervalnow = xP(:,allint(i,2):allint(i,3));
      if max_flag                               % MAX
         [amax,bmax] = max(sum(intervalnow,2)); % MAX
         target  = intervalnow(bmax,:);         % MAX
         xT(allint(i,2):allint(i,3)) = target;  % MAX
      else                                      % MAX
         target  = xT(:,allint(i,2):allint(i,3));
      end
      missind     = ~all(isnan(intervalnow),2);
      if ~all(isnan(target)) && sum(missind) ~= 0
         [cosh_interval,loc_ind]              = coshifta(target,intervalnow(missind,:),0,n,[options(1:3),BLOCKSIZE]);
         xCS(missind,allint(i,2):allint(i,3)) = cosh_interval;
         ind(missind,i)                       = loc_ind';
      else
         xCS(:,allint(i,2):allint(i,3)) = intervalnow;
      end
      %         top_int(i) = max(max(intervalnow));
      home
      
   end
   
end

% If n == 'whole' and some values are missing, reconstruct signal
if frag
   
   Xn  = NaN(nP,mP);
   for i_sam = 1:nP
      
      for i_seg = 1:size(InOr,1),
         
         Xn(i_sam,InOr(i_seg,1):InOr(i_seg,2)) = xCS(i_sam,InOr(i_seg,3):InOr(i_seg,4));
         if loc_ind(i_sam) < 0
            if flag_nan(i_seg,1,i_sam), Xn(i_sam,InOr(i_seg,1):InOr(i_seg,1) - loc_ind(i_sam)) = NaN; end
         elseif loc_ind(i_sam) > 0
            if flag_nan(i_seg,2,i_sam), Xn(i_sam,InOr(i_seg,2) - loc_ind(i_sam) + 1:InOr(i_seg,2)) = NaN; end
         end
         
      end
      
   end
   if options(1) == 2 % if plots are activated
      xP  = xPBU;
      xT  = xTBU;
   end
   xCS = Xn;
   
end
target = xT;
%% PLOT COMPARISON rigid shift
if flag_coshift %Plots the actual raw data
   %     xP = xPo;
   %     clear xPo
   sp1 = 3; sp2 = 1;    %FRSA Subplot type
else sp1 = 2; sp2 = 1;   %FRSA Subplot type
end

if options(1)==2
   
   if flag_allsegs, Colour = [1 1 1];
   else Colour = [1.0000    0.8314    0.8314];
   end
   pbu = get(0,'units');
   set(0,'Units','pixels');
   Scrsiz = get(0,'ScreenSize');
   set(0,'Units',pbu);
   ScrLength = Scrsiz(3);
   ScrHight  = Scrsiz(4);
   figpos    = [0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];
   figure('unit','pix','Position',figpos);
   st        = ceil(mP/(2^14));
   x         = 1:st:mP;
   CM        = ExtCM(size(xP,1),1 - MyCM);
   top       = max(xP(:));
   bot       = min(xP(:));
   sh        = NaN(1,2);
   
   sh(sp2)     = subplot(sp1,1,sp2,'nextplot','add','colorord',CM,'color','w');
   if flag_coshift
      ph        = plot(Scal(x),xPo(:,x)'); %line([mi ma; mi ma], [0 0; top top], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2)
   else ph        = plot(Scal(x),xP(:,x)');
   end
   ax        = axis;
   ax(4)     = ax(4)*1.1;
   if ax(3) > bot
      ax(3) = bot*2;
   end
   if round(bot) >= 0
      ax(3) = -top/15;
   end
   if flag
      plot(Scal([mi ma; mi ma]), ax(3:4),':r','LineWidth',2);
      text((Scal(mi) + Scal(ma))/2,max(max(xP(:,inter))*1.1),'Ref.','Color','k', 'Rotation', 90, 'clip','on');
   else
      
      %         xpts = [(allint(1:end,2))' (allint(1:end,3))'; (allint(1:end,2))' (allint(1:end,3))'];
      XX   = [(allint(1:end,2))'; (allint(1:end,2))';(allint(1:end,3))';(allint(1:end,3))'];
      YY   = [ax(3);ax(4);ax(4);ax(3)]*ones(1,size(allint,1));
      if ~whole
         h = fill(Scal(XX),YY,Colour);
         set(h,'facealpha',0.4,{'tag'},cellstr(num2str(allint(:,1))),'buttondown',@DispInfo,'edgecolor','b','linewidth',1)
         %          plot(xpts, ax(3:4),'b','LineWidth',1);
         text((Scal(allint(:,2)) + Scal(allint(:,3)))/2,ones(size(allint,1),1)*top,num2str(allint(:,1)), 'Rotation', 90,'clip','on');
      end
      
   end
   axis tight
   title('Raw data');
   delete(ph)
   if flag_coshift
      ph1 = plot(Scal(x),xPo(:,x)'); %line([mi ma; mi ma], [0 0; top top], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2)
      avrg = nanmean(xPo);
      clear xPo;
   else ph1 = plot(Scal(x),xP(:,x)');
   end
   set(ph1,{'tag'},cellstr(num2str((1:nP)','Sam. #%i')),'buttondown',@DispInfo)
   
   % 2nd plot
   if flag_coshift
      sp2 = sp2+1;
      sh(sp2) = subplot(sp1,1,sp2,'nextplot','add','colorord',CM,'color','w');
      ph = plot(Scal(x),xP(:,x)');%line([mi ma; mi ma], [0 0; top top], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 2)
      if flag
         plot(Scal([mi ma; mi ma]), ax(3:4),':g','LineWidth',2);
         text((Scal(mi) + Scal(ma))/2,max(max(xP(:,inter))*1.1),'Ref.','Color','k', 'Rotation', 90, 'clip','on');
         title(['Reference aligned data',title_suff]);
      else
         
         if ~whole
            
            h = fill(Scal(XX),YY,Colour);
            set(h,'facealpha',0.5,{'tag'},cellstr(num2str(allint(:,1))),'buttondown',@DispInfo,'edgecolor','b','linewidth',1)
            %          plot(xpts, ax(3:4),'b','LineWidth',1);
            if flag_coshift
               title('pre-coshifted data');
            else title(['Interval-aligned data',title_suff]);
            end
            text((Scal(allint(:,2)) + Scal(allint(:,3)))/2,ones(size(allint,1),1)*top,num2str(allint(:,1)), 'Rotation', 90,'clip','on');
            
         else title('Whole spectra aligned data');
         end
         
      end
      axis tight
      delete(ph)
      ph2 = plot(Scal(x),xP(:,x)');%line([mi ma; mi ma], [0 0; top top], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 2)
      set(ph2,{'tag'},cellstr(num2str((1:nP)','Sam. #%i')),'buttondown',@DispInfo)
      linkaxes(sh,'xy')
   end
   
   % 2nd or 3rd plot
   sp2 = sp2+1;
   sh(sp2) = subplot(sp1,1,sp2,'nextplot','add','colorord',CM,'color','w');
   ph = plot(Scal(x),xCS(:,x)');%line([mi ma; mi ma], [0 0; top top], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 2)
   if flag
      plot(Scal([mi ma; mi ma]), ax(3:4),':g','LineWidth',2);
      text((Scal(mi) + Scal(ma))/2,max(max(xP(:,inter))*1.1),'Ref.','Color','k', 'Rotation', 90, 'clip','on');
      title(['Reference aligned data',title_suff]);
   else
      
      if ~whole
         
         h = fill(Scal(XX),YY,Colour);
         set(h,'facealpha',0.5,{'tag'},cellstr(num2str(allint(:,1))),'buttondown',@DispInfo,'edgecolor','b','linewidth',1)
         %          plot(xpts, ax(3:4),'b','LineWidth',1);
         title(['Interval-aligned data',title_suff]);
         text((Scal(allint(:,2)) + Scal(allint(:,3)))/2,ones(size(allint,1),1)*top,num2str(allint(:,1)), 'Rotation', 90,'clip','on');
         
      else title('Whole spectra aligned data');
      end
      
   end
   axis tight
   delete(ph)
   ph3 = plot(Scal(x),xCS(:,x)');%line([mi ma; mi ma], [0 0; top top], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 2)
   set(ph3,{'tag'},cellstr(num2str((1:nP)','Sam. #%i')),'buttondown',@DispInfo)
   linkaxes(sh,'xy')
   
   if flag_coshift
      subplot(sp1,1,1)
      h = plot (Scal(x),avrg(x), 'Color','r', 'LineWidth', 4);
      legend (h, 'average');
      subplot(sp1,1,2)
      h = plot (Scal(x),xT(x), 'Color','r', 'LineWidth', 4);
      legend (h, note);
   else
      if any(strcmpi(note,{'average','median', 'max'})) % MAX
         subplot(sp1,1,1)
         h = plot (Scal(x),xT(x), 'Color','r', 'LineWidth', 4);
         legend (h, note);
         a = sum(abs(ind),2);
         b = find(a == min(a));
         result = ['Sample(s) No. [' num2str(b') '] is/are the most similar to the ' note ' one'];
         msgbox(result)
      else subplot(sp1,1,1)
         h = plot (Scal(x),xT(x), 'Color','r', 'LineWidth', 4);
         legend (h, 'Target');
      end
   end
   
   if ~flag_ScalDir, set(sh,'xdir','rev'), end % Reverts X axis direction if Scal is monotonically decreasing
   for i_sam = 1:nP
      if flag_coshift
         setappdata(ph1(i_sam),'Extra',ph3(i_sam))
         setappdata(ph2(i_sam),'Extra',ph1(i_sam))
         setappdata(ph3(i_sam),'Extra',ph2(i_sam))
      else
         setappdata(ph1(i_sam),'Extra',ph3(i_sam))
         setappdata(ph3(i_sam),'Extra',ph1(i_sam))
      end
      setappdata(ph1(i_sam),'flag_AvoidEndlessRecursion',1)
      setappdata(ph3(i_sam),'flag_AvoidEndlessRecursion',1)
      if flag_coshift
         setappdata(ph2(i_sam),'flag_AvoidEndlessRecursion',1)
      end
   end
   
end
if flag_coshift,
   ind = ind + wint * ones(1,size(ind,2)); end

function [xW,ind,R] = coshifta(xT,xP,refW,n,options)
% function [xW,ind,R] = coshifta(xT,xP,refW,n,options)
% Based on Correlation optimized shifting (COshift)
% Author:
% Frans van den Berg* (FvdB) email: fb@life.ku.dk
% Improved by:
% Francesco Savorani* (FS)   email: frsa@life.ku.dk
% Giorgio Tomasi*§    (GT)   email: gto@life.ku.dk
%
% * Department of Food Science Quality & Technology - Spectroscopy and Chemometrics group - www.models.life.ku.dk
% § Department of Basic Science and Environment - Soil and Environmental Chemistry group - www.igm.life.ku.dk
% Faculty of Life Sciences - University of Copenhagen - Denmark
%
% loosely based on : V.G. van Mispelaar et.al. 'Quantitative analysis of target components by comprehensive two-dimensional gas chromatography'
%                   J. Chromatogr. A 1019(2003)15-29

% HISTORY
% 1.00.00 2008 Jun 19 First documented version
% 2.00.00 2009 Nov 24 Blocking in case of large X (32MB blocks by default - options(4))


% ERRORS CHECK
if (isempty(refW)), refW = 0; end
if all(refW >= 0), rw = length(refW);
else rw = 1;
end
options_def = [ones(1,3) 2^25];
options(end + 1:length(options_def)) = options_def(length(options) + 1:end);
if options(2) == 1, Filling = -inf;
elseif options(2) == 0, Filling = NaN;
else error('options(2) must be 0 or 1');
end

if  strcmpi(xT,'average'), xT = nanmean(xP); end
[nT,mT] = size(xT);
[nP,mP] = size(xP);
[nR,mR] = size(refW);

if (mT ~= mP), error('Target "xT" and sample "xP" must be of compatible size (same vectors, same matrices or row + matrix of rows)'); end
if any(n <= 0),error('Shift(s) "n" must be larger than zero'); end
if (nR ~= 1), error('Reference windows "refW" must be either a single vector or 0'); end
if rw > 1 && (min(refW)<1) || (max(refW)>mT), error('Reference window "refW" must be a subset of xP'); end
if (nT ~= 1), error('Target "xT" must be a single row spectrum/chromatogram'); end

auto = 0;
if  strcmpi(n,'b')
   auto     = 1;                      % switch for the best automatic serach on
   if rw ~=1
      n        = max(fix(0.05 * mR),10); % change here the first searching point for the best "n"
      src_step = fix(mR * 0.05);         % Change here to define the searching step
   else
      n        = max(fix(0.05 * mP),10); % change here the first searching point for the best "n"
      src_step = fix(mP * 0.05);         % Change here to define the searching step
   end
   try_last = 0;                          %FRSA
elseif strcmpi(n,'f')
   
   auto = 1; % switch for the fast automatic serach on %FRSA modified from 1 to 0
   if rw ~=1
      n        = mR-1; %floor(mR/2);     %FRSA speed up fast option
      %         src_step = round(mR/2)-1;        %FRSA
   else
      n        = mP-1; %floor(mP/2);     %FRSA
      %         src_step = round(mP/2)-1;        %FRSA
   end
   try_last = 1;                          %FRSA speed up fast option
end
if (nT ~= 1), error('ERROR: Target "xT" must be a single row spectrum/chromatogram'); end
% End ERRORS CHECK

xW                             = NaN(nP,mP); % GT # Insert NaNs
ind                            = zeros(1,nP);
BLOCKSIZE                      = options(4);
nBlocks                        = whos('xP');
nBlocks                        = ceil(nBlocks.bytes/BLOCKSIZE);
SamxBlock                      = floor(nP/nBlocks);
indBlocks                      = SamxBlock(1,ones(1,nBlocks));
indBlocks(1:rem(nP,SamxBlock)) = SamxBlock + 1;
indBlocks                      = [0,cumsum(indBlocks)];
   
% try_last = 0;

if auto == 1
   
   while auto == 1 % Automatic search for the best "n" for each interval
      
      if Filling == -inf,    xtemp = [xP(:,ones(1,n)) xP xP(:,mP(1,ones(1,n)))];
      elseif isnan(Filling), xtemp = [NaN(nP,n),xP,NaN(nP,n)];
      end
      if rw == 1, refW = 1:mP; end
      ind = NaN(nP,1);
      R   = [];
      for i_block = 1:nBlocks
         block_indices = indBlocks(i_block) + 1:indBlocks(i_block + 1);
         [dummy,ind(block_indices),Ri] = ... %FFT Co-Shifting
            CC_FFTShift(xT(1,refW),xP(block_indices,refW),[-n n 2 1 Filling]); %#ok<ASGLU> 
            R = cat(1,R,Ri);
      end
      temp_index = -n:n;
      for i_sam = 1:nP
         index       = find(temp_index == ind(i_sam));
         xW(i_sam,:) = xtemp(i_sam,index:index + mP - 1);
      end
      if (max(abs(ind)) == n) && try_last ~= 1 %to check
         if n >= size(refW,2), break, end
         n = n+src_step ;
         continue
      elseif (max(abs(ind)) < n) && n+5 < size(refW,2) && try_last ~= 1
         n = n+src_step ;
         try_last = 1;
         continue
      else
         auto = 0;
         if options(1) ~= 0
            fprintf('\n Best shift allowed for this interval = %g \n',n);
         end
      end
      
   end
   
else % Uses the same defined "n" for each interval
   
   if Filling == -inf,    xtemp = [xP(:,ones(1,n)) xP xP(:,mP(1,ones(1,n)))];
   elseif isnan(Filling), xtemp = [NaN(nP,n),xP,NaN(nP,n)];
   end
   if rw == 1, refW = 1:mP; end
   ind = NaN(nP,1);
   R   = [];
   for i_block = 1:nBlocks
      block_indices = indBlocks(i_block) + 1:indBlocks(i_block + 1);
      [dummy,ind(block_indices),Ri] = ... %FFT Co-Shifting
         CC_FFTShift(xT(1,refW),xP(block_indices,refW),[-n n 2 1 Filling]); %#ok<ASGLU>
      R = cat(1,R,Ri);
   end
   temp_index = -n:n;
   for i_sam = 1:nP
      index       = find(temp_index == ind(i_sam));
      xW(i_sam,:) = xtemp(i_sam,index:index + mP - 1);
   end
   if (max(abs(ind)) == n) && options(1) ~= 0
      disp('Warning: Scrolling window size "n" may not be enough wide because extreme limit has been reached')
   end
   
end

function [inter] = defints(xP,interv,opt)
% Defines a vector "inter" [startint1 endint1 startint2 endint2...startintN endintN]
% for splitting a dataset "xP" into intervals having an imposed wideness of "interv" datapoints.
%
% SYNTAX
% [inter] = defints(xP,interv,opt)
%
% INPUT
% Xp (n x m): matrix to be split in regular intervals "interv" points wide
% interv (1 x 1): number of datapoints defining intervals' wideness
% opt (1 x 1): triggers warnings on/off
%
% OUTPUT
% inter [startint1 endint1 startint2 endint2...startintN endintN]: vector
% defining intervals boundaries
%
% Author: Francesco Savorani
%         frsa@life.ku.dk
% 11/11/08 (FrSa)
%
[nP,mP] = size(xP);
sizechk = mP/interv-round(mP/interv);
plus = (mP/interv-round(mP/interv))*interv;
mbx=['Warning: the last interval will not fulfill the selected intervals size "inter" = ', num2str(interv)];
if plus >= 0
   mbx2=['Size of the last interval = ', num2str(plus)];
else
   mbx2=['Size of the last interval = ', num2str(interv+plus)];
end
mbx3={mbx;mbx2};
if opt(1)==2 && (sizechk ~= 0)
   uiwait(msgbox(mbx3,'Warning','warn'));
elseif opt(1)~=0 && (sizechk ~= 0)
   fprintf('\r Warning: the last interval will not fulfill the selected intervals size "inter"=%g. \r Size of the last interval = %g \r',interv,plus);
end
t = cat(2,0:interv:mP,mP);
if t(end) == t(end - 1), t(end) = []; end
t     = cat(1,t(1:end - 1) + 1,t(2:end));
inter = t(:)';

function [Xwarp,Shift,Values] = CC_FFTShift(T,X,Options)
% Calculate optimal rigid shift using Fast Fourier Transform
% Signals are expected to be on the rows (horizontal slabs).
%
% SYNTAX
% [Xwarp,Shift,Values] = CC_FFTShift(T,X,Options)
%
% INPUT
% T      : target
% X      : signals
% Options: 1 x 4 vector. If element is NaN, the default value is used
%          (1) lower bound for shift
%          (2) upper bound for shift
%          (3) mode along which shift is to be determined (the last by default)
%          (4) 0 -> find maximum mean of cross-correlation for multivariate signals
%              1 -> find maximum product of cross-correlation for multivariate signals (default)
%          (5) fill in value (default: NaN), -inf uses the previous/following point
%
% OUTPUT
% Xwarp  : aligned signals
% Shifts : (i) shift for the i-th signal
% Values : (i + 1,j) cost function for the i-th signal and shift (j),
%          the first row contains the investigated shifts
%
% Author: Giorgio Tomasi
%         giorgio.tomasi@gmail.com
%
% Created      : August 2006
% Last modified: 24 March, 2009; 01:15

% HISTORY
% 1.00.00 August 2006    -> First working version
% 1.01.00 23 March, 2009 -> Fixed problem with single signals
% 2.00.00 23 March, 2009 -> Included handling of leading and trailing NaNs
% 2.00.01 23 Nov,   2009 -> Minor M-lint stuff

dimX = size(X);
dimT = size(T);
if nargin < 3, Options = []; end

% Default options, max shift is half the target's length
OptionsDefault                          = [-fix(dimT(end) * 0.5),fix(dimT(end) * 0.5) ndims(T) 1 NaN];
Options(end + 1:length(OptionsDefault)) = OptionsDefault(length(Options) + 1:length(OptionsDefault));
Options(isnan(Options))                 = OptionsDefault(isnan(Options));
if Options(1) > Options(2), error('Lower bound for shift is larger than upper bound'), end
TimeDim = Options(3);
if not(isequal(dimX([2:TimeDim - 1,TimeDim + 1:end]),dimX([2:TimeDim - 1,TimeDim + 1:end])))
   error('Target and signals do not have compatible dimensions')
end

% Normalise to length 1 along time dimension
ord        = [TimeDim,2:TimeDim - 1,TimeDim + 1:ndims(X),1]; % Leaves the sample mode as last
X_fft      = permute(X,ord);
X_fft      = reshape(X_fft,dimX(TimeDim),prod(dimX(ord(2:end))));
X_fft      = full(X_fft / sparse(1:prod(dimX(ord(2:end))),1:prod(dimX(ord(2:end))),sqrt(SumwNaN(X_fft.^2))));
T          = permute(T,ord);
T          = reshape(T,dimT(TimeDim),prod(dimT(ord(2:end))));
T          = Normalise(T);
[nP,mP]    = size(X_fft);
nT         = size(T,1);

% Remove leading NaN's
flag_miss = any(isnan(X_fft(:))) || any(isnan(T(:)));
if flag_miss
   
   if ndims(X) > 2, error('Multidimensional handling of missing not implemented, yet'), end
   MissOff = NaN(1,mP);
   for i_signal = 1:mP
      
      Limits = RemoveNaN([1,nP],X_fft(:,i_signal)',@all);
      if ~isequal(size(Limits),[1,2]), error('Missing values can be handled only if leading or trailing'); end
      if any(cat(2,Limits(1),mP - Limits(2)) > max(abs(Options(1:2)))),
         error('Missing values band larger than largest admitted shift'),
      end
      MissOff(i_signal) = Limits(1);
      if MissOff(i_signal) > 1, X_fft(1:Limits(2) - Limits(1) + 1,i_signal)  = X_fft(Limits(1):Limits(2),i_signal); end
      if Limits(2) < nP, X_fft(Limits(2) - Limits(1) + 1:nP,i_signal) = 0; end
      
   end
   Limits                            = RemoveNaN([1,nT],T',@all);
   T(1:Limits(2) - Limits(1) + 1,:)  = T(Limits(1):Limits(2),:);
   T(Limits(2) - Limits(1) + 1:nP,:) = 0;
   MissOff                           = MissOff(1:mP) - Limits(1);
   
end

% Zero pad to avoid pollution (cf. Press & Teukolsky pg. 540 and 545)
X_fft = cat(1,X_fft,zeros(max(abs(Options(1:2))),prod(dimX(ord(2:end)))));
T     = cat(1,T,zeros(max(abs(Options(1:2))),prod(dimT(ord(2:end)))));

% Calculate cross-correlation
len_fft = max(size(X_fft,1),size(T,1));
Shift   = Options(1):Options(2);
if Options(1) < 0 && Options(2) > 0
   ind = [len_fft + Options(1) + 1:len_fft,1:Options(2) + 1];
elseif Options(1) < 0 && Options(2) < 0
   ind = len_fft + Options(1) + 1:len_fft + Options(2) + 1;
elseif Options(1) < 0 && Options(2) == 0
   ind = [len_fft + Options(1) + 1:len_fft + Options(2),1];
else % Options(1) >= 0 && Options(2) > 0
   ind = Options(1) + 1:Options(2) + 1;
end
X_fft = fft(X_fft,len_fft,1);
T_fft = fft(T,len_fft,1);
T_fft = conj(T_fft);
T_fft = reshape(T_fft(:,:,ones(dimX(1),1)),size(X_fft));
cc    = ifft(X_fft .* T_fft,len_fft,1);
cc    = reshape(cc(ind,:),[Options(2) - Options(1) + 1,prod(dimX(ord(2:end - 1))),dimX(1)]);
if Options(4) == 0
   cc = squeeze(mean(cc,2))';
elseif Options(4) == 1
   cc = squeeze(prod(cc,2))';
else
   error('Invalid options for correlation of multivariate signals')
end
[nil,pos]  = max(cc,[],2);
Values     = cat(1,Shift,cc);
Shift      = Shift(pos);
if flag_miss, Shift = Shift + MissOff; end
Xwarp      = NaN([dimX(1),dimT(2:end)]);
[ind,indw] = deal(repmat({':'},ndims(X),1));
for i_X = 1:dimX(1)
   
   ind{1}  = i_X;
   indw{1} = i_X;
   if Shift(i_X) >= 0
      
      ind{TimeDim}  = Shift(i_X) + 1:dimX(TimeDim);
      indw{TimeDim} = 1:(dimX(TimeDim) - Shift(i_X));
      if Options(5) == -inf
         ind{TimeDim}  = cat(2,ind{TimeDim},dimX(TimeDim(ones(1,abs(Shift(i_X))))));
         indw{TimeDim} = cat(2,indw{TimeDim},dimX(TimeDim) - Shift(i_X) + 1:dimX(TimeDim));
      end
      
   elseif Shift(i_X) < 0
      
      ind{TimeDim}  = 1:(dimX(TimeDim) + Shift(i_X));
      indw{TimeDim} = -Shift(i_X) + 1:dimX(TimeDim);
      if Options(5) == -inf
         ind{TimeDim}  = cat(2,ones(1,-Shift(i_X)),ind{TimeDim});
         indw{TimeDim} = cat(2,1:-Shift(i_X),indw{TimeDim});
      end
      
   end
   Xwarp(indw{:}) = X(ind{:});
   
end

function DispInfo(varargin)
% Displays the 'tag' (or the 'userdata' if a cell vector of strings).
% The tag remains visible on the plot for 5 seconds (the time can be
% changed setting the DispInfoDelay property of the figure). It changes
% linewidth and color for the same amount of time as to highlight the object.
%
% DispInfo(Handle);
%
% Inputs:
% Handle: handle of the object whose tag is to be displayed (optional). Callback object is used if Handle
%         is not specified
%
% Outputs:
% None
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Created      : 05 March, 2003; 18:37
% Last modified: 09 March, 2009; 14:18

if nargin < 1
   h = gcbo;
else
   h = varargin{1};
end
s = get(h,'tag');
g = get(h,'userdata');
ex = getappdata(h,'Extra');
if ~isempty(g) && isa(g,'cell') && all(cellfun(@(x) isa(x,'char'),g))
   s = g;
end
a = get(gcf,'units');
set(gcf,'units','points')
t = get(gcf,'currentpoint');
set(gcf,'units',a)
a = axes('units','points','position',[t,50,10],'box','on','xtick',[],'ytick',[]);
b = text(.51,.5,s,'verticalalignment','middle','horizontalalignment','center','fontsize',9);
u = get(b,'units');
set(b,'units','points')
c = get(b,'extent');
set(b,'units',u);
set(a,'position',[t,c(3) * 1.05 c(4) * 1.2])
if isequal(get(h,'type'),'line')
   
   d = get(h,'linewidth');
   e = get(h,'color');
   if isequal(e,'r') || isequal(e,[1 0 0])
      col = [0 0 0.8];
   else
      col = 'r';
   end
   set(h,'linewidth',d+1,'color',col)
   
elseif isequal(get(h,'type'),'patch')
   
   d = get(h,'linewidth');
   e = get(h,'edgecolor');
   if isequal(e,'r') || isequal(e,[1 0 0])
      col = [0 0 0.8];
   else
      col = 'r';
   end
   set(h,'edgecolor',col,'linewidth',d + 1);
   
end
if isempty(getappdata(gcf,'DispInfoDelay'))
   Delay = 5;
else
   Delay = getappdata(gcf,'DispInfoDelay');
end
t = timer('StartDelay',Delay,'TimerFcn',@RemInfoAx,'userdata',{d e h,a},'stopfcn',@Delt);
start(t);
setappdata(h,'flag_AvoidEndlessRecursion',0)
if ~isempty(ex)
   a = getappdata(ex,'flag_AvoidEndlessRecursion');
   if isempty(a) || a, arrayfun(@DispInfo,ex); end
end

function RemInfoAx(varargin)
o = get(varargin{1},'userdata');
[d e h,a] = deal(o{:});
if isequal(get(h,'type'),'line')
   set(h,'linewidth',d,'color',e)
elseif isequal(get(h,'type'),'patch')
   set(h,'edgecolor',e,'linewidth',d);
end
setappdata(h,'flag_AvoidEndlessRecursion',1)
delete(a)

function Delt(varargin)
delete(varargin{1})

function [An,flag] = RemoveNaN(B, Signal, Select)
% Rearrange segments so that they do not include NaN's
%
% [Bn] = RemoveNaN(B, Signal, Select)
%
% INPUT
% B     : (p × 2) Boundary matrix (i.e. [Seg_start(1) Seg_end(1); Seg_start(2) Seg_end(2);...]
% Signal: (n × 2) Matrix of signals (with signals on rows)
% Select: (1 × 1) function handle to selecting operator
%                 e.g. @any (default) eliminate a column from signal matrix
%                                     if one or more elements are missing
%                      @all           eliminate a column from signal matrix
%                                     if all elements are missing
%
% OUTPUT
% Bn  : (q × 2)     new Boundary matrix in which NaN's are removed
% flag: (q × 2 × n) flag matrix if there are NaN before (column 1) or after (column 2)
%                   the corresponding segment in the n signals.
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Created      : 25 February, 2009
% Last modified: 23 March, 2009; 18:02

% HISTORY
% 1.00.00 09 Mar 09 -> First working version
% 2.00.00 23 Mar 09 -> Added output for adjacent NaN's in signals
% 2.01.00 23 Mar 09 -> Added Select input parameter

error(nargchk(2,3,nargin))
if nargin < 3, Select = @any; end
C      = NaN(size(B));
count  = 0;
Signal = isnan(Signal);
for i_el = 1:size(B,1)
   
   ind = B(i_el,1):B(i_el,2);
   in  = Select(Signal(:,ind),1);
   if any(in)
      
      p = diff([0 in],1,2);
      a = find(p < 0);
      b = find(p > 0) - 1;
      if ~in(1), a = cat(2,1,a); else b = b(2:end); end
      if ~in(end), b = cat(2,b,length(ind)); end
      a = unique(a);
      b = unique(b);
      C(count + 1:count + length(a),:) = ind(cat(2,a(:),b(:)));
      count = count + length(a);
      
   else
      count      = count + 1;
      C(count,:) = B(i_el,:);
   end
   
end
An = C;
if nargout > 1
   
   flag            = false(size(C,1),2,size(Signal,1));
   Cinds           = C(:,1) > 1;
   Cinde           = C(:,2) < size(Signal,2);
   flag(Cinds,1,:) = Signal(:,C(Cinds,1) - 1)';
   flag(Cinde,2,:) = Signal(:,C(Cinde,2) + 1)';
   
end

function Cols = ExtCM(n,CM)
% Extend colormap to have n colors
if nargin < 2
   CM = MyCM;
end
if isa(CM,'char')
   CM = eval(CM);
end
if isa(CM,'function_handle')
   CM = CM();
end
if ~isequal(size(CM),[64,3]) || any(CM(:) < 0 | CM(:) > 1)
   error('CM argument is not a colormap')
end
if n < 0 || ~isfinite(n)
   error('n must be positive and finite')
end
Cols           = interp1((1:64)',CM,linspace(1,64,n)');
Cols(Cols < 0) = 0;
Cols(Cols > 1) = 1;

function Colours = MyCM
% Purple to orange colour map
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
Colours(:,1) = linspace(0,1,64)';
Colours(:,2) = abs(linspace(-.5,.5,64)');
Colours(:,3) = flipud(Colours(:,1));

function [S] = SumwNaN(X, Dim, flag)
% Sum with missing values
%
% [S] = SumwNaN(X, Dim, flag)
% 
% INPUT
% X   : data array
% Dim : dimension along which the sum is to be performed
% flag: true if any missing / missingness of X
% 
% OUTPUT
% S: sum along Dim of X ignoring missing values
% 
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
% 
% Created      : 10 March, 2009; 15:59
% Last modified: 24 November, 2009; 16:03
if nargin < 1,error('not enough input arguments'); end
if nargin < 2,Dim = 1;end
if nargin < 3,flag = [];end
if isempty(flag) || (numel(flag) == 1 && flag)
   flag = isnan(X);
   N    = true;
elseif isequal(size(flag),size(X)), N = true;
elseif numel(flag) == 1 && ~flag,   N = false;
else error('Invalid flag'),
end
if N, X(flag) = 0; end
S = sum(X,Dim);

function [S] = MeanwNaN(X, Dim, flag)
% Sum with missing values
%
% [S] = MeanwNaN(X, Dim, flag)
% 
% INPUT
% X   : data array
% Dim : dimension along which the mean is calculated
% flag: true if any missing / missingness of X
% 
% OUTPUT
% S: sum along Dim of X ignoring missing values
% 
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
% 
% Created      : 24 November, 2009; 09:04
% Last modified: 24 November, 2009; 09:25

% HISTORY
% 1.00.00 2009 Nov 24 -> Dunction derived from SumwNaN
error(nargchk(1,3,nargin))
if nargin < 2,Dim = 1;end
if nargin < 3,flag = [];end
if isempty(flag) || (numel(flag) == 1 && flag)
   flag = isnan(X);
   N    = true;
elseif isequal(size(flag),size(X)), N = true;
elseif numel(flag) == 1 && ~flag,   N = false;
else error('Invalid flag'),
end
if N, 
   X(flag) = 0; 
   flag    = ~flag;
   c       = sum(flag,Dim);
   c(~c)   = NaN;
else
   c = size(X,Dim);
end
S = sum(X,Dim) ./ c;

function [Xn] = Normalise(X, Flag)
% Column-wise normalise matrix
% NaN's are ignored
%
% [Xn] = Normalise(X, Flag)
%
% INPUT
% X   : Marix
% Flag: true if any NaNs are present (optional - it saves time for large matrices)
%
% OUTPUT
% Xn: Column-wise normalised matrix
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Created      : 09 March, 2009; 13:18
% Last modified: 09 March, 2009; 13:50
if nargin < 2
   Patt = ~isnan(X);
   Flag = any(~Patt(:));
else
   if Flag, Patt = ~isnan(X); end
end
[M,N] = size(X);
Xn    = NaN(M,N);
if Flag
   
   for i_n = 1:N
      n = norm(X(Patt(:,i_n),i_n));
      if ~n, n = 1; end
      Xn(Patt(:,i_n),i_n) = X(Patt(:,i_n),i_n) / n;
   end
   
else
   
   for i_n = 1:N
      n = norm(X(:,i_n));
      if ~n, n = 1; end
      Xn(:,i_n) = X(:,i_n) / n;
   end
   
end

function [XSeg,SegNew] = ExtractSegments(X, Segments)
% Extract segments from signals
%
% [XSeg] = ExtractSegments(X, Segments)
%
% INPUT
% X       : (n × p) data matrix
% Segments: (s × 2) segment boundary matrix
%
% OUTPUT
% XSeg: (n × q) data matrix in which segments have been removed
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Created      : 23 March, 2009; 07:51
% Last modified: 23 March, 2009; 15:07

% HISTORY
% 0.00.01 23 Mar 09 -> Generated function with blank help
% 1.00.00 23 Mar 09 -> First working version

% Some initialisations
[n,p] = size(X);
Sd    = diff(Segments,1,2) + 1;
q     = sum(Sd + 1);
[s,t] = size(Segments);

% Check input
error(nargchk(2,2,nargin,'struct'))

flag_si = t ~= 2;
flag_in = any(Segments(:) ~= fix(Segments(:)));
flag_ob = any(Segments(:,1) < 1) || any(Segments(:,2) > p);
flag_ni = any(diff(Segments(:,1)) < 0) || any(diff(Segments(:,2)) < 0);
flag_ab = any(Sd < 2);
if flag_si, error('Segment boundary matrix must have two columns'), end
if flag_in, error('Segment boundaries must be integers'), end
if flag_ob, error('Segment boundaries outside of segment'), end
if flag_ni, error('Segments boundaries must be monotonically increasing'), end
if flag_ab, error('Segments must be at least two points long'), end

% Initialise output
if nargout > 0, XSeg = NaN(n,q); else return, end

% Calculate Segments
cdim = cat(1,0,cumsum(Sd));
for i_seg = 1:s
   XSeg(:,cdim(i_seg) + 1:cdim(i_seg + 1)) = X(:,Segments(i_seg,1):Segments(i_seg,2));
end
if nargout > 1, SegNew = cat(2,cdim((1:s)') + 1,cdim((2:s + 1)')); end

function pts = scal2pts(ppmi,ppm,prec)
% Transforms scalars in data points
%
% pts = scal2pts(values,scal)
%
% INPUT
% values: scalars whose position is sought
% scal  : vector scalars
% prec  : precision (optional) to handle endpoints
%
% OUTPUT
% pts   : position of the requested scalars (NaN if it is outside of 'scal')
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Created      : 12 February, 2009; 17:43
% Last modified: 11 March, 2009; 15:14

% HISTORY
% 1.00.00 12 Feb 09 -> First working version
% 1.01.00 11 Mar 09 -> Added input parameter check
if nargin < 3, prec = min(abs(unique(diff(ppm)))); end
dimppmi = size(ppmi);
ppmi    = ppmi(:);
ppm     = ppm(:);
rev     = ppm(1) > ppm(2);
if rev, ppm = ppm(end:-1:1); end
ubound  = (ppmi - ppm(end)) < prec & (ppmi - ppm(end)) > 0;
lbound  = (ppm(1) - ppmi) < prec & (ppm(1) - ppmi) > 0;
ppmi(ubound) = ppm(end);
ppmi(lbound) = ppm(1);
if nargin < 2, error('Not enough input arguments'); end
if length(ppmi) > length(ppm), warning('icoshift:scal2pts','ppm vector is shorter than the value''s'), end
[xxi,k]              = sort(ppmi(:));
[nil,j]                = sort([ppm(:);xxi(:)]);
r(j)                 = 1:length(j);
r                    = r(length(ppm) + 1:end) - (1:length(ppmi));
r(k)                 = r;
r(ppmi == ppm(end))  = length(ppm);
ind                  = find((r > 0) & (r <= length(ppm)));
ind                  = ind(:);
pts                  = Inf(size(ppmi));
pts(ind)             = r(ind);
ptsp1                = min(length(ppm),abs(pts + 1));
ptsm1                = max(1,abs(pts - 1));
ind                  = find(isfinite(pts));
dp0                  = abs(ppm(pts(ind)) - ppmi(ind));
dpp1                 = abs(ppm(ptsp1(ind)) - ppmi(ind));
dpm1                 = abs(ppm(ptsm1(ind)) - ppmi(ind));
pts(ind(dpp1 < dp0)) = pts(ind(dpp1 < dp0)) + 1;
pts(ind(dpm1 < dp0)) = pts(ind(dpm1 < dp0)) - 1;
if isempty(pts), pts = []; end
pts(not(isfinite(pts))) = NaN;
if rev, pts = length(ppm) - pts + 1; end
if ~isequal(size(pts),ppmi), pts = reshape(pts,dimppmi); end

function I = dscal2dpts(d,ppm,varargin)
% Translates an interval width from scal to the best approximation in sampling points.
%
% I = dppm2dpts(Delta,scal,prec)
%
% INPUT
% Delta: interval widths in scale units
% scal : scale
% prec : precision on the scal axes
%
% OUTPUT
% I: interval widths in sampling points
%
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
%
% Last modified: 21st February, 2009
%
if isempty(d), I = []; return; end
if nargin < 2, error('Not enough input arguments'); end
if d <= 0; error('Delta must be positive'); end
if ppm(1) < ppm(2), I = scal2pts(ppm(1) + d,ppm,varargin{:}) - 1;
else I = length(ppm) - scal2pts(ppm(end) + d,ppm,varargin{:}) + 1;
end