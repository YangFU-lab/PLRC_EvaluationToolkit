%% Procedure 2: Spectral Combination for Effective Spectral Reflectance
% Run this script from the directory containing Procedure2_Input.m.
% INPUTS: r_DHlamp, r_solsim, r_pseudo, i_source, each [nm, value].
% Reflectances must already be dark/reference corrected, dimensionless,
% and compared at compatible instrumental spectral resolution.
% Optional third column for reflectance inputs: 1 = valid, 0 = rejected.
% Invalid rows are not interpolated across. No extrapolation is permitted.
% Quality/linearity/stability and negligible PL outside the declared
% emission band must be established experimentally, not inferred here.
clear; clc; close all;
run Procedure2_Input.m

%% 1. User-defined physical boundaries and acceptance limits
cfg.lambdaMin = 250;
cfg.lambdaMax = 2500;
cfg.emStart = 450;       % Replace with experimentally determined PL start
cfg.emEnd = 800;         % Replace with experimentally determined PL end
cfg.solarReliableMin = 350; 
cfg.solarReliableMax = 1000;
cfg.stepNm = 2;          % Sampling interval, NOT instrumental resolution
cfg.rmseLimit = 0.02;    % Absolute reflectance: 2 percentage point
cfg.seamLimit = 0.02;    % Proposed operational limit; validate experimentally
cfg.mode = 'auto';      % 'auto', 'direct', or 'mean'
% auto: try direct combination; if seams fail, try regional-mean combination.
% mean: regional-mean combination
%       DH over [lambdaMin, solarReliableMin);
%       DH + Solar mean over [solarReliableMin, emStart]; Solar over (emStart, emEnd);
%       DH + Solar + Pre mean over [emEnd, solarReliableMax];
%       Pre over (solarReliableMax, lambdaMax].

%% 2. Consistency check, seam assessment and spectral combination
out = combine_procedure2(r_DHlamp, r_solsim, r_pseudo, i_source, cfg);
if out.accepted
    fprintf('Consistency-check confirmed with width of %.1f nm\n', out.checkWidthNm);
else
    error('Procedure2:Rejected', '%s', out.message);
end
lambda = out.lambda;
R_comb = out.R_comb;
I_sc = out.I_sc;
ESR = out.ESR;

%% 3. Visualize calibrated inputs, merged spectrum, and junctions
figure('Name', 'Procedure 2: spectral combination');
plot(r_DHlamp(:,1), r_DHlamp(:,2), 'b--'); hold on;
plot(r_solsim(:,1), r_solsim(:,2), 'Color', [0.85 0.55 0]);
plot(r_pseudo(:,1), r_pseudo(:,2), 'k:');
plot(lambda, R_comb, 'r-', 'LineWidth', 1.4);
xline(cfg.emStart, '--'); xline(cfg.emEnd, '--');
if strcmp(out.method, 'mean')
    xline(cfg.solarReliableMin, ':');
    xline(cfg.solarReliableMax, ':');
end
xlabel('Wavelength (nm)'); ylabel('Effective spectral reflectance');
legend('DH lamp', 'Solar simulator', 'Pre-monochromator', ...
    'Combined spectrum', 'Location', 'best');
title('Combination of Spectral Reflectance'); grid on;

fprintf('ESR over %.1f-%.1f nm: %.6f\n', cfg.lambdaMin, cfg.lambdaMax, ESR);
% i_source defines the weighting spectrum. It is AM1.5 only if that is
% actually the supplied spectrum. Values of R_eff > 1 are retained.

function out = combine_procedure2(rDH, rSolar, rPre, iSource, cfg)
%COMBINE_PROCEDURE2 Validate and combine independently calibrated spectra.
% No scale factors are estimated from PL samples. A rejected consistency
% check requires diagnosis and, where appropriate, calibrated gray-standard
% measurements. Apply independently validated corrections before rerunning.
%
% Each input has two columns [wavelength_nm, value], or three columns with
% a validity flag (0/1). Nonfinite values or zero validity flags split the
% spectrum into separate valid runs. Interpolation never crosses such a run
% boundary. Missing intervals must be explicitly marked by invalid rows.
% Inputs must be resolution-compatible before calling this function.

limits = [cfg.lambdaMin cfg.solarReliableMin cfg.emStart cfg.emEnd ...
    cfg.solarReliableMax cfg.lambdaMax];
assert(all(isfinite(limits)) && all(diff(limits)>0), ...
    'Require lambdaMin < solarReliableMin < emStart < emEnd < solarReliableMax < lambdaMax.');
assert(isscalar(cfg.stepNm) && isfinite(cfg.stepNm) && cfg.stepNm>0);
assert(isscalar(cfg.rmseLimit) && isfinite(cfg.rmseLimit) && cfg.rmseLimit>=0);
assert(isscalar(cfg.seamLimit) && isfinite(cfg.seamLimit) && cfg.seamLimit>=0);
mode = validatestring(cfg.mode, {'auto','direct','mean'});
% Include exact boundaries even if not integer multiples of stepNm.
lambda = unique([(cfg.lambdaMin:cfg.stepNm:cfg.lambdaMax)'; limits(:)]);
R = [validInterp(rDH,lambda), validInterp(rSolar,lambda), ...
    validInterp(rPre,lambda)];
% Explicitly exclude unreliable solar measurements beyond 1000 nm (or the
% independently established limit provided by the user).
R(lambda>cfg.solarReliableMax,2) = NaN;

out = struct('accepted',false,'status','unchecked','message','', ...
    'method','','lambda',lambda,'R_inputs',R,'R_comb',[], ...
    'I_sc',[],'ESR',NaN,'checkWidthNm',0,'cfg',cfg);
out.consistency = table(); out.directSeams = table();
out.meanSeams = table(); out.selectedSeams = table();
out.spectrum = table();

%% A. Three pairwise RMSE values on the SAME common valid domain
mask = lambda>=cfg.emEnd & lambda<=cfg.solarReliableMax ...
    & all(isfinite(R),2);
out.checkMask = mask;
% Only adjacent pairs of valid nodes contribute: do not bridge invalid gaps.
edge = mask(1:end-1) & mask(2:end);
dl = diff(lambda);
width = sum(dl(edge));
out.checkWidthNm = width;
if width<=0
    out.status = 'no_common_interval';
    out.message = ['No common valid non-emissive interval is available. ' ...
        'Extend the measured coverage or revise the method; do not extrapolate.'];
    return
end
pairs = [1 2; 1 3; 2 3];
names = ["DH-Solar"; "DH-Pre"; "Solar-Pre"];
rmse = zeros(3,1);
for k=1:3
    d = R(:,pairs(k,1))-R(:,pairs(k,2));
    q = (d(1:end-1).^2+d(2:end).^2)/2;
    rmse(k) = sqrt(sum(dl(edge).*q(edge))/width);
end
out.consistency = table(names,rmse,rmse<=cfg.rmseLimit, ...
    'VariableNames',{'Pair','RMSE','Pass'});
if any(rmse>cfg.rmseLimit)
    out.status = 'consistency_failed';
    out.message = ['At least one pairwise RMSE exceeds the limit. Check ' ...
        'reference normalization, geometry, stability, signal quality and ' ...
        'residual PL. Use certified non-fluorescent gray standards to ' ...
        'diagnose bias; apply only independently validated corrections.'];
    return
end

%% B. Direct method seam mismatch at the SAME wavelength, not adjacent samples
boundaryR = [validInterp(rDH,limits(2:5)'), ...
    validInterp(rSolar,limits(2:5)'), validInterp(rPre,limits(2:5)')];
% Rows: solar reliable start, PL start, PL end, solar reliable end.
directLeft = [boundaryR(2,1); boundaryR(3,2)];
directRight = [boundaryR(2,2); boundaryR(3,3)];
out.directSeams = seamTable(limits(3:4)',directLeft,directRight,cfg.seamLimit);

%% C. Regional-mean method has FOUR boundaries; test each of them
meanLeft = [boundaryR(1,1); mean(boundaryR(2,1:2)); ...
    boundaryR(3,2); mean(boundaryR(4,:))];
meanRight = [mean(boundaryR(1,1:2)); boundaryR(2,2); ...
    mean(boundaryR(3,:)); boundaryR(4,3)];
out.meanSeams = seamTable(limits(2:5)',meanLeft,meanRight,cfg.seamLimit);
if strcmp(mode,'auto')
    if all(out.directSeams.Pass)
        method = 'direct';
    else
        method = 'mean';
    end
else
    method = mode;
end
out.method = method;
if strcmp(method,'direct')
    out.selectedSeams = out.directSeams;
else
    out.selectedSeams = out.meanSeams;
end
if ~all(out.selectedSeams.Pass)
    out.status = 'seam_failed';
    out.message = ['The selected combination has an unacceptable or ' ...
        'unmeasurable boundary mismatch. Equal averaging does not guarantee ' ...
        'continuity. Check calibration/coverage or revise the overlap windows.'];
    return
end

%% D. Source assignment: no extrapolation, scaling, clipping or smoothing
weights = zeros(numel(lambda),3);
if strcmp(method,'direct')
    % Direct method
    weights(lambda<cfg.emStart,1) = 1;
    weights(lambda>=cfg.emStart & lambda<=cfg.emEnd,2) = 1;
    weights(lambda>cfg.emEnd,3) = 1;
else
    % Regional-mean method
    weights(lambda<cfg.solarReliableMin,1) = 1;
    low = lambda>=cfg.solarReliableMin & lambda<=cfg.emStart;
    weights(low,1:2) = 1/2;
    weights(lambda>cfg.emStart & lambda<cfg.emEnd,2) = 1;
    high = lambda>=cfg.emEnd & lambda<=cfg.solarReliableMax;
    weights(high,:) = 1/3;
    weights(lambda>cfg.solarReliableMax,3) = 1;
end
out.weights = weights;
missing = any(weights>0 & ~isfinite(R),2);
out.missingRequiredWavelengths = lambda(missing);
if any(missing)
    out.status = 'coverage_failed';
    out.message = sprintf(['Required valid measured data are missing at ' ...
        '%d output wavelengths (first: %.3f nm). No extrapolation is allowed.'], ...
        nnz(missing),lambda(find(missing,1)));
    return
end
% Only zero-weight unavailable inputs may be discarded.
Rused = R; Rused(weights==0) = 0;
merged = sum(weights.*Rused,2);
I = validInterp(iSource,lambda);
if any(~isfinite(I)) || any(I<0) || trapz(lambda,I)<=0
    out.status = 'irradiance_failed';
    out.message = 'The weighting spectrum must be valid over the full integration range, nonnegative, and have positive integral.';
    return
end
out.R_comb = merged;
out.I_sc = I;
out.ESR = trapz(lambda,merged.*I)/trapz(lambda,I);
out.spectrum = table(lambda,merged,I,weights(:,1),weights(:,2),weights(:,3), ...
    'VariableNames',{'Wavelength_nm','R_eff','I_source','Weight_DH', ...
    'Weight_Solar','Weight_Pre'});
out.accepted = true; out.status = 'accepted';
out.message = 'Consistency, seam and coverage checks passed.';
end

function T = seamTable(wavelength,left,right,limit)
J = abs(left-right);
pass = isfinite(J) & J<=limit;
T = table(wavelength,left,right,J,pass, ...
    'VariableNames',{'Wavelength_nm','LeftValue','RightValue','Mismatch','Pass'});
end

function yq = validInterp(data,xq)
assert(isnumeric(data) && ismatrix(data) && size(data,1)>=2 ...
    && any(size(data,2)==[2 3]), 'Inputs must be N-by-2 or N-by-3 numeric arrays.');
assert(all(isfinite(data(:,1))), 'Wavelength values must be finite.');
[x,order] = sort(data(:,1)); data = data(order,:);
assert(all(diff(x)>0), 'Wavelengths must be unique.');
valid = isfinite(data(:,2));
if size(data,2)==3
    assert(all(ismember(data(:,3),[0 1])), 'Validity flags must be 0 or 1.');
    valid = valid & data(:,3)==1;
end
start = find(diff([false;valid])==1);
stop = find(diff([valid;false])==-1);
yq = NaN(size(xq));
for k=1:numel(start)
    ind = start(k):stop(k);
    if numel(ind)==1
        yq(xq==x(ind)) = data(ind,2);
    else
        select = xq>=x(ind(1)) & xq<=x(ind(end));
        yq(select) = interp1(x(ind),data(ind,2),xq(select),'linear');
    end
end
end
