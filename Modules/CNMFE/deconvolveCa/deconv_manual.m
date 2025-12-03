function [C, S, C_denoised, kernel_pars, noise_sn] = deconv_manual(C_raw, deconv_options, use_parallel, method_noise)
% Deconvolve manually extracted traces using CNMF-E's solver
%
% Inputs:
%   C_raw          – K x T raw traces (or neuropil-corrected traces)
%   deconv_options – structure with CNMF-E options.deconv_options
%   use_parallel   – true/false (default true)
%   method_noise   – 'psd' (default) or 'histogram'
%
% Outputs:
%   C              – denoised (calcium) trace, K x T
%   S              – sparse inferred spikes, K x T
%   C_denoised     – C_raw with estimated baseline removed (like CNMF-E’s obj.C_raw)
%   kernel_pars    – AR(1)/AR(2) time constants for each ROI
%   noise_sn       – noise ? for each ROI
%
% Author: rewritten from CNMF-E Sources2D.deconvTemporal for manual ROI mode

if nargin < 3 || isempty(use_parallel)
    use_parallel = true;
end
if nargin < 4 || isempty(method_noise)
    method_noise = 'psd';
end

[K, T] = size(C_raw);

if K == 0
    warning('deconv_manual: No neurons provided for deconvolution.');
    C = [];
    S = [];
    C_denoised = [];
    kernel_pars = [];
    noise_sn = [];
    return
end

% Prepare containers
C = zeros(K,T);
S = zeros(K,T);
C_denoised = zeros(K,T);
kernel_pars = zeros(K,3);     % typical AR(1)/AR(2) kernels
noise_sn = zeros(K,1);

C_cells = mat2cell(C_raw, ones(K,1), T);

fprintf('Deconvolving %d traces...\n', K);
num_per_row = 80;
for m = 1:K, fprintf('|'); if mod(m,num_per_row)==0, fprintf('\n'); end, end
fprintf('\n');

% ----- Parallel version -----
if use_parallel
    parfor k = 1:K
        ck_raw = C_cells{k};

        % Handle NaNs
        if any(isnan(ck_raw))
            C(k,:) = zeros(1,T);
            S(k,:) = zeros(1,T);
            C_denoised(k,:) = zeros(1,T);
            noise_sn(k) = 0;
            kernel_pars(k,:) = [0 0 0];
            continue
        end

        % Estimate noise
        if strcmpi(method_noise, 'histogram')
            [~, sn_k] = estimate_baseline_noise(ck_raw);
        else
            sn_k = GetSn(ck_raw);
        end
        noise_sn(k) = sn_k;

        % Deconvolution
        [ck, sk, tmp] = deconvolveCa(ck_raw, deconv_options, 'sn', sn_k);

        if sum(abs(ck)) == 0
            ck = ck_raw;
        end

        % Save
        C(k,:) = ck;
        S(k,:) = sk;
        C_denoised(k,:) = ck_raw - tmp.b;
        kernel_pars(k,:) = tmp.pars(:)';
    end

% ----- Serial version -----
else
    for k = 1:K
        ck_raw = C_cells{k};

        if any(isnan(ck_raw))
            C(k,:) = zeros(1,T);
            S(k,:) = zeros(1,T);
            C_denoised(k,:) = zeros(1,T);
            noise_sn(k) = 0;
            kernel_pars(k,:) = [0 0 0];
            continue
        end

        if strcmpi(method_noise, 'histogram')
            [~, sn_k] = estimate_baseline_noise(ck_raw);
        else
            sn_k = GetSn(ck_raw);
        end
        noise_sn(k) = sn_k;

        [ck, sk, tmp] = deconvolveCa(ck_raw, deconv_options, 'sn', sn_k);

        if sum(abs(ck)) == 0
            ck = ck_raw;
        end

        C(k,:) = ck;
        S(k,:) = sk;
        C_denoised(k,:) = ck_raw - tmp.b;
        kernel_pars(k,:) = tmp.pars(:)';

        fprintf('.');
        if mod(k,num_per_row)==0, fprintf('\n'); end
    end
end

fprintf('\nDone.\n');
end


