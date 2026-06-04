function [snr_fft_dB, snr_fft_lin, crlb_m_fft, crlb_m_persample] = snr_after_fft(snr_persample_dB, L_prs, BW_prs, c)
% SNR_AFTER_FFT  Post-coherent-integration SNR and CRLB (paper §II-B, eq. 7).
%
%   Cross-correlating the received signal against the L_prs-sample PRS
%   reference provides a coherent integration gain of L_prs, lifting the
%   effective SNR well above the per-sample link-budget value:
%
%       SNR_fft = L_prs × SNR_per_sample
%
%   The CRLB on range error expressed via SNR_fft is numerically identical
%   to the per-sample form that keeps L_prs inside the radical:
%
%       sigma_CRLB = c*sqrt(3) / (pi * BW * sqrt(2 * SNR_fft))
%                 = c*sqrt(3) / (pi * BW * sqrt(2 * L_prs * SNR_sample))
%
%   Both forms are returned so callers can verify |crlb_m_fft - crlb_m_persample| ~ 0.
%
%   Note on the paper's gamma factor: in eq.(10) SNR = gamma * SNR_inj where
%   gamma = M / (Delta_f * T_int) is the PRS duty cycle inside the DLL
%   integration window T_int.  In this simulation T_int equals exactly one
%   PRS burst, so gamma = 1 and SNR_eff = SNR_inj.  The coherent integration
%   gain (L_prs) captures the same physics for the matched-filter receiver.
%
%   Inputs
%     snr_persample_dB   per-sample SNR [dB], scalar or vector
%     L_prs              matched-filter (PRS waveform) length [samples]
%     BW_prs             active PRS bandwidth [Hz]
%     c                  speed of light [m/s]
%
%   Outputs
%     snr_fft_dB         post-FFT SNR [dB]  = snr_persample_dB + 10*log10(L_prs)
%     snr_fft_lin        post-FFT SNR [linear]
%     crlb_m_fft         CRLB range sigma computed from SNR_fft [m]
%     crlb_m_persample   CRLB range sigma computed via per-sample + L_prs [m]

snr_persample_lin = 10.^(snr_persample_dB ./ 10);

snr_fft_lin = L_prs .* snr_persample_lin;
snr_fft_dB  = snr_persample_dB + 10*log10(L_prs);

crlb_m_fft       = c * sqrt(3) ./ (pi * BW_prs * sqrt(2 .* snr_fft_lin));
crlb_m_persample = c * sqrt(3) ./ (pi * BW_prs * sqrt(2 .* L_prs .* snr_persample_lin));
end
