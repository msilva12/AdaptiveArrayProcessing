%% Optimal Costrained Beamformer
%%
% <html><h3>Array and Source Parameters</h3></html>
%
% * Number of array elements: $$L=16$
% * Element spacing: $$d=\frac{\lambda}{2}m$
% * Propagation speed: $$c=3*10^8\frac{m}{s}$
% * Signal to Noise Ratio: $$SNR=[3, 10]dB$
% * Signal power: $$ps=1$
% * Noise power: $$\sigma_n^2=[0.5, 0.1]$
% * Number of samples: $$N=[10,100]$
% * Source direction: $$\theta_k=30^o$
% * Steering Vector in the Look direction: $$\textbf{S}_0=[exp({j2{\pi}f_0{\tau}_l}),\ldots,exp({j2{\pi}f_0{\tau}_L})]^T$
%
% Where $$\tau_l(\theta_k)=\frac{d}{c}(l-1)sin(\theta_k)$
%
% <html><h3>Optimal Costrained Beamformer</h3></html>
%
% Weight vector: $$\hat{\textbf{w}}=\frac{R^{-1}_N\textbf{S}_0}{\textbf{S}_0^HR^{-1}_N\textbf{S}_0}$
%
% <html><h3>Correlation Matrix</h3></html>
%
% $$R = p_s\textbf{S}_0\textbf{S}_0^H+p_I\textbf{S}_I\textbf{S}_I^H+\sigma^2_nI$
%
% <html><h3>Array output</h3></html>
%
% $$y(t)=\hat{\textbf{w}}^H\textbf{x}(t)$
%
% <html><h3>Simulation Algorithm</h3></html>
%
% # Generate signal for L element array
% # Use Coventional Beamformer weight vector
% # Calculate the array output y
% # Calculate output SNR
% # Calculate array Gain
%% Signal + Noise Only
%
% <html><h3>Array input</h3></html>
%
% $$\textbf{x}(t)=exp({j2{\pi}f_0t})\textbf{S}_0+\textbf{n}(t)$$
%
% <html><h3>Correlation Matrix</h3></html>
%
% $$R=ps*\textbf{S}_0\textbf{S}_0^H+{\sigma}_n^2{I}$$
%
% <html><h3>Signal to Noise Ratio</h3></html>
%
% $$SNR=\frac{ps}{\sigma_n^2}L$$
%
% <html><h3>Array Gain</h3></html>
%
% $$G=L$$
%
% <html><h3>Simulations</h3></html>
%
optmbeam_snr
%% References
% [1] L. Godara, Smart antennas. Boca Raton: CRC Press, 2004.