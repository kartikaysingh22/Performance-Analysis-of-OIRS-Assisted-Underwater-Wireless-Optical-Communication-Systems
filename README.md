# Performance-Analysis-of-OIRS-Assisted-Underwater-Wireless-Optical-Communication-Systems
This repository contains the MATLAB-based simulation framework and research findings for a comparative study between Mirror Array-based OIRS and Metasurface-based OIRS in Underwater Wireless Optical Communication (UWOC) systems.
Optical Intelligent Reflecting Surfaces (OIRS) are a promising solution for Non-Line-of-Sight (NLOS) underwater links, providing high-bandwidth, low-latency alternatives to traditional acoustic and RF communication

# Key Features
1. Detailed Channel Modeling: Per-element modeling of incidence, irradiance, and reflection characteristics.
2. Environmental Realism: Includes exponential underwater attenuation and Exponential-Generalized Gamma (EGG) turbulence models.
3. Performance Metrics: Evaluates Spectral Efficiency (SE), Signal-to-Noise Ratio (SNR), and Bit Error Rate (BER).
4. Architectural Comparison: Contrast between dynamic per-element steering (Mirror Array) and fixed-geometry steering (Metasurface)

# System Model & Parameters
The simulation assumes a transmitter and receiver separated by 3 meters at a depth of 1 meter below an OIRS plane.
- OIRS Grid Size                13×13 (169 elements) 
- Element Size                  0.08×0.08 m 
- Field of View (FOV)           90∘ 
- Water Attenuation (cλ​)        0.056 m−1 
- LED Lambertian Order          m=2 
- Modulation                    On-Off Keying (OOK) 
- Turbulence Model              Multi-EGG (Weak, Moderate, Strong)
- Photodiode Responsivity (Re​)  0.5 
- Noise Variance (σi2​)          10^(−12) 


# Key Findings
Mirror Array Dominance: Consistently outperforms Metasurfaces due to its ability to dynamically align every micro-mirror element to satisfy the law of reflection.
Spectral Efficiency: Mirror Arrays achieved a peak SE of approximately $30~bits/s/Hz$, significantly higher than Metasurfaces and flat mirrors.
Turbulence Resilience: Mirror Arrays maintain better BER performance in turbulent water because their beams are more concentrated.
Alignment Sensitivity: SE drops rapidly with receiver misalignment; however, the Mirror Array decays slower due to superior steering precision
