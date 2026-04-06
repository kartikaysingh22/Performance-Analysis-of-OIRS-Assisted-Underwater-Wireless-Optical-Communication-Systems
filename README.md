# Performance-Analysis-of-OIRS-Assisted-Underwater-Wireless-Optical-Communication-Systems
This repository contains the MATLAB-based simulation framework and research findings for a comparative study between Mirror Array-based OIRS and Metasurface-based OIRS in Underwater Wireless Optical Communication (UWOC) systems.
Optical Intelligent Reflecting Surfaces (OIRS) are a promising solution for Non-Line-of-Sight (NLOS) underwater links, providing high-bandwidth, low-latency alternatives to traditional acoustic and RF communication

# Key Features
Detailed Channel Modeling: Per-element modeling of incidence, irradiance, and reflection characteristics.
Environmental Realism: Includes exponential underwater attenuation and Exponential-Generalized Gamma (EGG) turbulence models.
Performance Metrics: Evaluates Spectral Efficiency (SE), Signal-to-Noise Ratio (SNR), and Bit Error Rate (BER).
Architectural Comparison: Contrast between dynamic per-element steering (Mirror Array) and fixed-geometry steering (Metasurface)

# System Model & Parameters
The simulation assumes a transmitter and receiver separated by 3 meters at a depth of 1 meter below an OIRS plane.
Parameter,Value
OIRS Grid Size,13×13 (169 elements) 
Field of View (FOV),90∘ 
Water Attenuation (cλ​),0.056 m−1 
Modulation,On-Off Keying (OOK) 
Turbulence Model,"Multi-EGG (Weak, Moderate, Strong) +1"

