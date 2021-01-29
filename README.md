# seaice_benchmark 

stores files to reproduce comparison runs for Mehlmann et al 

directories 16000, 8000, 4000, 2000 (16km, 8km, 4km, 2km), output: snapshots every 30min, timePhase = 1800.
- runfe00: deltaT = 1800., SEAICEscaleSurfStress = .FALSE., SEAICEreplPressFac = 0., SEAICEadvScheme = 77, (default) (diag frequency = 1, strange results)
- runfe01: deltaT = 120., SEAICEadvScheme = 77, (default); this run is used for the paper
- runfe02: deltaT = 120., SEAICEadvScheme =  7, 
- runfe03: deltaT = 120., SEAICEadvScheme =  1, (only for 2000)
