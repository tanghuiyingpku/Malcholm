function dens = CalcFLDens(i,P)
global Fluid
refDens = Fluid.fluid{i}.dens;
refP = Fluid.fluid{i}.refp;
cmpr = Fluid.fluid{i}.compr*1e6;
dens = refDens*(1+cmpr*(P-refP));
end