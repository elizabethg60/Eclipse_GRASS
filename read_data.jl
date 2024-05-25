using MyProject
using RvSpectMLBase
using EchelleInstruments

path = "/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Data/neidL1_20231014T152618.fits"

spectrum = NEID.read_data(path)

NEID.ChunkOfSpectrum(spectrum)
