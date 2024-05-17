using MyProject
using RvSpectMLBase
using EchelleInstruments

path = "/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/Spectrum/Data/neidL1_20231014T152618.fits"

spectrum = NEID.read_data(path)

NEID.ChunkOfSpectrum(spectrum)


#TO BE CLEAR
#eclipse grass = src: plots (remake and finalize plots) + tests (rerun and get figures for tests)
#grass = figures (remake iag comparison for NEID)
#grass = scripts (update for array wavelength and turning for variability and rerun my scripts)
#grass = src: gpu (finalize and clear) + structures (GPU ones need to be finalized)
#grass once gpu done run with gpu + clear loose in src and loose outside src

