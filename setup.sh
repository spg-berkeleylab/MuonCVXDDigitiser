trackperf=/global/cfs/cdirs/atlas/sferrar2/TrackMuC/MuonCVXDDigitiser

export PATH=${trackperf}/install/bin:$PATH
export LD_LIBRARY_PATH=${trackperf}/install/lib:${trackperf}/install/lib64:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=${trackperf}/install/include:$ROOT_INCLUDE_PATH
export PYTHONPATH=${trackperf}/install/python:$PYTHONPATH
export CMAKE_PREFIX_PATH=${trackperf}/install:$CMAKE_PREFIX_PATH
