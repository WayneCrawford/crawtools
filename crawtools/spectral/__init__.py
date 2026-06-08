"""
Spectral toolbox

Probably obsoleted by tiskitpy
"""
from .PSD import PSD, PSDs
# from .coherence import Coherence, Coherences
from .transfer_functions import TransferFunctions
from .spectral_density import SpectralDensity
from .data_cleaner import (DataCleaner, DCTFs, DCTF,
                           remove_str, strip_remove_str, strip_remove_one)
from .Peterson_noise_model import PetersonNoiseModel