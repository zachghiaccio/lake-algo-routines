# lake-algo-routines

This file provides a brief summary of each routine listed.

## ATL06
atl06_chunking_sub.m - Divides windowed ATL06 parameters into smaller chunks.

atl06_windowing_sub.m - Divides ATL06 parameters into iterative windows.

--------------------------------------------------------------------------------
## ATM
autoatm5.m - The "master" algorithm for ATM lake depth retrieval.

atm_chunking_sub.m - The "chunking" subroutine for ATM.

atm_windowing_sub.m - The "windowing" subroutine for ATM.

atm_sfc_detect_sub.m - Uses the chunked data to find lake surfaces and lake beds.

--------------------------------------------------------------------------------
## ATL03
autoatl.m - The "master" algorithm for IS2 lake depth retrieval. Also implements ATL06 data for quick geophysical analysis.

is2_class_merge.m - Condenses the "signal_conf_ph" parameter in ATL03 from 5 x N to 1 x N, to remove dependence on surface type.

is2_chunking_sub.m - The "chunking" subroutine for ATL03.

is2_windowing_sub.m - The "windowing" subroutine for ATL03.

is2_sfc_detect_sub.m - Uses the chunked data to find lake surfaces. Unlike ATM, lake beds are not found in this routine.

is2_ph_dist.m - Uses the windowed data to find lake beds. Based on statistical tools used in the ATL06 surface-finding algorithm.

is2_polyfit_sub.m - Adds a 3rd-order polynomial curve to lake beds, to fill the gaps in attenuated regions. A separate subroutine was needed to account for multiple lakes in a scene.

is2_raw_plotting.m - Quickly plots raw ATL03 data for easy visualization. Currently unused.

is2_sig_count_sub.m - Uses ATL03 photon counts to approximate signal-to-noise ratios. Currently unused in favor of ATL06 parameters.
