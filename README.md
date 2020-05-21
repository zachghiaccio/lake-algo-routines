# lake-algo-routines

This file provides a brief summary of each routine listed.

--------------------------------------------------------------------------------
## ATM
autoatm5.m - The "master" code for ATM lake depth retrieval. Runs through ILATM1B files one at a time.

atm_windowing_sub.m - The "windowing" subroutine for ATM. Breaks an ILATM1B track into discrete windows.

atm_chunking_sub.m - The "chunking" subroutine for ATM. Breaks an ILATM1B window into smaller chunks to resolve small lakes.

atm_sfc_detect_sub.m - Uses the chunked data to find lake surfaces and lake beds.

--------------------------------------------------------------------------------
## ATL03
autoatl.m - The "master" algorithm for IS2 lake depth retrieval. Currently only accepts one ATL03 file per code execution.

is2_class_merge.m - Condenses the "signal_conf_ph" parameter in ATL03 from [5 x N] to [1 x N], to remove dependence on surface type.

is2_windowing_sub.m - The "windowing" subroutine for ATL03. Breaks an ATL03 track into discrete windows.

is2_chunking_sub.m - The "chunking" subroutine for ATL03. Breaks an ATL03 window into smaller chunks to resolve smaller lakes.

is2_sfc_detect_sub.m - Uses the chunked data to find lake surfaces. An initial guess is also given for the lake bed.

is2_ph_dist.m - Uses photon refinement procedures (see Smith et al., 2019 for more details) to improve lake bed estimation. 

is2_polyfit_sub.m - Adds a 3rd-order polynomial curve to lake beds, to fill the gaps in attenuated regions. A separate subroutine was needed to account for multiple lakes in a scene.

manuatl.m - An alternative to "autoatl." If the coordinates of a lake are known, then this algorithm may be used for quicker analysis of singular lakes. Photon refinement is optional in this routine.
