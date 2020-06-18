# lake-algo-routines

A series of MATLAB routines designed to perform supraglacial lake depth retrievals with two laser altimeters: the Airborne Topographic Mapper (ATM) and the Ice, Cloud, and Land Elevation Satellite-2 (ICESat-2).

MATLAB 2016a or newer is required, otherwise the routines will not run.

--------------------------------------------------------------------------------
## ATM
REQUIRED DATA - IceBridge ATM L1B Elevation and Return Strength (ILATM1B). The narrow-swath equivalent (ILNSA1B) is also valid. The ILATM1B "Q-Fit" product does NOT work for this routine.

autoatm5.m - The "master" code for ATM lake depth retrieval. Currently accepts a list (in .txt format) of ILATM1B files, but it may be configured to accept one file.

atm_windowing_sub.m - The "windowing" subroutine for ATM. Breaks an ILATM1B track into discrete windows.

atm_chunking_sub.m - The "chunking" subroutine for ATM. Breaks an ILATM1B window into smaller chunks to resolve small lakes.

atm_sfc_detect_sub.m - Uses the chunked data to find lake surfaces and lake beds.

--------------------------------------------------------------------------------
## ATL03
REQUIRED DATA - ATL03 Geolocated Photon Data. The routine accepts both Version 002 and Version 003 data. 

autoatl.m - The "master" algorithm for IS2 lake depth retrieval. Currently only accepts one ATL03 file per code execution.

is2_class_merge.m - Condenses the "signal_conf_ph" parameter in ATL03 from [5 x N] to [1 x N], to remove dependence on surface type.

is2_windowing_sub.m - The "windowing" subroutine for ATL03. Breaks an ATL03 track into discrete windows.

is2_chunking_sub.m - The "chunking" subroutine for ATL03. Breaks an ATL03 window into smaller chunks to resolve smaller lakes.

is2_sfc_detect_sub.m - Uses the chunked data to find lake surfaces. An initial guess is also given for the lake bed.

is2_sfc_find.m - Gives a second guess for the location of a lake surface. Currently has issues with smaller lakes seen on Greenland.

is2_ph_dist.m - Uses photon refinement procedures (see Smith et al., 2019 for more details) to improve lake bed estimation. 

is2_polyfit_sub.m - Adds a 3rd-order polynomial curve to lake beds, to fill the gaps in attenuated regions. A separate subroutine was needed to account for multiple lakes in a scene.

depth_refrac_fix.m - Corrects lake depth retrievals for the refractive index of water. See Parrish et al., (2019) for more details.

manuatl.m - An alternative to "autoatl." If the coordinates of a lake are known, then this algorithm may be used for quicker analysis of singular lakes. Photon refinement is optional in this routine. Recommended if the lat/lon coordinates are known for a supraglacial lake.

--------------------------------------------------------------------------------
## Known Issues

The "is2_sfc_find.m" routine is designed to refine a lake surface assignment in a manner resembling ATL06 surface heights. However, it struggles with lakes smaller than 500 m in diameter and interferes with the lake bed estimation. Therefore, it is not recommended to apply surface refinement at this time when running "manuatl.m". The "autoatl.m" script currently does not incorporate this routine.

Very small (diameter < 200 m) and/or shallow (depth < 1 m) lakes may be difficult to resolve, especially with ICESat-2. In such cases, the photon cloud is often too noisy for the algorithm to notice the presence of a lake. 


--------------------------------------------------------------------------------
## Contact and Release Information

Have questions about the code? Contact Zachary Fair at zhfair@umich.edu.

The most recent version of the code was completed in October 2019. I have since moved on to other projects, so any future updates will be limited in scope. However, this project was envisioned as a proof-of-concept for lake depth retrievals, from which other researchers can refine and improve the technique. Therefore, users are welcomed to make changes and improvements to the code, provided they use the citation given below. 


--------------------------------------------------------------------------------
## Citation

If using either algorithm, please cite using the following study:

Fair Z., Flanner, M., Brunt, K. M., Fricker, H. A., and Gardner, A. S.: Using ICESat-2 and Operation IceBridge altimetry for supraglacial lake depth retrievals, The Cryosphere Discuss.,  https://doi.org/10.5194/tc-2020-136, in review, 2020.
