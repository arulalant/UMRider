# Change Log

## [v1.0.1](https://github.com/NCMRWF/UMRider/tree/v1.0.1) (2016-03-16)
[Full Changelog](https://github.com/NCMRWF/UMRider/compare/v1.0.0...v1.0.1)

**Implemented enhancements:**

- Set Centre Code from user defined config file [\#46](https://github.com/NCMRWF/UMRider/issues/46)
- Need to create date timestamp directories for loggin [\#35](https://github.com/NCMRWF/UMRider/issues/35)
- Fix unknown standard variable names in iris load [\#26](https://github.com/NCMRWF/UMRider/issues/26)
- we should include horizontal resolution in grib2 file name [\#24](https://github.com/NCMRWF/UMRider/issues/24)
- In out file name should give option for ensemble  [\#22](https://github.com/NCMRWF/UMRider/issues/22)
- Calculate surface upward longwave & softwave radiation flux [\#18](https://github.com/NCMRWF/UMRider/issues/18)
- Need to support for HYCOM  [\#17](https://github.com/NCMRWF/UMRider/issues/17)
- user specified ordered vars [\#15](https://github.com/NCMRWF/UMRider/issues/15)
- Need to support for IND REGION [\#14](https://github.com/NCMRWF/UMRider/issues/14)
- Need to support for OSF  [\#13](https://github.com/NCMRWF/UMRider/issues/13)
- Need to support for VSDB Input [\#12](https://github.com/NCMRWF/UMRider/issues/12)
- Need to support for custom grib2 table parameter options [\#69](https://github.com/NCMRWF/UMRider/issues/69)
- Make startdate, enddate as environmental variable  [\#63](https://github.com/NCMRWF/UMRider/issues/63)
- Need to support for IMD MFI Input [\#54](https://github.com/NCMRWF/UMRider/issues/54)
- support for call back script [\#53](https://github.com/NCMRWF/UMRider/issues/53)
- support for tar ball  [\#41](https://github.com/NCMRWF/UMRider/issues/41)
- Need to extract 3 hourly data and write into 3 hourly files [\#16](https://github.com/NCMRWF/UMRider/issues/16)
- create tar.bz2 in parallel using pbzip2 [\#61](https://github.com/NCMRWF/UMRider/pull/61) ([arulalant](https://github.com/arulalant))

**Fixed bugs:**

- Fix g2ctl.pl option & forecast reference, bounds time for analysis files [\#39](https://github.com/NCMRWF/UMRider/issues/39)
- sea ice concentration has full of zero [\#33](https://github.com/NCMRWF/UMRider/issues/33)

**Closed issues:**

- custom date formate in outfile name [\#11](https://github.com/NCMRWF/UMRider/issues/11)
- Grib2 to Grib1 conversion option [\#8](https://github.com/NCMRWF/UMRider/issues/8)
- support for custom out file names [\#7](https://github.com/NCMRWF/UMRider/issues/7)

**Merged pull requests:**

- update in docstring [\#67](https://github.com/NCMRWF/UMRider/pull/67) ([arulalant](https://github.com/arulalant))
- cleanup [\#66](https://github.com/NCMRWF/UMRider/pull/66) ([arulalant](https://github.com/arulalant))
- 	bug fixes; update python path; update doc [\#65](https://github.com/NCMRWF/UMRider/pull/65) ([arulalant](https://github.com/arulalant))
- support for UMRIDER\_STARTDATE & UMRIDER\_ENDDATE environmental vars [\#64](https://github.com/NCMRWF/UMRider/pull/64) ([arulalant](https://github.com/arulalant))
- update on IMD MFI grib2 with -set\_radius as 0 [\#62](https://github.com/NCMRWF/UMRider/pull/62) ([arulalant](https://github.com/arulalant))
- out, tmp, log, paths update [\#60](https://github.com/NCMRWF/UMRider/pull/60) ([arulalant](https://github.com/arulalant))
- update on osf, hycom path [\#59](https://github.com/NCMRWF/UMRider/pull/59) ([arulalant](https://github.com/arulalant))
- 	added pressureLevels option [\#58](https://github.com/NCMRWF/UMRider/pull/58) ([arulalant](https://github.com/arulalant))
- added tar ball creation callbackscript [\#57](https://github.com/NCMRWF/UMRider/pull/57) ([arulalant](https://github.com/arulalant))
- Fixed bug \#16 [\#56](https://github.com/NCMRWF/UMRider/pull/56) ([arulalant](https://github.com/arulalant))
- 	added callBackScript option & imd mfi support [\#55](https://github.com/NCMRWF/UMRider/pull/55) ([arulalant](https://github.com/arulalant))
- Updated doc [\#51](https://github.com/NCMRWF/UMRider/pull/51) ([arulalant](https://github.com/arulalant))
- Update on README and docs [\#50](https://github.com/NCMRWF/UMRider/pull/50) ([arulalant](https://github.com/arulalant))
- Updated doc [\#42](https://github.com/NCMRWF/UMRider/pull/42) ([arulalant](https://github.com/arulalant))
- support to extract 3-hourly analysis data [\#40](https://github.com/NCMRWF/UMRider/pull/40) ([arulalant](https://github.com/arulalant))
- changed to ultra [\#38](https://github.com/NCMRWF/UMRider/pull/38) ([arulalant](https://github.com/arulalant))
- added	temporary sourcing umtid\_bashrc [\#37](https://github.com/NCMRWF/UMRider/pull/37) ([arulalant](https://github.com/arulalant))
- create log date directory [\#36](https://github.com/NCMRWF/UMRider/pull/36) ([arulalant](https://github.com/arulalant))
- Fixed bug sea ice missing value \#33 [\#34](https://github.com/NCMRWF/UMRider/pull/34) ([arulalant](https://github.com/arulalant))
- Support for HYCOM Input \#17 [\#32](https://github.com/NCMRWF/UMRider/pull/32) ([arulalant](https://github.com/arulalant))
- kept landsea mask at last [\#31](https://github.com/NCMRWF/UMRider/pull/31) ([arulalant](https://github.com/arulalant))
- Included grid resolution in filename [\#30](https://github.com/NCMRWF/UMRider/pull/30) ([arulalant](https://github.com/arulalant))
- Included grid resolution in filename [\#29](https://github.com/NCMRWF/UMRider/pull/29) ([arulalant](https://github.com/arulalant))
- Supports for OSF Input [\#28](https://github.com/NCMRWF/UMRider/pull/28) ([arulalant](https://github.com/arulalant))
- added ncum load rules before loading in iris cubes via callback  [\#27](https://github.com/NCMRWF/UMRider/pull/27) ([arulalant](https://github.com/arulalant))
- Supports for VSDB Input and INDIAN REGION [\#19](https://github.com/NCMRWF/UMRider/pull/19) ([arulalant](https://github.com/arulalant))
- user defined file names; convert to grib1 option; grib2ctl option; g2ctl option [\#9](https://github.com/NCMRWF/UMRider/pull/9) ([arulalant](https://github.com/arulalant))
- Supports for 12 UTC 5 days long forecast [\#1](https://github.com/NCMRWF/UMRider/pull/1) ([arulalant](https://github.com/arulalant))
- fixed bugs in callback scripts of osf, hycom [\#72](https://github.com/NCMRWF/UMRider/pull/72) ([arulalant](https://github.com/arulalant))
- Fixed bugs in bsub scripts [\#71](https://github.com/NCMRWF/UMRider/pull/71) ([arulalant](https://github.com/arulalant))
- added setGrib2TableParameters option in setup.cfg file [\#70](https://github.com/NCMRWF/UMRider/pull/70) ([arulalant](https://github.com/arulalant))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*