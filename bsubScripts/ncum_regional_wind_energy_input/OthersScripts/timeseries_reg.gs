'reinit'
'sdfopen final.nc'
'set grads off'
ifrm = 1
while(ifrm <= 25)
*'clear'
*'reset'
ret=read('latlon.csv')
    code=sublin(ret,2)
    stn=subwrd(code,1)
    yvl=subwrd(code,2)
    xvl=subwrd(code,3)
'set t 1'
'q time'
s1 = sublin(result,1)
s2 = subwrd(s1,3)
s5 = substr(s2,4,9)
s3 = substr(s2,1,2)
if (s3 =13)
s4 = 12Z
else
s4 = 00Z
endif
outfile='wind50m_'%stn%'_IC_'%s4%''%s5%'.csv'
'!rm -f 'outfile
dummy=write(outfile,'Date/Time ws(m/s) wd(deg)',append)
tt=72
'set lat '%yvl
'set lon '%xvl

n=1
while(n <= tt)
'set t 'n''
r1 = sublin(result,1)
r2 = subwrd(r1,5)

'set gxout print'
'set prnopts %13.2f 10 1'
'define ws = mag(ugrd_50mabovegr,vgrd_50mabovegr)'
'define wd = 57.3*atan2(ugrd_50mabovegr,vgrd_50mabovegr)+180'
'd ws'
l1 = sublin(result,2)
wss = subwrd(l1,1)
'd wd'
l1 = sublin(result,2)
wdd = subwrd(l1,1)

dummy=write(outfile,''%r2%' 'wss' 'wdd,append)
n = n +1
endwhile
ifrm = ifrm + 1
endwhile
cls=close('latlon.csv')

'quit'

