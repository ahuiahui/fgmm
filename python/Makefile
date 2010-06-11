

_fgmm.so : fgmm_wrap.o 
	g++ -shared fgmm_wrap.o ../_build_/default/libfgmm.a  -o _fgmm.so

fgmm_wrap.o : fgmm_wrap.cxx 
	 g++ -c -fpic -I/usr/include/python2.6  -I../src fgmm_wrap.cxx

fgmm_wrap.cxx : fgmm.i 
	swig -c++ -python -classic fgmm.i