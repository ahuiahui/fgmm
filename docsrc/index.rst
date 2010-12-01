.. fgmm documentation master file, created by
   sphinx-quickstart on Wed Dec  1 16:26:50 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to fgmm's documentation!
================================

fgmm aims at being a really fast Gaussian Mixture model implementation. It 
has virtually no depency at all, and is pretty fast now. 
fgmm is written in :doc:`C <rawc>`, but as wrappers for :doc:`c++ <cpp>`,
:doc:`python <python>` and matlab/octave.


Compilation
===========

the build system used by fgmm is `waf <http://code.google.com/p/waf/>`_ an
awesome building system written in `python <http://www.python.org>`_ .  You
don't need to install it, since it is provided with the source code,
but you do need python, so go install it, and become a really cool
computer user.  If you don't want to, write your own build script,
this should not be a big deal anyway.

compiling is as simple as::

  ./waf configure (-h shows you compilation options ) 
  ./waf 
  ./waf test
  ./waf install 


Contents:
=========
.. toctree::
   :maxdepth: 2


   c api  <rawc>
   python	
   c++ <cpp>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

