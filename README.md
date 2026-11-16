Directional selection can drive the evolution of modularity in complex traits
=============================================================================

[![Join the chat at https://gitter.im/lem-usp/evomod](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/lem-usp/evomod?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Codes and scripts for Melo, D., & Marroig, G. (2014). Directional selection can drive the evolution of modularity in complex traits. Proceedings of the National Academy of Sciences of the United States of America, 2015 112 (2) 470-475; doi:10.1073/pnas.1322632112



http://www.pnas.org/content/early/2014/12/25/1322632112.abstract

Requirements
============

1. GSL http://www.gnu.org/software/gsl/
2. IPython notebook http://ipython.org/notebook.html
3. ggplot2 http://ggplot2.org/

Python version
==============

`py/pop.py` is a Python 3 port of the C model. Its environment is managed with
[pixi](https://pixi.sh):

    cd py
    pixi install
    pixi run sim -n 2500 -t 1000 -d 0.01 -r 42   # see `pixi run sim -h`

Output goes to `py/dats/`.
