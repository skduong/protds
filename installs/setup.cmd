@echo off
call C:\Anaconda3\Scripts\activate.bat
pip install -r requirements.txt
::icn3dpy reqs
conda install -c conda-forge nodejs
jupyter labextension install jupyterlab_3dmol
::pdb query
pip install git
pip install git+git://github.com/williamgilpin/pypdb