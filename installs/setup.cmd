@ECHO OFF
call C:\ProgramData\Anaconda3\Scripts\activate.bat
pip install -r requirements.txt
::icn3dpy reqs
conda install -c conda-forge nodejs
jupyter labextension install jupyterlab_3dmol
::ipynb shortcuts
python3 -m pip install nbopen
python3 -m nbopen.install_win
::C:\Users\SUL\anaconda3