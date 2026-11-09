# source this: runtime for the sanitizer task (QGIS 4.2.3 built from this source tree)
export QGIS_VENV=$HOME/Desktop/VSCode/qgis_venv
export QGIS_BUILD=$HOME/Desktop/VSCode/QGIS/build
export PATH=$QGIS_VENV/bin:$PATH
export QT_QPA_PLATFORM=offscreen PROJ_DATA=$QGIS_VENV/share/proj
export PYTHONPATH=$QGIS_BUILD/output/python:$QGIS_VENV/share/qgis/python
export DYLD_FALLBACK_LIBRARY_PATH=$QGIS_BUILD/output/Contents/Frameworks:$QGIS_BUILD/output/lib:$QGIS_VENV/lib
