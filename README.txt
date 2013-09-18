***********************************************************
*****************                        ******************
**************        Open Chemera           **************
**************           Library             **************
*****************                        ******************
***********************************************************
       Open source library for Chemera and BiGGER
                     Version: 0.12


***********************************************************
******                   DISCLAIMER                  ******

This software is provided "as is", with no explicit or 
implied warranties. None whatsoever. I mean it. I would 
even write all this in caps if only it didn't look so ugly.

Use this software at your own risk.

***********************************************************
******                   Copyright                   ******

This software is public domain. You have the right to copy, 
distribute, reuse, modify, improve and debug it. Especially
the last one.

The *.pdb files in the data/ligands folder are from the
PDBeChem ligand dictionary. For more information, see

  http://www.ebi.ac.uk/msd-srv/chempdb 

***********************************************************
******                   Compiling                   ******

Using the Lazarus IDE:

1-  Install Lazarus and FPC, which are available here:
	http://www.lazarus.freepascal.org/

2-  Install the TOpenGLControl LCL control in the package
    lazarus/components/opengl/lazopenglcontext.lpk
	(in the Lazarus IDE, select "Package->Install/uninstall
	packages")
	Note for Debian distributions: these distros do not 
	include this package precompiled. Please compile 
	Lazarus from source on Debian systems.

3-  Open the project file (e.g. chemera/chemera.lpi)

4-  Press F9
		
***********************************************************
******                   Installing                  ******

Copy the files and subfolders in the oclibrary/data folder
to your configuration folder.

Windows 7
  [user]\AppData\Local\oclibrary

Ubuntu
  ~/.oclibrary

***********************************************************
******               Acknowledgements                ******

This work was partially supported by Portuguese National 
funds through Fundação para a Ciência e Tecnologia (FCT) 
under project CREMA PTDC/EIA-CCO/115999/2009

***********************************************************

Source code available at:
https://github.com/lkrippahl/Open-Chemera

Ludwig Krippahl
