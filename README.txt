***********************************************************
*****************                        ******************
**************        Open Chemera           **************
**************           Library             **************
*****************                        ******************
***********************************************************
       Open source library for Chemera and BiGGER
                     Version: 0.15


***********************************************************
******                   DISCLAIMER                  ******

This software is provided "as is", with no explicit or 
implied warranties. None whatsoever. I mean it. I would 
even write all this in caps if only it didn't look so ugly.

The author and contributors disclaim any express or implied
warranties, including, but not limited to, the implied 
warranties of merchantability and fitness for a particular 
purpose. In no event shall the author or contributors be
liable for any direct, indirect, incidental, special, 
exemplary, or consequential damages (including, but not 
limited to, procurement of substitute goods or services; 
loss of use, data, or profits; or business interruption)
however caused and on any theory of liability, whether in 
contract, strict liability, or tort (including negligence 
or otherwise) arising in any way out of the use of this 
software, even if advised of the possibility of such damage.

In other words, use this software at your own risk. 

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
******                  Contents                     ******

bigger:
  The implementation of the BiGGER algorithm

cclassif:
  Contacts classifier
  
chemera
  visualization
  
dockprep
  Utility for preparing docking order files for bigger
  
pdbtools
  Utility for processing pdb structures and sequences
  
cclassif
  Classifier for contact predictions

  
***********************************************************
******               Acknowledgements                ******

This work was partially supported by Portuguese National 
funds through Fundação para a Ciência e Tecnologia (FCT) 
under project CREMA PTDC/EIA-CCO/115999/2009

***********************************************************

Source code available at:
https://github.com/lkrippahl/Open-Chemera

Ludwig Krippahl