# FAQ September 2020

### Is Graynorm Python 3 compatible?

Release 1.2 is Python 3 compatible.

### I'm stuck with Python 2, can I still use Graynorm?

Yes, you can, use release 1.1


## Trouble shooting

Please be aware of following essentials that are often overlooked causing failure in generating an output. 

* all header lines of the inputfile should start with "#", also the
    non-essential ones
* make sure there is exact correspondence in spelling and capitalization
    between the essential info in the # header lines and the headers of
    the data table (e.g., gene names)  
* make sure you select the correct working directory in the top right
    corner of the command window, and that the graynorm.py file and input
    file are in that folder. 
* do not run the script by trying to click the green arrow from the editor
    window (top), run the script by pressing the return key after entering
    the command line in the command window.
* make sure that Microsoft Excel is set to use "." and not "," as decimal
    separator (you can change this in the "file>options>advanced" menu
    of Excel).
