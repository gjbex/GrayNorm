GrayNorm
========

What is it?
-----------
GrayNorm is an algorithm that helps the researcher to identify suitable reference genes that are minimally influenced by the experiment.  Gene expression levels are determined by quantitative reverse transcription PCR (RT-qPCR) and one is interested in up- or down regulation of various genes as a consequence of varying experimental conditions.

Requirements
-----------
GrayNorm has been implemented in Python 2.7.x, using only Python's standard library. Hence it will run as a command line program on any platform on which Python is installed. In the following its usage in the “Enthought Canopy” Python platform (free version available from https://store.enthought.com/) is outlined.

How to use it?
--------------
GrayNorm (as the “graynorm.txt” or the `graynorm.py` file) can be opened in Canopy in the editor window (top window). Then the Graynorm procedure on the data can be called from the Canopy command line as follows:
```
%run graynorm.py -in data.csv -out output.csv
```
OR 
```
%run graynorm.txt -in data.csv -out output.csv
```
The input file (in the above command line e.g., `data.csv`) is a file in CSV format, and it will be documented in the next section. The output of the program is stored in an CSV file as well (in the above example e.g., `output.csv` in the example above). The GrayNorm script (`graynorm.py`) and the input file should be in the same working directory that can be chosen as follows: choose “change working directory” from the drop down menu in the top right corner of the lower command window in Canopy. 

Input data format
-----------------
The input data is in CSV format, and can be prepared in a spreadsheet application such as Microsoft Excel if desired.  The file consists of two parts, the header and the body. The header contains meta-information, each line starting with a `#` character, while the body contains the actual experimental data, arranged in rows and columns. Each row corresponds to an individual sample and contains a sample ID for identification, the values of the variables for the sample, and the gene expression levels for each of the genes measured in the experiment.

The header should contain at least the following information:

1) the name of the column containing the sample ID, e.g.,
```
 # sampleid: samplenr
```
2) the names of the columns that contain the expression levels of the reference genes, separated by commas, e.g.,
```
 # refgenes: Msd1, At5g15710, At5g08290, At2g28390, ACT2
```
3) the list of variables, and their control values, separated by commas, e.g.,
```
 # controls: exposure = 0, time = 0, genotype = wildtype
```

Additional header lines (starting with `#`) are ignored.

The body starts with a line containing the column names, and hence should at least have the sample ID, the variables names, as well as the names of the reference genes as field values.

Each additional line contains the experimental data for an individual sample. Hence the body will have as many lines as there are samples obtained from the experiment, plus one header line.

Blanks lines are ignored.

Refer to the Supplemental File 7 that contains three examples of data input files. Make sure that there is correspondence in writing (e.g., capital letters) between header and body information.

Output format
-------------
The output is an CSV file, suitable for further analysis in a spreadsheet application such as Microsoft Excel. it starts with a row containing the column titles, further rows list the metrics for each set of candidate normalization genes. In Excel, a text-to-columns transformation using `,` as a column separator would format the data in an appropriate table.

Contact information
-------------------
For questions on the GrayNorm algorithm in general, please contact:
Tony Remans <tony.remans@uhasselt.be>

For questions regarding the actual implementation, for support or to report bugs, please contact:
Geert Jan Bex <geertjan.bex@uhasselt.be>

Copyright
---------
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
