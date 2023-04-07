# VevoRawReader
A MATLAB Class to Read and Organize Raw Data from Vevo LAZR-X (2100 and 3100) and Convert to NIfTI Files

## Quick Setup
Clone the repository (or download the VevoRawReader.m file) and store in a directory on your computer. Add that directory to MATLAB path interactively via Home/Environment/Set Path.  
Then create an instance of the class whenever you want to run a method (function) from this class by:
```
VevoRawReader_instance = VevoRawReader;
```

## Organizing Files
### Read Header Info from .xml File
Simply pass the filename without the extension to the following command.
```
file_attributes = VevoRawReader_instance.read_vevo_xml(filename);
```
### Table Containing Header Information
To create a MATLAB Table containing essential header information from individual scans in a given directory, run the following code on the directory_of_interest location:
```
scan_info = VevoRawReader_instance.scan_information_table(directory_of_interest);
```

## Reading Files into NIfTI Format
Use either of the `vevo_nifti_write_bmode`, `vevo_nifti_write_oxyhemo`, or `vevo_nifti_write_color` methods to read corresponding raw files and convert them to NIfTI. 
For example:
```
VevoRawReader_instance.vevo_nifti_write_bmode(fnameBase, scanmode, output_filename);
```
where:


**fnameBase** is the name of the file without extension.  
**scanmode** is an arbitrary description that goes into the name of the output (e.g. the bmode image could be pure_bmode, or bmode_color, or bmode_oxyhemo, depending on what type of other image type it accompanies)  
**output_filename** is the desired output nifti file name without extension.

## There is more!
Refer to the code for more information on other useful methods in the class.
