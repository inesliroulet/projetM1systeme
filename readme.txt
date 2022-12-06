==========================
SAM FILE READER & ANALYSER
==========================

Prerequisites:

Python 3.7 (https://www.python.org/downloads/). Not compatible with Python 1 or 2.  
We did not test on Python releases anterior to 3.7 and cannot vouch for proper functioning on those releases.

Instructions:

1. Download the folder to the directory of your choice.

2. Check that the folder conntains both SAM_reader.py and flag_table.csv.

3. Add your SAM file to the folder.

4. Open the command prompt of your computer and go to the directory of the folder, using the "cd" command.

5. Give permission to execute the program by typing "chmod +x ./SAM_reader.py" and pressing 'Enter'.

6. To start the SAM reader, type the following command line: "./SAM_reader.py" and press 'Enter'.

The program will then ask for the name of the desired SAM file, that you previously put in the folder.
Type the name in, with its extension, and press 'Enter'. If the name is incorrect, the program will stop,
and you will have to launch it again.

After successfully identifying the file, the program will ask whether or not you want to change the default
minimum quality value. You must answer by "y" for yes, or "n" for no. If you answer incorrectly, the program
will ask the question again until you give an appropriate answer.

The program will then analyse your file and give you a summary of its data, including:
- the total number of reads
- the number of reads with a quality above the default or given value
- the number of matched reads
- the number of mapped reads
- the number of unmapped reads
- the number of properly mapped reads

The program will also create new SAM files containing only the reads that were extracted during the analysis.
The type of extracted reads (matched, mapped, etc.) will be indicated in the new file's name.
The program will also display the time taken by the analysis, in seconds.