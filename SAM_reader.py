#!/usr/bin/python3
#-*- coding : utf-8 -*-


__authors__ = ("Inès Liroulet", "Daphné Navratil")
__contact__ = ("ines.liroulet@etu.umontpellier.fr","daphne.navratil@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/06/2022"
__licence__ = "This program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\nGNU General Public License for more details.\nYou should have received a copy of the GNU General Public License\nalong with this program. If not, see <https://www.gnu.org/licenses/>."


############### IMPORTING MODULES ###############


import os, sys, re, time


############### DEFINING FUNCTIONS ###############


# 1 - Checking if the file is adequate

def check():
    print("\n") #Skipping a line
    fichier = input("Please provide a SAM file to analyse: ") #Asking for file
    if os.path.isfile(fichier): #Checking if path is a file
        if fichier.split(".")[1] == "sam": #Checking if extension of the file is sam
            if os.stat(fichier).st_size > 0: #Checking if size is greater than 0 bytes
                return fichier #All conditions are verified : file name is returned
            else:
                print("Your file is empty! Please launch the program again.") #Printing error message
                print("\n") #Skipping a line
                exit() #Exiting program
        else:
            print("Your file is not a SAM file! Please launch the program again.") #Printing error message
            print("\n") #Skipping a line
            exit() #Exiting program
    else:
        print("What you provided is not a file! Please launch the program again.") #Printing error message
        print("\n") #Skipping a line
        exit() #Exiting program

# 2 - Reading and storing file content

def readAndStore(fichier): #Parameter: given sam file
    f = open(fichier, "r") #Opening sam file in reading mode
    header = [] #Creating list to store header lines
    data = {} #Creating dictionary to store reads
    for ligne in f:
        if ligne[0] == "@": #Focusing on header lines
            header.append(ligne) #Adding each line to the list
        if ligne[0] != "@": #Excluding header lines
            chaine = ligne.split("\t") #Splitting the string (denominator = tabulation)
            key = chaine[0] #Getting read name to use as key
            if (key[-2:] == "/1") or (key[-2:] == "/2") or (key[-2:] == ".1") or (key[-2:] == ".2") or (key[-2:] == "_1") or (key[-2:] == "_2"): #For new SAM files : removing numbers from name of read + its mate
                key = key[:len(key)-2] #Removing last two caracters
            if key in data.keys(): #Checking if read key already in dictionary
                data[key].append(chaine[1:]) #If yes: adding other read (mate) to the key
            else :
                data[key] = [chaine[1:]] #If not: creating key and its corresponding read
                #List indices: 0 is flag, 1 is ref, 2 is position, 3 is quality, 4 is CIGAR, 5 is "=", 6 is mate reference, 7 is mate position, 8 is sequence length, 9 is quality sequence, 10+ is tags.
    f.close() #Closing file
    return header, data #Returning list of header lines and dictionary of all reads

# 3 - Analysing file content

#Asking minimum quality value to user

def question():
    print("\n") #Skipping a line
    qualityValue = 20 #Defining default value
    answer = input("The default minimum quality value is 20. Do you want to change it? (y/n) ").lower() #Asking and storing answer
    if answer == "y": #If user wants to change value
        qualityValue = int(input("Please provide a minimum quality value: ")) #Asking wanted value and storing it
        print("The minimum quality value will be set to:", qualityValue) #Confirming change to user
    elif answer == "n": #If not
        print("The minimum quality value will not be changed.") #Confirming answer to user
    else: 
        print("Please answer by typing 'y' or 'n'!") #Telling user to answer correctly if they don't
        question() #Asking question again (launching function again)
    return qualityValue #Returning minimum quality value

#Extracting reads with a minimal quality value

def extractQualityReads(data, minValue): #Parameters: (main) dictionary, miminum quality value
    dataQuality = {} #Creating dictionary to store wanted reads
    for read in data.items(): #Going through (main) dictionary
        for i in range(len(read[1])): #Going through each read stored in key
            if int(read[1][i][3]) >= minValue: #Checking if quality value above minimum value
                if read[0] in dataQuality.keys(): #If yes: checking if read key already in dictionary
                   dataQuality[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                else:
                     dataQuality[read[0]] = [read[1][i]] #If not: creating key and its corresponding read
    return dataQuality #Returning new dictionary with only wanted reads

#Extracting reads matching reference sequence (using CIGAR)

def testCigar(data): #Parameter: (main) dictionary
    dataCigar = {} #Creating dictionary to store wanted reads
    for read in data.items(): #Going through (main) dictionary
        for i in range(len(read[1])): #Going through each read stored in key
            if re.match("[0-9]+M", read[1][i][4]): #Checking if CIGAR is a series of digits followed by 'M'
                if read[0] in dataCigar.keys(): #If yes: checking if read key already in dictionary
                    dataCigar[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                else:
                    dataCigar[read[0]] = [read[1][i]] #If not: creating key and its corresponding read
    return dataCigar #Returning new dictionary with only wanted reads

#Checking if a flag value is in given flag

def flagSearch(flag, valFlag): #Parameters: searched flag value, given flag
    fichier_flag = open("flag_table.csv", "r") #Opening csv file with table of flag values in reading mode
    t = fichier_flag.readlines() #Reading each line and assigning them to a variable
    b = False #Setting boolean to False by default
    i = len(t)-1 #Setting i to start from last line of table
    while i >= 0 and b == False: #As long as i has not reached 0 and boolean is False
        value = int(t[i].split(",")[0]) #Select value in table line
        if value <= valFlag: #Checking if flag value inferior to given flag
            if value == flag: #If yes: checking if value is equal to searched flag
                b = True #If yes: boolean changes to True
            else:
                valFlag = valFlag - value #If no: substracting value from given flag
        i = i-1 #Decrease i by one
    return b #Return boolean value

#Extracting mapped reads and unmapped reads

def mapped_unmapped(data): #Parameter: (main) dictionary
    dataMapped = {} #Creating dictionary to store wanted reads (mapped)
    dataUnmapped = {} #Creating dictionary to store wanted reads (unmapped)
    listGood = [] #Creating list to store flags corresponding to mapped reads ("good list")
    listBad = [] #Creating list to store flags corresponding to unmapped reads ("bad list")
    for read in data.items(): #Going through (main) dictionary
        for i in range(len(read[1])): #Going through each read stored in key
            if (read[1][i][0] in listBad) == False: #If read flag not in bad list
                if read[1][i][0] in listGood: #If not in bad list: checking if read flag in good list
                    if read[0] in dataMapped.keys(): #If in good list: checking if read key already in mapped dictionary
                        dataMapped[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                    else:
                        dataMapped[read[0]] = [read[1][i]] #If not: creating key and its corresponding read
                elif flagSearch(4,int(read[1][i][0])) == False: #If not in good list: checking if read is mapped
                    listGood.append(read[1][i][0]) #If yes: adding read flag to good list
                    if read[0] in dataMapped.keys(): #And checking if read key already in mapped dictionary
                        dataMapped[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                    else:
                        dataMapped[read[0]] = [read[1][i]]  #If not: creating key and its corresponding read
                else: #If read flag not in bad list and not mapped
                    listBad.append(read[1][i][0]) #Adding read flag to bad list
                    if read[0] in dataUnmapped.keys(): #Checking if read key already in unmapped dictionary
                        dataUnmapped[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                    else:
                        dataUnmapped[read[0]] = [read[1][i]] #If not: creating key and its corresponding read
            else: #If read flag in bad list
                if read[0] in dataUnmapped.keys(): #Checking if read key already in unmapped dictionary
                    dataUnmapped[read[0]].append(read[1][i]) #If yes: adding other read (mate) to key
                else:
                    dataUnmapped[read[0]] = [read[1][i]] #If not: creating key and its corresponding read            
    return dataMapped, dataUnmapped #Returning new dictionaries with only wanted reads (mapped and unmapped)

#Extracting properly mapped reads

def properlyMapped(data): #Parameter: (main) dictionary
    dataProperlyMapped = {} #Creating dictionary to store wanted reads
    listGood_1 = [] #Creating list to store flags corresponding to properly mapped reads ("list 1")
    listGood_2 = [] #Creating list to store flags corresponding to properly mapped mates ("list 2")
    for read in data.items(): #Going through (main) dictionary
        if len(read[1]) == 2: #Checking if key contains two reads (read + mate)
            bol = False #Setting default boolean value to false
            i = 0 #Setting i to 0
            while bol == False and i < len(listGood_1): #As long as boolean is False and i is inferior to length of list 1
                if read[1][0][0] == listGood_1[i] and read[1][1][0] == listGood_2[i]: #Checking if read flag in list 1 and mate flag in list 2
                    dataProperlyMapped[read[0]] = read[1] #If yes : adding read and mate to dictionary
                    bol = True #Changing boolean value to True
                i = i+1 #Increasing i by one
            if bol == False: #Checking if boolean still False after i goes over list length in previous loop
                if ((flagSearch(2,int(read[1][0][0])) and flagSearch(16,int(read[1][0][0]))) or (flagSearch(2,int(read[1][0][0])) and flagSearch(32,int(read[1][0][0])))) and ((flagSearch(2,int(read[1][1][0])) and flagSearch(32,int(read[1][1][0]))) or (flagSearch(2,int(read[1][1][0])) and flagSearch(16,int(read[1][1][0])))): #Checking if read and mate are properly mapped using flags
                    dataProperlyMapped[read[0]] = read[1] #If yes: adding read to dictionary
                    listGood_1.append(read[1][0][0]) #Adding read flag to list 1 
                    listGood_2.append(read[1][1][0]) #Adding mate flag to list 2
    return dataProperlyMapped #Returning new dictionary with only wanted reads

#Counting reads in a dictionary

def readCounter(data): #Parameter: any of our dictionaries
    counter = 0 #Setting counter variable to 0
    for read in data.values(): #Going through reads stored in each key of dictionary
        counter += len(read) #Adding number of reads to counter
    return counter #Returning counter

#Finding percentage

def percentage(dataLen, newDataLen): #Parameters: main dictionary length, new dictionary length
    percentage = round((newDataLen/dataLen)*100, 3) #Calculating percentage
    return percentage #Returning percentage
 
#Printing data summary

def summary(qValue, data, dataQuality, dataCigar, dataMapped, dataUnmapped, dataProperlyMapped): #Parameters: minimum quality value, main dictionary, new dictionaries with extracted reads
    dataLen = readCounter(data) #Assigning main dictionary length to a variable
    dataQualityLen = readCounter(dataQuality) #Assigning quality dictionary length to a variable
    dataCigarLen = readCounter(dataCigar) #Assigning matching dictionary length to a variable
    dataMappedLen = readCounter(dataMapped) #Assigning mapped dictionary length to a variable
    dataUnmappedLen = readCounter(dataUnmapped) #Assigning unmapped dictionary length to a variable
    dataProperlyMappedLen = readCounter(dataProperlyMapped) #Assigning properly mapped dictionary length to a variable
    print("\n") #Skipping a line
    print("Your file contains a total number of", dataLen, "reads.") #Printing total number of reads in main dictionary
    print("The number of reads above", qValue, "is", dataQualityLen, "reads, which corresponds to", percentage(dataLen, dataQualityLen), "percents of total reads.") #Printing number and % of reads in data dictionary
    print("The number of matching reads", dataCigarLen, "reads, which corresponds to", percentage(dataLen, dataCigarLen), "percents of total reads.") #Printing number and % of reads in matching dictionary
    print("The number of mapped reads is", dataMappedLen, "reads, which corresponds to", percentage(dataLen, dataMappedLen), "percents of total reads.") #Printing number and % of reads in mapped dictionary
    print("The number of unmapped reads is", dataUnmappedLen, "reads, which corresponds to", percentage(dataLen, dataUnmappedLen), "percents of total reads.") #Printing number and % of reads in unmapped dictionary
    print("The number of properly mapped reads is", dataProperlyMappedLen, "reads, which corresponds to", percentage(dataMappedLen, dataProperlyMappedLen), "percents of mapped reads and", percentage(dataLen, dataProperlyMappedLen), "percents of total reads.") #Printing number and % of reads in properly mapped dictionary
    print("\n") #Skipping a line

#Writing SAM file of dictionary of reads

def writeSAMfile(fichier, Type, header, data): #Parameters: original sam file, type of extraction, sam file header, dictionary of reads
    fileName = fichier.rstrip(".sam")+"_"+Type+".sam" #Creating name of new file
    with open(fileName, 'w') as f: #Opening new file (and closing after loop is completed)
        for item in header: #Going through each line stored in header list
            f.write(item) #Writing header line
        for read in data.items(): #Going through dictionary
            for i in range(len(read[1])): #Going through each read stored in key
                f.write(read[0]+"\t"+"\t".join(read[1][i])) #Writing key (read name) and each read information, separated by tabulation


############### DEFINING MAIN ###############


def main():
    fichier = check() #Storing file name in variable (after checking)
    response = question() #Asking user wanted minimum quality value
    start = time.time() #Starting timing of program
    header, dico = readAndStore(fichier) #Reading and storing header and file content in variables
    dicoQuality = extractQualityReads(dico, response) #Creating and storing quality dictionary
    dicoCigar = testCigar(dico) #Creating and storing matching dictionary
    dicoMapped, dicoUnmapped = mapped_unmapped(dico) #Creating and storing mapped and unmapped dictionaries
    dicoProperlyMapped = properlyMapped(dicoMapped) #Creating and storing properly mapped dictionary
    summary(response, dico, dicoQuality, dicoCigar, dicoMapped, dicoUnmapped, dicoProperlyMapped) #Printing data summary
    
    writeSAMfile(fichier, "quality_above_"+str(response), header, dicoQuality) #Creating file with quality reads
    writeSAMfile(fichier, "matching_reads", header, dicoCigar) #Creating file with matching reads
    writeSAMfile(fichier, "mapped_reads", header, dicoMapped) #Creating file with mapped reads
    writeSAMfile(fichier, "unmapped_reads", header, dicoUnmapped) #Creating file with unmapped reads
    writeSAMfile(fichier, "properly_mapped_reads", header, dicoProperlyMapped) #Creating file with properly mapped reads

    end = time.time() #Stopping timing of program
    print("Time taken:", end - start, "seconds.") #Print the time
    print("\n") #Skipping a line


############### LAUNCHING THE SCRIPT ###############


if __name__ == "__main__": #Checking if file runs as a script
    main() #If yes: launching program