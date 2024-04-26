#!/usr/bin/python3
"""
Created on Tue Apr 16 15:32:15 2024

@author: ariannaalamshahi
"""
'''
sparrows_Alamshahi.py - modified and expanded version of Joel's original script (sparrows.py)

Description: This script will parse VCF input for the hybrid italian sparrow, spanish sparrow, house sparrow, and tree sparrow. 
Using user-provided cutoffs, it will designate alleles per locus as spanish specific, house specific, or tree specific. 
Then designate as ancestral or derived. Includes "skip-rules" which specifies instances that will be 
overlooked/skipped for the sake of our research

See README for in-depth description of skip rules

User-defined functions: designate_alleles

Non-standard modules: none
    
Required inputs: VCF file and tsv file of bird IDs and associated population
Optional inputs: tmin & tmax (defaults if user doesn't profide), tsv with alternate chromosome names

Output: TSV that shows locus, allele frequency (House, Spanish, Tree), 
and designations (House specific/Spanish specific, Ancestral/Derived) for each sample 

Usage: Python3 sparrows_alamshahi.py -v input.vcf -b popmap.tsv -c 0.75,0.8,0.85,0.90,0.95,0.99,1.00 -n chromNames.tsv -o sparrowVCF_filtered

Date: April 25, 2024
Name: Arianna Alamshahi, original author: Joel
'''
#%% Import hub

# Use to run code from terminal
import sys

# Regex used in file format verification
import re

#%% Error handling hub - for custom exceptions! 

class WrongNumberArguments(Exception):
    pass

class WrongFormat(Exception):
    pass

class ArgumentNotSpecified(Exception):
    pass

#%%% Allow program to be ran from the terminal
            
# will only run if minimum 9 args provided, max 15 args, and the VCF and population tsv are both files
if 9 <= len(sys.argv) <= 15: #and Path(sys.argv[1:2]).is_file():
    
    sparrows_python = sys.argv[0] # sparrows_Alamshahi.py
    
    # Initialize a template dictionary to store flags and user provided arguments
    sys_dict = {'-h':'-v: the path to the VCF file\n-b: the path to the tsv file that assigns individual birds to groups/populations\n-c: a comma-delimited string of cut-offs\n-tmin: tree sparrows\' minimum allele count when analyzing triallelic loci, default 2\n-tmax: tree sparrows\' max minor allele count when analyzing biallelic loci, default 1\n-n: (optional) path to tsv file which has chromosome name used in vcf and corresponding chromosome name using chromosme numbers (e.g. chr1)\n-o: the prefix of the output file',
                '-v':None, '-b':None, '-c':None, '-tmin':2,'-tmax':1, '-n':None, '-o':None}
    
    # Add every sys.argv from index 1 onwards to this list
    sys_list = sys.argv[1:] 
    
    # Add arguments provided by user to sys_dict
    flag = None # set flag to None
    
    for argument in sys_list: # for each argument in sys_list
    
        if argument.startswith('-'): # if it starts with a '-'
            flag = argument # then we know it is a flag
            
        else:
            if flag != None: # if flag isn't empty
                sys_dict[flag] = argument # then we know the argument corresponds to the flag, so save it in the dict
                flag = None # reset flag to None
    
    # If any of the dictionary values are still none (except for -n flag), exit program
    # For example, if user puts dash before path, it won't be saved, leaving an empty entry
    try:
        for key, value in sys_dict.items(): # loop through the dictionary key (flag) & value (argument)
        
            if value == None: # if a flag doesn't have a specified argument
                if key != '-n': # and if that flag is not the optional -n flag
                    raise ArgumentNotSpecified # then raise an error
                    
    except ArgumentNotSpecified:
        print("\nOh no, it looks like you didn't specify the " + key + " flag!\n\nHere is the manual in case you need a refresh:\n" + sys_dict['-h'] + '!')
        sys.exit()
    
    # Assign path to VCF file as input_VCF
    input_VCF = sys_dict['-v']
    
    #Error handling! Don't allow a non-VCF file
    try:        
        vcf = open(input_VCF,'r') #open up the file that was entered for sys.argv[1]       
        headervcf = vcf.readline() # read the first line
             
        if headervcf.startswith("##fileformat=VCFv4.2") == False: # if the first line in [1] doesn't start with "##fileformat=VCFv4.2" then raise WrongFormat error
            raise WrongFormat    
            
    except FileNotFoundError: # this is a built in python exception
        print("Error! Cannot find the input VCF in the specified directory.")
        sys.exit()
        
    except WrongFormat: 
        print("Error! The input VCF doesn't look like a file that I know how to parse!")
        sys.exit()
        
    vcf.close() #close the file

    # Assign path to popmap file as input_popmap
    input_popmap = sys_dict['-b'] # popmap.
    
    
    #Error handling! Don't allow an incorrectly formatted birdID-population file in position [2]
    try:        
        popmap_file = open(input_popmap,'r')
        popmap = popmap_file.read() # read the contents of the file in sys.argv[2] into popmap  
        
        if re.match(r'\S*\t\S*\n*',popmap): # should follow format of ID, tab, population, newline
            pass
        else:
            raise WrongFormat() 
            
    except FileNotFoundError: # this is a built in python exception
        print("Error! Cannot find the pop input file in the current directory.")
        sys.exit()
        
    except WrongFormat: 
        print("Error! Make sure your birdID-population file is tab separated with each birdID-population pair on a new line")
        sys.exit()
        
    popmap_file.close() #close the file
    
    # Assign more arguments to variables
    cutoffs = sys_dict['-c']
    
    tmin = float(sys_dict['-tmin'])
    print(tmin)
    
    tmax = float(sys_dict['-tmax'])
    
    # Optional to provide tsv of alt chrom names!
    chromNames = None # still need to intialize chromNames variable (to check if empty or not down the line)
    
    if sys_dict['-n'] != None: # if there is an argument specified for chromNames
        chromNames = sys_dict['-n'] # then store that argument in chromNames variable
    
    # Assign output prefix argument to variable
    outputPrefix = sys_dict['-o']
    
# Error handling! If more than 15 or less than 9 arguments are provided
try:
    if len(sys.argv) > 15 or len(sys.argv) < 9:
        raise WrongNumberArguments
        
except WrongNumberArguments:
    print("Error! It looks like you have provided the wrong number of arguments. Goodbye!")
    sys.exit()

#%% Function Hub

''' The designate_alleles function will output 2 dictionaries. When we call on this function in the script, 
we will name the output directories as allele_designations and afs

allele_frequencies: uses alleles (0, 1, and 2) as keys and minor allele frequency ([housemaf], [spanishmaf], [treemaf]) as values (tuples)
allele_designations: uses alleles (0, 1, and 2) as keys and designation (ancestral, derived, neither) as values (tuples)


Per Joel: "It ignores missing calls, so if many calls are missing for e.g. Spanish sparrows, it is very unlikely 
that an allele can be designated 'Spanish', since the division does not take into account the missing ones, 
but assumes zero missingness. I think it's best this way"
'''

def designate_alleles(cut,alleles,private_allele):
    
    # Create an empty dictionary store allele frequencies
    allele_frequencies	= {}
    
    # Create an empty dictionary to store allele designation
    allele_designations	= {}
    

    # Loop through the alleles dictionary
    for allele in alleles: # allele = 0, 1, or 2 (these are the keys of this dictionary)

		# Determine allele frequencies
        house_af = alleles[allele]['House']/(2*(Nhouse-hm))
        # alleles[allele]['House'] -> pulls counts for each House allele (0, 1, and 2)
        # (Nhouse-hm) -> will subtract the number of missing House genotypes from the number of House individuals
        # Multiply by 2 because diploid
                
        spanish_af = alleles[allele]['Spanish']/(2*(Nspanish-sm))
        
        tree_af = alleles[allele]['Tree']/(2*(Ntree-tm))
        
        
        # Append the allele frequencies to the dictionary
        allele_frequencies[allele]=(str(house_af),str(spanish_af),str(tree_af)) # key is alelle (0, 1, 2) and values are tuples of strings with the frequencies

        
		# Is this allele Spanish or House? Does not depend on the number of alleles
        if spanish_af > cut and house_af < (1-cut): # if spanish allele freq is higher than the cutoff and house allele freq is lower than 1-cutoff
            spanish = True # we designate it as spanish
            house = False # we DO NOT designate it as house
            
        
        elif house_af > cut and spanish_af	< (1-cut): # if house allele freq is higher than the cutoff and spanish allele freq is lower than 1-cutoff
            spanish = False # we DO NOT designate it as spanish
            house = True # we designate it as house
            
        
        else:  # If neither of the above sets of conditions are met, then we cannot designate as spanish nor house
            spanish = False
            house = False
            
		
        # If biallelic, one REF and one ALT
        if nalleles == 2:
			
            # If the current allele has more counts than the maximum allowed minor allele count, it must be the major allele,
			# which implies it should be classified Ancestral. The reverse is true for the Derived allele.
            if	alleles[allele]['Tree'] > tmax: # implies ancestral
                tree = True
                
            else: # implies derived
                tree = False
        
        
        # If triallelic, one REF and two ALT
        elif nalleles == 3:
            
			# We can assume, because of the skipping rules, that the Tree sparrow has two alleles,
			# one of which is private. The not private one is thus Ancestral, and the other Derived.
			# The private allele will be called Tree, denoted with a 'T'. It's S/H value will be '-'.
            if allele == private_allele:
                tree = 'private'
                
            elif alleles[allele]['Tree'] == 0: # implies derived
                tree = False
                
            else: # implies ancestral
                tree = True
        
        # Error handling, can improve
        else:
            print('Something\'s very wrong...',file=sys.stderr) ; sys.exit()


		# Bring previous steps together with by appending to allele_designations dictionary
        if tree =='private':
            
            allele_designations[allele] = ('neither','tree') # then dict should indicate neither house nor spanish
            
            
        elif tree: # if tree == True, meaning the allele is found in tree sparrow
            
            if spanish and not house: # if spanish == True and house == False
            
                allele_designations[allele] = ('spanish','ancestral') # then the allele is designated spanish ancestral
                
            elif house and not spanish: # if house == True and spanish == False
            
                allele_designations[allele] = ('house','ancestral') # then the allele is designated house ancestral
                
            else: # if house == False and spanish == false
                allele_designations[allele] = ('neither','ancestral') # then the allele is ancestral without a parent designation
            
                
        else: # if tree == False, meaning the allele is not found in the tree sparrow
            
            if spanish and not house: # if spanish == True and house == False
            
                allele_designations[allele] = ('spanish','derived') # then allele is designated spanish derived
                
            elif house and not spanish:	 # if house == True and spanish == False
            
                allele_designations[allele] = ('house','derived') # then allele is designated house derived
                
            else: # if house == False and spanish == False
            
                allele_designations[allele] = ('neither','derived') # then the allele is derived without a parent designation
    
                
    # The function will output the dictionaries allele_designations and allele_frequencies     
    return allele_designations,allele_frequencies # When we call this function in the code, we will put allele_frequencies output directly in afs dictionary (shorter name, easier to call)

#%% Parse input - cutoffs

# Separate cutoffs by comma and store each cutoff in a list
strcuts = cutoffs.split(',')

# Initiate an empty dictionary
cutoffs_dict = {}

# Now a dict where values are lists to be filled
for cutoff in strcuts:
    
    cutoff_float = float(cutoff) # store the value as a float
    
    cutoffs_dict[cutoff_float] = [] # store cutoffs as keys and empty lists as values

#%% Parse input - birdID & associated population

# Create sets for existence checks later on
threelilbirbs = {'House','Spanish','Tree'}
italian = {'Sardinia','Malta','Corsica','Crete'}

# Create empty dictionaries for existence checks later on (when using the list instead of dict)
TLBindices = {} # TLB = threelilbirbs
Iindices = {} # I = Italian

# Create an empty dictionary which will store birdID as key and group/population as value
groups = {}

# Set total number of respective individuals, to simplify later calcs
Nhouse = 0 # number of house individuals etc.
Nspanish = 0
Ntree = 0


with open(input_popmap, 'r') as birbPops: # open the file containing birdID and associated population
    
    # loop through birbPops line by line
    for line in birbPops:
        
        birb,group = line.split() # split each line into birb (birdID) and group (population), based on whitespace
        
        groups[birb] = group # append to the dictionary we initialized, birb (birdID) is key and group (population) is value
        
        if group =='House': # if the group is "House"
            Nhouse += 1 # then increase the "House" count by 1
            
        elif group =='Spanish': # if the group is "Spanish"
            Nspanish += 1 # then increase the "Spanish" count by 1
            
        elif group =='Tree': # if the group is "Tree"
            Ntree += 1 # then increase the "Tree" count by 1

# Total number of birbs, we will use this to calculate missingness
Nbirbs = Nhouse + Nspanish + Ntree

#%% Parse input - pull the #CHROM header line
with open (input_VCF, 'r') as VCF:
    
    for header in VCF:
        
        if header.startswith("#CHROM"):
            
            break
        
# Create a list of just birdIDs
birds = header.split()[9:] # splits the line by whitespace and stores chunks from 9th index onwards (after FORMAT field, so just birdIDs)

# Loop through birds list
for i,birb in enumerate(birds): # enurmerate birds, returns an index (i) and associated birdID (birb)

    if groups[birb] in threelilbirbs: # use birb (birdID) as a key for dictionary and if the associated value is in threelilbirbs (House, Tree, Spanish)
        
        TLBindices[i] = birb # then append to the TLBindices dictionary where key is index and birb (birdID) is value
        
    elif groups[birb] in italian: # use birb (birdID) as a key for dictionary and if the associated value is in italian
        
        Iindices[i]=birb # then append to Iindices dictionary where key is index (i) and birb (birdID) is value

#%%
if chromNames != None: # if chromNames is not empty
    chromNames_dict = {} # initialize empty dictionary
    
    with open(chromNames, 'r') as chromNames: # open tsv
    
        for line in chromNames: # look through it line by line
        
            longChrom,numberChrom = (line.strip()).split('\t') # remove whitespace split the line by tab and save first entry to original, second entry to number
            
            chromNames_dict[longChrom] = numberChrom # update the dictionary with these entries so original is key and number is value

#%% Parse input - the VCF

with open(input_VCF, "r") as VCF:
    # Iterate through each line in the file
    for line in VCF:
        # Check if the line starts with #
        if not line.startswith('#'):
            
            # Separate the line by fields and store in a list
            fields = line.split() # splits "line" based on whitespaces and save each section into a list of strings called "fields"

            # Create a string which stores locus, chromosome_position_refAllele_altAllele
            locus = '_'.join([fields[0], fields[1], fields[3], fields[4]]) # pull chromosome, position, refAllele, altAllele and join by underscore
            
            if chromNames != None: # if chromNames is not empty
                updated_locus = locus
                for longChrom, numberChrom in chromNames_dict.items():
                    updated_locus = updated_locus.replace(longChrom, numberChrom)


            # Determine the number of alleles at this locus
            nalleles = fields[4].count(',') + 2 # counts the number of commas present in the ALT allele section and then adds 2 (will defintiley have 1 value for REF and 1 for ALT)
            # e.g. if there is one comma, means there are 2 ALTs and 1 Ref -> this calc will output 3

            # If it is a multiallelic site:
            if nalleles > 3:
                continue # Then skip-rule 4
            
            # Initiate an empty list to store genotypes at that locus
            locusdata = []

            # Loop over birbs' FORMAT fields
            for data in fields[9:]:
                
                # Store the genotype at this locus for each birb as a list of genotypes (e.g [0,1]) within a list
                locusdata.append(data.split(':')[0].split('/')) # separate the string based on colons, then pull the first chunk (genotype) and split it at the slash

            # Skip-rules 1, 2, and 3.
            threshold = 0.5 # establish a missingness threshold

            # Number of birbs of each type missing genotype at this locus, integer value
            hm = 0  # House, missing genotype at locus
            sm = 0  # Spanish, missing genotype at locus
            tm = 0  # Tree, missing genotype at locus

        
            for x, y in zip(locusdata, birds): # locusdata is x and birds [sampleID] is y
            
                if '.' in x: # if there is a period in the locus data, indicates misisng allele
                
                    # Use y (birdID) as key to access value in groups dict
                    if groups[y] == 'House': # if the y is associated with value "House" in groups dict
                        hm += 1 # then increase number of missing House genotypes by 1
                        
                    elif groups[y] == 'Spanish':
                        sm += 1
                        
                    elif groups[y] == 'Tree':
                        tm += 1

            if hm/Nhouse > threshold: # if (number of missing House missing genotypes at this locus) / (number of House) is greater than threshold
                continue # then start new iteration of loop
                
            elif sm/Nspanish > threshold: 
                continue
            
            elif tm/Ntree > threshold: 
                continue
                
            # Establish a nested allele count dictionary
            # alleles (0, 1, and 2) as keys and nested within the values are keys for House/Spanish/Tree with associate counts as values
            if nalleles==2: 
                alleles={'0':{'House':0,'Spanish':0,'Tree':0},'1':{'House':0,'Spanish':0,'Tree':0}}
            elif nalleles==3: 
                alleles={'0':{'House':0,'Spanish':0,'Tree':0},'1':{'House':0,'Spanish':0,'Tree':0},'2':{'House':0,'Spanish':0,'Tree':0}}
            
            # Loop over birbs in order
            for i,genotype in enumerate(locusdata): # i is the index for the genotype, genotype = genotype in locusdata
            # Will create a list called genotype with two entries
            
                # Loop through dictionary of indexes and TLB birdIDs
                if i in TLBindices: # if ID belongs to Spanish/House/Tree
                    
                    if '.' in genotype: # no genotype for this birb at this locus
                        continue # ignore missing birbs, start next iteration of loop
                    
                    
                    for allele in genotype: # for each allele in the genotype list
                        
                        # Adjust the allele count
                        alleles[allele][groups[TLBindices[i]]] += 1
                        # access alleles dictionary with TLBindices[i] as key to pull birdID
                        # use that birdID as key in groups dictionary to pull group (population)
                        # for each 0,1,2 allele present -> increase allele count for that allele for that population
                        

    		# The remaining skip-rules
            if nalleles == 2: # if biallelic
                
                private_allele = False # then no private allele for tree
                
    			# Skip-rule 5
                # If the lowest alelle count out of allele 0 and 1 for Tree is less than tmax
                if min(alleles['0']['Tree'],alleles['1']['Tree']) > tmax:
                    continue # then start next iteration of loop
                
            if nalleles == 3: # if triallelic, one REF two ALTs
    			
                # First find out if there is a private allele (just one, because two implies Spanish == House)
                # Keep track of the number of private alleles (0, 1, or 2)
                pvtcount = 0 
                count = 0
                
                for allele in alleles: # loop through the alleles dict, allele (0, 1, 2) as key
                    
                    # if the count of the allele for House and Spanish are 0, and the count for Tree is not 0
                    if alleles[allele]['House'] == alleles[allele]['Spanish'] == 0 and alleles[allele]['Tree']!=0 :
                        
                        private_allele = allele # then that allele is a private allele
                        
                        pvtcount += 1 # and we should increase the private allele count by 1
                        
                    elif alleles[allele]['Tree'] != 0: # if the allele count for tree is not 0, but the other conditions aren't met
                        
                        count += 1 # then increase the non-private allele count by 1

    			# Skip-rule 6:
                if not pvtcount == 1 or not count == 2: # if the private allele count is one, or the non-private count is not 2
                    continue # then skip-rule 6
                    
    			#Skip-rule 7:
                if alleles[private_allele]['Tree'] < tmin: # use the private allele as key for alleles, if the count of Tree for that allele is less than tmin
                    continue # then skip-rule 7


    		# Create long lists per user supplied cutoff
            for cut in cutoffs_dict: # "cut" pulls cutoffs as keys
            
                
                # Run the function designate_alleles(cuttof value,alleles dictionary,private_allele(true or false))
                allele_designations, afs = designate_alleles(cut,alleles,private_allele) # save outputs (allele_designatations and allele_frequencies) into allele designations and afs respectively

                # Create a list to eventually store all items we want to print in a row for output
                if chromNames != None: # if chromNames is not empty, then we will include the updated locus
                    row = [locus,updated_locus,str(nalleles)] # first entry is locus, second is number of alleles (string)
                    
                elif chromNames == None: # if chromNames is empty, then there is no alternate chrom name and no updated locus
                    row = [locus,str(nalleles)]
    			
                # Pull values for the first three allele frequency columns
                if nalleles==2: # if biallelic
                    housemaf = afs['0'][0]+'/'+afs['1'][0] # extracts allele frequencies for house sparrow as a string separated by "/"
                    
                    spanishmaf = afs['0'][1]+'/'+afs['1'][1]
                    
                    treemaf = afs['0'][2]+'/'+afs['1'][2]
                    
                if nalleles==3: # if triallelic
                    housemaf = afs['0'][0]+'/'+afs['1'][0]+'/'+afs['2'][0] # extracts allele frequencies for house sparrow as a string separated by "/"
                    
                    spanishmaf = afs['0'][1]+'/'+afs['1'][1]+'/'+afs['2'][1]
                    
                    treemaf = afs['0'][2]+'/'+afs['1'][2]+'/'+afs['2'][2]
                    
                
                # Append the house minor allele frequency string to the row list
                row.append(housemaf)
                row.append(spanishmaf)
                row.append(treemaf)
                
                # The birbs themselves
                skip = True # set skip equal to True to begin with
                
                for i,genotype in enumerate(locusdata): # loop through locusdata with i as index and genotype as the current genotype in the list
                    
                    # Check the index of the genotype against the index of italian sparrows
                    if i in Iindices: # if this is an italian sparrow
                        
                        if '.' in genotype: # if the genotype is missing
                            
                            row.append('-') # then append a "-" to the row list
                            
                            continue # This birb has nodata, move directly on to the next birb
                        
                        
                        # Initiate an empty list
                        cats = [] # container for the 'AB/CD' thing
                        
                        for allele in genotype: # for alleles in genotype (of italian sparrows)
                            
                            if 'spanish' == allele_designations[allele][0]: # use allele as key to access allele_designations dict, [0] index will pull the house/spanish/neither designation

                                cats.append('S') # append the letter S to the cats list
                                
                                skip = False # set skip equal to false, skip-ruel
                            
                            elif 'house' == allele_designations[allele][0]:
                                
                                cats.append('H') # append the letter H to the cats list
                                
                                skip = False # set skip equal to false
                            
                            else: # when allele was designated 'neither'
                            
                                cats.append('-') # append a dash to the cats list
                            
                            if 'ancestral' == allele_designations[allele][1]: # use allele as key to access allele_designations dict, [1] index will pull the ancestral/derived/tree designation
                                cats.append('A') # append the letter A to the cats list
                            
                            elif 'derived' == allele_designations[allele][1]:
                                cats.append('D') # append the letter D to the cats list
                            
                            else:  # when allele was designated 'tree'
                                cats.append('T') # append the letter T to the cats list
                                
                            
                        cats.insert(2,'/') # append "/" at index [2] of the list to separate the designation for each allele
                        
                        row.append(''.join(cats)) # join all entries currently in cats as one string, then append to row
                            
                if skip: # if skip is still equal to true
                    continue # Skip-rule 9
                
                # use the cuttoff as a key, and append the joined row (all entries of list on one line) with a tab separator
                cutoffs_dict[cut].append('\t'.join(row))
                
#%% Write output

# Set the header for the output file
if chromNames != None: # if chromNames is not empty
    header = ['#Locus','locus_AltFormat','num_alleles','House_af','Spanish_af','Tree_af']+[bird for bird in birds if groups[bird] in italian]
elif chromNames == None: # if chromNames is empty
    header = ['#Locus','num_alleles','House_af','Spanish_af','Tree_af']+[bird for bird in birds if groups[bird] in italian]


# Loop through the cuttofs dictionary and the list of cutoffs
for cut,cutstr in zip(cutoffs_dict,strcuts):
    
    with open(outputPrefix+'.designated_cutoff_'+cutstr+'.tsv','w') as f:
            
        f.write('\t'.join(header)+'\n')
            
        for row in (cutoffs_dict)[cut]:
            
            f.write(row+'\n')

