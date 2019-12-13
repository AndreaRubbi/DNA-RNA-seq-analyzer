__author__ = 'Andrea Rubbi'

import string
import random
import turtle
import urllib.request

dnaletters=['a','t','g','c','A','C','G','T']
           
rna2amino = {'UUU':'F','UUC':'Y','UUA':'L','GGG':'G',
             'UUG':'L','CUU':'L','CUC':'L','CUA':'L',
             'CUG':'L','AUU':'I','AUC':'I','AUA':'I',
             'AUG':'M','GUU':'V','GUC':'V','GUA':'V',
             'GUG':'V','UCU':'S','UCC':'S','UCA':'S',
             'UCG':'S','CCU':'P','CCC':'P','CCA':'P',
             'CCG':'P','ACU':'T','ACC':'T','ACA':'T',
             'ACG':'T','GCU':'A','GCC':'A','GCA':'A',
             'GCG':'A','UAU':'Y','CAU':'H','CAC':'H',
             'CAA':'Q','CAG':'Q','AAU':'N','AAC':'N',
             'AAA':'K','AAG':'K','GAU':'D','GAC':'D',
             'GAA':'E','GAG':'E','UGU':'C','UGC':'C',
             'UGG':'W','CGU':'R','CGC':'R','CGA':'R',
             'CGG':'R','AGU':'S','AGC':'S','AGA':'R',
             'AGG':'R','GGU':'G','GGC':'G','GGA':'G',}
symbols={'r':('a','g'),'y':('c','t'),'k':('g','t'),
             'm':('a','c'),'s':('c','g'),'w':('a','t'),
             'b':('c','g','t'),'d':('a','g','t'),'h':('a','g','t'),
             'v':('a','c','g'),'n':('a','c','g','t')}


"""dictionary rna2amino contains all the codon:amino-acid key:value pairs where:
        codon: a three-letter string object representing a codon
        amino-acid: an upper case one-letter string object representing the coded amino-acid
"""

def gcContent(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    ratio=(dna.count('c')+dna.count('g'))/len(dna)
    return ratio
   
"""Returns the GC content ratio of a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a real number between 0 and 1
    Example: gcContent("atcgttcaag") = 0.4
"""

def countCodon(dna, codon):
    if type(dna)!=str or type(codon)!=str:
        return ('String type required')
    if len(dna)==0 or len(codon)==0:
        return('No DNA or codon was insert')
    for i in dna:
        if i not in dnaletters:
            return 'This DNA contains invalid characters'
    for c in codon:
        if c not in dnaletters:
            return 'This Codon contains invalid characters'
    dna=dna.lower()
    codon=codon.lower()
    l=counter=0
    while l<len(dna):
        if dna[l:l+3]==codon:
            counter+=1
        l+=3
    return counter

"""Returns the number of (non-overlapping) occurrences of a codon in a DNA sequence
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value: the integer number of instances of the target codon in dna
    Example: countCodon("aaaaaaaa", "aaa") = 2
"""

def countACG(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    counter=len(dna)-dna.count('t')
    return counter
                   
    """Returns the number of nucleotides that are not T in a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: the integer number of nucleotides in dna that are not T
    Example: countACG("atcgttcaag") = 7
    """

def printCodons(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return ('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    while len(dna)>=3:
        print (dna[:3],end=' ')
        dna=dna[3:]
      
"""Prints the sequence of non-overlapping codons in dna, 
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printCodons("ggtacactgta") would print: ggt aca ctg
    """
def kmp_matcher(T,P):
    prf=[0]
    L = []
    u=i=j=0
    m=len(P)
    n=len(T)
    for i in range(1,m):
        while u>0 and not P[u] == P[i]:   
            u = prf[u-1]
        if P[u] == P[i]:    
            u += 1
        prf.append(u)
    while i < n:
        while j > 0 and not P[j] == T[i]:    
            j = prf[j-1]
        if P[j] == T[i]:   
            if j == m-1:      
                L.append(i-j)
                j = prf[j-1]     
            else:
                j += 1
                i += 1
        else:
            i+=1            
    return L

def findCodon(dna, codon):
    if type(dna)!=str or type(codon)!=str:
        return ('String type required')
    if len(dna)==0:
        return ('No DNA was insert')
    if len(codon)==0:
        return('No codon was insert')
    for j in codon:
        if j not in dnaletters:
            return('The codon contains invalid characters')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    return kmp_matcher(dna.lower(),codon.lower())

    
    
    """Returns the index of the first occurrence of codon in dna
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value:
        (if codon found): the integer index of the first occurrence of codon in dna
        (if codon not found): None
    Example: findCodon("ggtacactacgta", "tac") = 2 
    """

def findATG(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    #dna=dna.lower()
    return kmp_matcher(dna.lower(),'atg')
    """Returns a list of all the positions of the codon ATG in dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a list of integer numbers
    Example: findATG("gatgtatgta") = [1,5] 
    """

def printReadingFrames(dna):
    if type(dna)!=str:
        return('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    l=0
    while len(dna[l:])>=3:
        a=dna[l:l+3]
        b=dna[l+1:l+4]
        c=dna[l+2:l+5]
        if len(dna[l+1:])<3 and len(dna[l+2:])<3:
            print(a)
        elif len(dna[l+2:])<3 and len(dna[l+1:])>=3:
            print(a, b,sep='\t')
        else:
            print (a,b,c,sep='\t')
        l+=3
        
    """Prints the sequences of non-overlapping codons in the dna
    with the three possible reading frames in separate columns
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printReadingFrames("aggcctggc") should print
        agg    ggc    gcc
        cct    ctg    tgg
        ggc
    """

def firstSSR(dna, seq):
    if type(dna)!=str:
        return('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    for c in seq:
        if c not in dnaletters:
            return 'This sequence contains invalid characters'
    #dna=dna.lower()
    matching=kmp_matcher(dna.lower(),seq.lower())
    n=0
    for c in matching:
        if n!=0:
            if c-cp>len(seq):
                break
        n+=1
        cp=c
    return n
        

"""Returns the length (number of repeats) of the first SSR in dna
    that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer number of seq repeats in the first SSR
        (if seq not found): 0
    Example: firstSSR("aggcctggcggcggc", "ggc") = 1
"""

def longestSSR(dna, seq):
    if type(dna)!=str:
        return('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    for c in seq:
        if c not in dnaletters:
            return 'This sequence contains invalid characters'

    matching=kmp_matcher(dna.lower(),seq.lower())
    n=0
    SSR=[]
    if len(matching)==0:
        return
    for c in matching:
        if n!=0:
            if c-cp>len(seq) or c==matching[-1]:
                SSR+=[n]
                n=1
        n+=1
        cp=c
    n=sorted(SSR,reverse=True)
    return n[0]
    
    
    """Returns the length of the longest SSR in dna that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer length of longest SSR in dna repeating seq
        (if seq not found): 0
    Example: longestSSR("aggcctggcggcggc", "ggc") = 3
    """

def longestSSRdin(dna):
    if type(dna)!=str:
        return('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    pos=('a','c','g','t')
    pairs=[]
    def thisdinucleotide(dna,din):
        i=count=0
        longest=0
        while i+len(din)<=len(dna):
            if dna[i:i+len(din)]==din:
                while dna[i:i+len(din)]==din:
                    count+=1
                    i+=len(din)
                if count>longest:
                    longest,count=count,0
            i+=1
            count=0
        return longest
    for p in pos:
        sett=[p+x for x in pos]
        for q in sett:
            count=thisdinucleotide(dna,q)
            pairs+=[(q,count)]
    def second ( t ) :
        return t [ 1 ]
    winner=sorted(pairs, key=second ,reverse=True)
    return winner[0]
    
    """Finds the longest SSR in dna for all the possible dinucleotides
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: (if len(dna)>1): a pair (din, len)
        din: a two-letter string object representing the dinucleotide with longest SSR in dna
        len: integer representing the length of the longest SSR of din
                  (if len(dna)<2): None
    Example: longestSSRdin("ctctctgcgccacacaca") = ("ca", 4)
    """

def complement(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    compl=''
    complementary={'a':'t','t':'a','g':'c','c':'g'}
    for c in dna:
        if c in complementary:
            compl+=complementary[c]
    return compl
        
    """Returns the complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the complement of the DNA sequence
    Example: complement("acgtac") = "tgcatg"
    """

def reverseComplement(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    return complement(dna)[::-1]
       
    """Returns the reverse  complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the reverse complement of the DNA sequence
    Example: reverseComplement("acgtac") = "gtacgt"
    """
    
def palindrome(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    return dna.lower()==reverseComplement(dna)
        
    """Returns true if dna is the same as its reverse complement
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: Bool
    Example: palindrome("atat") = True
    """

def dna2rna(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    return dna.replace('t','u')
           
    """Returns a copy of dna in which every "t" has been replaced by a "u"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: dna2rna("actgat") = "acugau"
    """

def transcribe(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    return dna2rna(reverseComplement(dna))

    """Returns the RNA equivalent of the reverse complement of dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: transcribe("acgtac") = "guacgu"
    """

def clean(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    dna=dna.lower()
    elegant=''
    for c in dna:
        if c in dnaletters:
            elegant+=c
        else:
            elegant+='n'
    return elegant

    """Returns a new DNA string in which every character in dna
    that is not an "a", "c", "g", or "t" is replaced with an "n"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing a clean DNA sequence
    Example: clean("goat") = "gnat"
    """

def fix(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    dna=dna.lower()
    fixed=''
    for n in dna:
        if n in symbols:            
            fixed+=random.choice(symbols[n])
        else:
            fixed+=n
    return fixed
        
    """Returns a DNA string in which each ambiguous symbol is replaced
    with one of the possible bases it represents, each with equal probability
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a string object representing a fixed DNA sequence
    Example: fix("arta") returns either "aata" or "agta"
        (with probability 1/2 for each one of these two possible returns)
    """

def fixAll(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    fixed=['']
    for n in dna.lower():
        if n in symbols:
            fixed=[x+y for y in symbols[n] for x in fixed]
        else:
            fixed=[x+n for x in fixed]     
    return fixed                            
    """Returns the list containing all possible DNA strings in which each
    ambiguous symbol in dna is replaced with all of the possible bases it
    represents; all combinations of replacements should be considered and
    returned within the list of fixed DNA strings
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a list of string objects representing fixed DNA sequences
    Examples:   fixAll("arma") = ["aaaa", "aaca", "agaa", "agca"]
                fixAll("agata") = ["agata"]
    """

def readFASTA(filename):
    fasta=open(filename,'r')
    sequence=''
    read=fasta.readlines()
    for line in range(1,len(read)):
        sequence+=read[line][:-1]
    fasta.close()
    return sequence
          
    """Reads a FASTA file and removes the header and all the newline characters
    Returns the DNA sequence contained in the file as a string
    Parameter:
        filename: the name of file containing a DNA sequence in the FASTA format
    Return value: a string object representing the DNA sequence in the file
    """

def getFASTA(id):
    prefix = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
    suffix = '&rettype=fasta&retmode=text'
    url = prefix + id + suffix
    read_file = urllib.request.urlopen(url)
    first_line = read_file.readline()
    dna_sequence = ''
    for line in read_file:
        line = line.decode('utf-8')
        dna_sequence += line[:-1]
    read_file.close()
    return dna_sequence
    """Fetches the DNA sequence with the given id from the NCBI database
    and returns it as a string (header and newline characters must be removed)
    Parameter:
        id: a string object representing the identifier (NCBI accession number) of a DNA sequence
    Return value: a string object containing the dna sequence
    """

def findseq(dna, seq):
    if len(dna)==0 or len(seq)==0:
        return('No DNA or seq was insert')
    for i in dna:
        if i not in dnaletters:
            return 'This DNA contains invalid characters' 
    for c in seq:
        if c not in dnaletters:
            return 'This sequence contains invalid characters'
    #dna=dna.lower() 
    return kmp_matcher(dna.lower(),seq.lower())

    """Returns the list of indexes (starting positions) of all the
    occurrences of seq in dna
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a sequence of DNA
    Return value: a list of integer numbers
    """

def mark(dna):
    if type(dna)!=str:
        return ('String type required')
    if len(dna)==0:
        return('No DNA was insert')
    for i in dna:
        if i not in dnaletters:
            return 'The DNA sequence contains invalid characters'
    dna=dna.lower()
    marked=''
    x=0
    y=0
    while x+3<=len(dna):
        if y==0 and dna[x:x+3]=='atg':
            marked+='>>>'
            x+=3
            y=1
        elif y==1 and (dna[x:x+3]=='taa' or dna[x:x+3]=='tag' or dna[x:x+3]=='tga'):
            marked+='<<<'
            x+=3
            y=0
        else:
            marked+=dna[x:x+3]
            x+=3         
    return marked
    
    """Returns a new DNA string in which every start codon (atg) in dna
    is replaced with ">>>" and every stop codon (taa, tag, or tga) is
    replaced with "<<<"
    The function does not consider overlapping codons (but just the
    reading frame starting from offset 0)
    If two subsequent start codons are found (without a stop codon in
    between), only the first one should be replaced with ">>>", since
    codon atg also codes for the Methionine amino-acid
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object represnting a marked sequence of DNA
    Example: mark("ttgatggagatgcattagaag") = "ttg>>>gagatgcat<<<aag"
    """

def proteins(marked_dna):
    if type(marked_dna)!=str:
        return ('String type required')
    if len(marked_dna)==0:
        return('No DNA was insert')
    for i in marked_dna:
        if i not in dnaletters and i!='<' and i!='>' :
            return 'The DNA sequence contains invalid characters'
    rna=marked_dna.replace('t','u')
    rna=rna.upper()
    x=0
    protein=[]
    while x+3<=len(rna):
        if rna[x:x+3]=='>>>':
            x+=3
            i = ''
            while rna[x:x+3]!= '<<<' and x+3<=len(rna):
                if rna[x:x+3]=='>>>':
                    i=''
                else:
                    i+=rna2amino[rna[x:x+3]]
                x+=3
            protein.append(i)     
        else:
            x+=3
    return protein

"""Returns the list of proteins traduced from marked_dna
    Proteins are represented as strings of amino-acids
    Proteins are obtained from the RNA sequences obtained from
    the dna enclosed between the markers ">>>" and "<<<"
    Parameter:
        marked_dna: a string object representing a marked sequence of DNA
    Return Value: a list of string objects representing proteins
    Example: proteins("ttg>>>gagcat<<<aagcag>>>aca>>>caccaacag<<<aga") = ["EH","HQQ"]
    """



""" The following lines define the base settings and functions for plotting through the turtle module
You do not need to add any code to the functions plot and bar!
"""

width = 1200		# width of the window
cols = width // 6	# number of columns of text
height = 600		# height of the window
rows = height // 100	# number of rows of text



def plot(tortoise, index, value, window):
	"""Plots GC fraction value for window ending at position index."""
	
	if (index == window) or (index - window + 1) // cols != (index - window) // cols:
		tortoise.up()	
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)
		tortoise.down()
	else:
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)

		
def bar(tortoise, index, rf):
	"""Draws a colored bar over codon starting at position index in
	   reading frame rf. Puts the turtle's pen up and down to
	   handle line breaks properly."""
	   
	tortoise.up()
	tortoise.goto(index % cols, index // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 1) % cols, (index + 1) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 2) % cols, (index + 2) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)


"""
Complete the following two functions to accomplish the required tasks
"""

def orf(dna, rf, tortoise):
    if type(dna) != str:
        return 'The parameter you have inserted is of an invalid type.'
    if dna == '':
        return 'You have inserted an empty DNA sequence.'
    dna = dna.lower()
    tortoise.color('red')
    i=rf
    while i+3 <= len(dna):
        m=dna[i:i+3]
        if m!='atg':
            bar(tortoise,i,rf)
            tortoise.color('red')
            i+=3
        else:
            while m not in ('taa','tga','tag') and i+3 <= len(dna):
                tortoise.color('blue')
                bar(tortoise,i,rf)
                i+=3

    """Finds and draws all ORFs in the reading frame rf
    Blue bars begin and end on start and stop codons
    Parameters:
        dna: a string object representing a sequence of DNA
        rf: reading frame offset (it's value can be either 0, 1, or 2)
        tortoise: the drawing turtle
    Return value: None
    """
	   
    # YOUR CODE GOES HERE
	
    # to place a bar in the current color over the codon starting at
    # position index in reading frame rf, call
    # bar(tortoise, index, rf)

	

def gcFreq(dna, window, tortoise):
    """Computes and plots the GC frequency in dna over a sliding window
    Parameters:
        dna: a string object representing a sequence of DNA
        window: integer size of the sliding window
        tortoise: the drawing turtle
    Return value: None
    """

    # draws red lines at 0.5 above the sequence:
	
    tortoise.pencolor('red')
    for index in range(len(dna) // cols + 1):
    	tortoise.up()
    	tortoise.goto(0, index + 0.825)
    	tortoise.down()
    	if index < len(dna) // cols:
    		tortoise.goto(cols - 1, index + 0.825)
    	else:
    		tortoise.goto((len(dna) - window) % cols, index + 0.825)
    tortoise.up()
    tortoise.pencolor('blue')
    k=0
    while k+window<len(dna):
            plot(tortoise, window+k,  gcContent(dna[k:k+window]), window)
            k+=1
	
    # YOUR CODE GOES HERE

    # get initial window count
	
    # get subsequent window counts and plot them
    # to plot a fraction for the window ending at position index,
    # call plot(tortoise, index, fraction, window)



"""
The viewer functuion calls the functions you have coded to build the final plot!
"""

def viewer(dna):
    """Displays GC content and ORFs in 3 forward reading frames."""

    dna = dna.upper()   # makes everything upper case

    tortoise = turtle.Turtle()
    screen = tortoise.getscreen()
    screen.setup(width, height)     # makes a long, thin window
    screen.setworldcoordinates(0, 0, cols, rows) # scales coord system so 1 char fits at each point
    screen.tracer(100)
    tortoise.hideturtle()
    tortoise.speed(0)
    tortoise.up()

    # Prints the DNA string in the window:
	
    for index in range(len(dna)):
    	tortoise.goto(index % cols, index // cols)
    	tortoise.write(dna[index], font = ('Courier', 9, 'normal'))
		
    # Finds ORFs in forward reading frames 0, 1, 2:
	
    tortoise.width(5)
    for rf in range(3):
    	orf(dna, rf, tortoise)
		
    # Plots GC frequency:
	
    tortoise.width(1)
    gcFreq(dna, 5, tortoise)

    screen.update()
