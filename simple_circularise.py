
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="[Input]: list of  sequences to be circularised")
parser.add_argument("output_fasta_file", help ="[Output]: circularised sequences")  
parser.add_argument("-p","--probability",type=float, help="custom probability (default = 0.005)")
parser.add_argument("-max","--maximum_size",type=int,help="maximum size of output sequence")
parser.add_argument("-min","--minimum_size",type=int,help="minimum size of output sequence")
parser.add_argument("-r","--repeat",type=int,help="change behaviour to optimise repeat size (default: optimise genome size)")
args = parser.parse_args()
if args.probability:
    print ("custom probability turned on") 
if args.maximum_size:
    print("maxium_size turned on")
if args.minimum_size:
    print("minimum_size turned on")    
if args.repeat:
    print("behaviour set to optimise size of repeat")
print("Circularising...")
import math

##calculates the probability of a string of length L occuring in a DNA seq.
def seq_prob (L):
    P = 0.25 ** L
    return P 

##determine the number of events in a sequence of size S with a sub string size of L
def events (S,L):
    P = S - (L - 1)
    return P 

#calculate probability of exactly b successes
# ** raises number to a power. 
##where: a = number of events, b = number of successes, c = probability of success     

def binomial_prob(a,b,c):
    d = 1 - c
    P = math.factorial(a)/ (math.factorial(a-b)* math.factorial(b)) * c**b * d**(a-b) 
    return P

# Use the poisson distribution to approximate 
def poisson(lmbda,b):
    e = 2.71828
    P = (lmbda ** b * e ** -lmbda) / math.factorial(b) 
    return P 

def lambd_a (eventn,sprob):
    lmbda = (eventn * sprob)
    return lmbda
    
#calculate probability of at least two occurances of sequence 
#  where a = number of events and c = probability of success
def cooccur_prob(a,c,lmbda,SeqL):
    d = 1 - c
    noccur = (d ** a)
    occur = 1 - noccur
    if SeqL > 9999:
        once = poisson(lmbda,1)
    else:
        once = binomial_prob(a,1,c)
    P = occur - once
    return P

##returns the probability that the linear contig is circular. That is the probability of a n base region being repeated
# twice or more
def pcircular(SeqL,SLength): 
    sprob = seq_prob(SLength)
    eventn = events(SeqL,SLength)
    lmbda = lambd_a(eventn,sprob)
    #print(sprob)
    #print(eventn)
    co_prob = cooccur_prob(eventn,sprob,lmbda,SeqL)
    #print(co_prob)
    return co_prob


##for a sequence of size S how big does a sub string L have to be for it to be >99.5% certain that duplication wasn't due to chance
##This is inefficient and only works for sequence lengths up to 10000. Above that values are too large to compute
# in factorial functions
# Above 10KB the possion distribution is used as an approximation of the binomial

def framesize(SeqL):
 y = 0
 co_prob = 1
 
 if args.probability:
    cusprob = args.probability
    while co_prob > cusprob:
          y = y+1
          co_prob = pcircular(SeqL,y)
 else:           
    while co_prob > 0.005:
           y = y+1 
           co_prob = pcircular(SeqL,y)
 print('Probability of repeat ',co_prob)
 return y 
        
 ##Script will return index values of duplicated elements in a list, modified from stack overflow user 'PaulMcG'
##Two methods:
#def list_duplicates_of(seq,item):
    #start at -1 nice way to include a +1 term in a loop
    #while True will loop indefinitely, in this case untill the seq.index returns an error
    #start_at = loc means that when a duplicate is found, next time the seq.index runs it will consider all the elements after that element
    # but it will know their true index position
    #start_at = -1
    #locs = []
    #while True:
        #try:
            #loc = seq.index(item,start_at+1)
        #except ValueError:
            #break
        #else:
            #locs.append(loc)
            #start_at = loc
    #return locs
    
#numpy useful scientific package, esp. for linear algebra etc. 
#nonzero - find all instances of conditions that are true
#flat nonzero - give terms in 'flattened form' ie it levels out a nested list ([a,b],[c,d]) = (a,b,c,d)
#this function will find all elements > than limit and set them to 0 
def limit_size(size,limit,method):
 import numpy as np
 x = np.array(size)
 if method == 1:
    a = list(np.flatnonzero(x>limit))
 if method == 0:
    a = list(np.flatnonzero(x<limit)) 
 #print(a)
 x[a] = 0
 x = np.array(x).tolist()
 return x        

#source = "AKTGRGRGEDSJAADSJSJWFAFDFAFK"
#print(list_duplicates_of(source, 'B'))

##Using defaultdict this will find all duplicated elements in a single stroke, and return them and their index positions.
#Default dict makes a dict in the form (repeated element,[index_1, index_2]). Default dict is used in place of dict, as default dict
#will not give an error if it encounters an element not in the dictionary, instead it saves it as a new key. 

#Using list as the default_factory, it is easy to group a sequence of key-value pairs into a dictionary of lists.
#enumerate is being used here to give both the count and the current item being enumerated, this will give the index for each item. 
#tally is a default dict, so when an item is repeated it is saved to the same key and updates the (i) element.
# the return statement returns key locs, produced by the for command. the for command only rpduces key,locs if len(locs) >1 i.e. if key is repeated.

#I've added lines to deduce the size of the contig

from collections import defaultdict

def list_duplicates(seq):
    #size = []
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items()
                            if len(locs)>1)

def largest_repeat(seq): 
    size = []
    #dup = []
    for key,locs in list_duplicates(seq):
        r = locs[1] - locs[0]
        size.append(r)
    if args.maximum_size:
        size = limit_size(size,args.maximum_size,1)
    if args.minimum_size:
        size = limit_size(size,args.minimum_size,0)
    max_value = max(size)
    if max_value == 0:
       max_value = []
    max_index = size.index(max_value)
    return max_index



# We don't need to worry about reading in every frame
#list_duplicates will have duplicates ordered by index pos of first apperance, if two sequences are same size, always take the lowest 
#this will happen autmatically with this script

##To add - can pass a sequence or a list of sequences
# can pass a minimum size of genome to circularise [optional - chance co-prob criteria untill a size of this genome is made]
# can pass a mandatory probability
# ln means it can deal with sequences of any length
#we need sp-1 because the length command gives true length i.e. two elements = 2 not 1(0,1)


ln = 1*10**7
#fasta_file = input("File name?")
with open(args.fasta_file, "rt") as data: #with open automatically calls file close
    #The enumerate command gives an index for where the > is found 
    text = data.read()
    sp = ([i for i, x in enumerate(text) if x == ">"]) #> demarks new sequences
    sequence_name = ['Blank']*len(sp)
    sequence = ['Blank']*len(sp)
    #print(sp)
    ##if there is only one sequence
    if sp == [0]: 
        query = text
        #print(query)
        sequence_name = query.splitlines()[0]
        #print(sequence_name)
        sequence = query.splitlines()[1:ln] #obviously there are fewer lines than the ln but this works
        #sequence = sequence[0] #removes double list
        #print(sequence)

   ##if there is more than one
    else:
        for ind in range(0,len(sp)-1):
            query = (text[sp[ind]:sp[ind+1]])
            #print(query)
            sequence_name[ind] = query.splitlines()[0]
            #print(sequence_name)
            sequence[ind] = query.splitlines()[1:ln] #obviously there are fewer lines than the len(sequence) but this works
            #print(sequence)
            sequence[ind] = sequence[ind][0] #removes double list
   
        
        query = (text[sp[len(sp)-1]:len(text)]) #deal with the last sequence outside of the for loop
        sequence_name[len(sp)-1] = query.splitlines()[0]
        sequence[len(sp)-1] = query.splitlines()[1:len(sequence)]
        sequence[len(sp)-1] = sequence[len(sp)-1][0]
        #print(query)
        #print(sequence_name)
        #print(sequence)

#if sp == 0 sequence must be contained in a list otherwise DNA read as first integer
circgenome = ['Blank']*len(sp)
count = range(0,len(sp))

for c,DNA in zip(count,sequence):
   print(sequence_name[c])
    #print(DNA)
        # sindex and gindex consider only the first and second element of the repeat
#DNA = ('CCCCTATTTTTCATGCCAATAGAGAGTAGTAGTAGGAAGAAAATCAAGAACATCCTCTGAAGAAGCTGAATTGCCAACAATTAAATCCATGATGTAAACATCACCACACCCTGGCTTCGCAGTGGTTGAGAAAGGAGATGAGGTCATCTGTCCACCAATTTCCTGATCGTTATAGTGGAGATTGCGATTGACTCCGTGCCATCGCTTGTAAGTCTGTACAACTCCGGTGTCGTTGCCTGATCGGATGTTTGTGACTTTGTCGTGCATGACCGAGATACGTGCCGTGTCGATAGGAGCAGTGATTGGATCAATCCAGTCGCGAGGTGTCGCACTGACATTATTGATGCCGAAGCCGCGAAACATGAACGCGTACACGTCCGCGACTAGTGGGAAAGGAAGAGCTGTGTTGACACGTTGATATTCAACCGTGCCGGTGCCATCGTCGATCTGAGAGAAAACGCGGAGTATGTCCGAAGTATCCGTGAAACCAGGAGGTAGACCTTTGAGTGTGAAGACGACGCGTCGCCAACGCCATGCATTACTAGTGTCCGTCCGGATGGTGATTGTCTCTTTGACTCCACGAATGAAGGGAGCTTGGCTGAGACGTTGGGTGTTGATGGGCATCGTGAATGCCGCGGTGGGAGTATCGGACAGGGATCTTGCCGTCGCGCACCACAAGACCAATGCCGGAACATCCGAGGTCACCGAGAGGGCACCAACGTTCGTAGGTAGTGGAGGAAACGTGTTTCCCGGAACCATCGTGTCCCGCTTTTTCGTGCTTGTCATGTTCAGGACCCTCTTTCGGCTCACCGAGCGTGGCGTCCTCCTTGTTGAGCGTCGGGGGCGGGATATTCTCTTTGTTCGAGTCGTCCTCCCATATCGTGGTGACCTTCGCGTACTTCGACGCGGGCGAGACCGCTTTCGGCCGTAAGCCATAATTTAATTTCGGGAAAAGAACTTGGCGGGAAAGCGGGGTAAAGGGGGTTTGATCTTGGCGTCTGCCTGGAGTGAAGTCGTCGTTTAACCCGTTTAGGAAGGGGGCGTCTCCGAGTGGCATGATGTTCGAGGTGAGGGGCAACAGGGGAGAGGAGAGAGTATTTATAGGGACGAGGTGTCCCTCTGTCCCTGGGCTATAATATTAGTTTGCCCAGGGACTTTTTTACAAATGCCGAATTTTGACATTCACTGTCGCTATGCTCTCATCACATACTCTCAGTGCGGAGATCTCTCCCCTACAGTCGTTGGAGAATTCTTTGAGGGCCGTGGATACAAGTCGATCATTGGACGAGAGAATCACGCGGACGGAGGCGTTCATCTACATTGCTTTGTCGACTTTGGAAGGAAGAGAAGGTTCCGGAGAGCTCGTTGCTTCGATATCGAAGGCCGTCACCCCAATATTGAGCCTTCTCGTGGAACACCAGAAAAGGGTTGGGACTATGCATGCAAGGACGGTGATGTCTGCTTTCAGTCCCTTGACCGTCCGGGGGAGAGCGGAGGAAGCAATGGCGGAACTCGTGATAAGTGGGCTGCGATCACGGGTGCGAGCGATCGAGAGTCGTTTTGGGATTTGGTCCATGAATTGGATCCAAAGAGCGCGGCTTGTTCTTTCACCCAACTTCAAAAGTACTGTGACTGGAAGTTCGCTCCTGTGCCTCCCGTCTATGCCACACCAGACGGAATCACTTTCGTCGGAGGAGATGTTGATGGAAGAGATGAATGGCTATTACAGTCTGGTATCGGAAGTGGAGAGGCACTCATAGGTTAGTTATTCCCATGGGCGGGGCGGAGGAAACTGTATAGTGGTTGCCGCATGTGCGCCTCGGGGGGACCCCCAGCCCCCCCCCTCCCTCGTTGCACAAGCTCACGTGTTTAGGCAGATGCATGTCTATATGCGTATACGGAGAGTCCAGAACCGGAAAAACACTATGGGCTCGATCTCTTGGCGAGCACATCTACTGTGTTGGGTTGGTTTCAGGAGATGAGTGTTTGAAGGCCCCAAACGCCGAGTACGCCGTATTCGACGATATACGTGGGGGGATTAAGTTTTTTCCTTCGTTCAAGGAATGGTTGGGTTGTCAGGCTTGGGTCACGGTTAAATGTCTTTACAGGGAGCCTAAATTGGTTAAGTGGGGTAAGCCATCAATTTGGTTGAGCAACACAGACCCGAGGGATTACATGGAGAACAGTGATATTGATTGGATGAACAAAAATTGTATTTTCGTGGACGTTAACGCCCCTATTTTTCATGCCAATAGAG')
   if args.repeat:
       f = args.repeat-1;
       error = 0
       while error == 0:
         SeqL = len(DNA)  
         f = f+1
         eventn = events(SeqL,f)
         repeats = ['Blank'] * eventn
         dupe = []
         for i in range(0,eventn):
               repeats[i] = DNA[i:f+i]
       
         for dup in list_duplicates(repeats):
               dupe.append(dup)
            
         try: 
               index = largest_repeat(repeats)
               #print(index)
         except ValueError:
               print ("Maximum repeat size exceeded")
               #circgenome[c] = ' '
               error = 1
               
               
         except TypeError:
               print ("Maximum repeat size exceeded")
               #circgenome[c] = ' '
               error = 1
               
             
         else:
               sindex = (dupe[index][1][0])
               gindex = (dupe[index][1][1])
               #print(sindex)
               #print(gindex)
               circgenome[c] = (DNA[sindex:gindex])
   
               print('Searching for repeats of minimum length',f)
               print('Circularising at',(dupe[index]))
               print('Circularised genome size is',(gindex - sindex))
     
   else:
       SeqL = len(DNA)  
       f = framesize(SeqL)
       print('Searching for repeats of minimum length',f) 
       eventn = events(SeqL,f)
       repeats = ['Blank'] * eventn
       dupe = []
       for i in range(0,eventn):
             repeats[i] = DNA[i:f+i]
       
       for dup in list_duplicates(repeats):
             dupe.append(dup)
            
       try: 
             index = largest_repeat(repeats)
             #print(index)
       except ValueError:
             print ("Oops!  This genome couldn't be circularised [Try increasing the probability]")
             circgenome[c] = ' '
       
       except TypeError:
             print ("Oops!  This genome couldn't be circularised [Try changing the min max boundaries]")
             circgenome[c] = ' '    
             
       else:
             sindex = (dupe[index][1][0])
             gindex = (dupe[index][1][1])
             #print(sindex)
             #print(gindex)
             print('Circularising at',(dupe[index]))
             print('Circularised genome size is',(gindex - sindex))
             circgenome[c] = (DNA[sindex:gindex])
        
#print(circgenome)
with open(args.output_fasta_file,'w') as output:

     for name,genome in zip(sequence_name,circgenome):
         #for f in range(0, len(circgenome)):
             output.write(name)
             output.write('\n')   
             output.write(genome)
             output.write('\n')   
             output.write('\n')     
    
    
    