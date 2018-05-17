# Simple-Circularise
For a sequence or list of sequences, simple-circularise will circularise a linear sequence that should be circular but was overlapped during assembly.

_Includes python script and geneious plugin extension._

**Geneious Plugin Installer:** Install the function as a geneious plugin on windows. The plugin will appear as 'Simple Circularise' in the 'Tools' menu on Geneious. When running for the first time running you will need to specify the location of the .exe file. It will be present in your geneious plugins folder which is in C:\ProgramData\Geneious\plugins on Windows 7 or C:\Users\AppData\Local\Geneious on Windows 10.

**Usage**: 
```shell 
!python simple_circularise.py [input.fasta] [output.fasta] [-p] [-min] [-max] [-r]
```
```[-min]```: set a minimum size for the output sequence. 

```[-max]```: set a maximum size for the output sequence.

```[-p]```: set the probability threshold for determining repeat size (default 0.005). 

```[-r]```: change the behaviour of the program to maximise repeat size (default is to maximise output sequence size). 


**Examples**:

```shell
!python simple_circularise.py linear.fasta circular.fasta -p 0.0001 -min 1000 -max 2000
``` 
Circularise a genome with a repeat size determined by a probability of co-occurance < 0.0001. Output the largest sequence between 1 - 2kb  that can be circularised.

```shell
!python simple_circularise.py linear.fasta circular.fasta -r 10 min 1000
```

Circularise a genome using largest possible repeat size. Start searching at repeat size 10 and increase until largest is found. Ensure output sequence is >1kb.
