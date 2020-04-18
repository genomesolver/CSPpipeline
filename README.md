# CSPpipeline

## Steps to run blast.py
1. Open google colab - https://colab.research.google.com/
2. Login with your gmail account to enable google drive access. We will use google drive to store and run our script.
3. Create a new notebook
4. Create the following two folders
    - output
    - tmp
5. Create an input file (input.txt) with a list of accession numbers
6. Copy paste the "blast.py" file to colab and save your notebook
7. Run blast

## Example terminal output
### Forward Blast QBX23717.1 (Limit to Bacteria)
|Query Cover % | E_Value | Accession Id | Subject Name|
|--|--|--|--|
99.72|0|WP_136022800|Streptococcus pyogenes|
99.72|0|WP_111679348|Streptococcus dysgalactiae|
99.72|0|WP_138083809|Streptococcus dysgalactiae|
99.72|0|WP_136284882|Streptococcus pyogenes|
99.72|0|WP_032460164|Streptococcus pyogenes|
99.72|0|WP_110410812|Streptococcus pyogenes|
99.72|0|WP_136263979|Streptococcus pyogenes|
99.72|0|WP_136298220|Streptococcus pyogenes|
96.36|0|WP_136037644|Streptococcus pyogenes|
96.36|0|WP_136132953|Streptococcus pyogenes|
**TOP HIT = WP_136022800 Streptococcus pyogenes**

---------------------------------------
### Reverse Blast WP_084828638.1 (Limit to Viruses)
|Query Cover % | E_Value | Accession Id | Subject Name|
|--|--|--|--|
99.72|0|WP_136022800|Streptococcus pyogenes|
99.72|0|QBX23717|Streptococcus phage Javan146|
99.72|0|WP_136112598|Streptococcus pyogenes|
95.84|0|WP_002990034|Streptococcus pyogenes|
95.84|0|WP_010922218|Streptococcus pyogenes|
99.72|0|WP_014635509|Streptococcus pyogenes|
95.84|0|WP_136040266|Streptococcus pyogenes|
95.84|0|WP_020833526|Streptococcus pyogenes|
95.84|0|WP_046177708|Streptococcus dysgalactiae|
99.72|0|WP_136111941|Streptococcus pyogenes|
**TOP HIT = WP_136022800 Streptococcus pyogenes**