# CSPpipeline

## Steps to Run Blast
1. Open google colab - https://colab.research.google.com/
2. Login with your gmail account to enable google drive access. We will use google drive to store and run our python script. All the results will be saved to google drive
3. Create the following two folders
    - output
    - tmp
4. Create an input file (input.txt) with a list of accession numbers
5. Copy paste the "blast.py" file to colab and save it
6. Run blast

## Example output
### Forward Blast AZA18259.1 (Limit to Bacteria)
| Percentage | Subject ID |E_Value|
|--|--|--|
100.000| AZA18259.1| 0.0|
94.540| WP_084828638.1*| 0.0|
91.379| WP_111679482.1| 0.0|
76.012| WP_155214194.1| 0.0|
75.434| WP_155214166.1| 0.0|
**The top hit = WP_084828638.1**

---------------------------------------
### Reverse Blast WP_084828638.1 (Limit to Viruses)
| Percentage | Subject ID |E_Value|
|--|--|--|
95.690|YP_001686804.1|0.0|
95.690|ARU14331.1|0.0|
95.115|NP_056680.1|0.0|
94.540|AZA24404.1|0.0|
95.402|ARU13760.1|0.0|
94.253|YP_003344853.1|0.0|
95.115|ARU14608.1|0.0|
94.253|AZF92090.1|0.0|
95.115|AYP29589.1|0.0|
94.253|AYP30036.1|0.0|
**The top hit =  YP_001686804.1**