# cpi-em
ChIP-seq Peak Intensity - Expectation Maximization (CPI-EM) algorithm for detecting cooperatively bound transcription factor pairs. Paper reference to be updated soon. 

# Manual
## Dependencies

cpiEm requires ```scipy```, ```numpy```, ```pandas``` and ```argparse``` modules to be installed. cpiEm runs with both python 2.7+ and python 3+. 

## Running cpi-em

The syntax for running cpiEm.py -- 

    ./cpiEm.py -h

    usage: cpiEm.py [-h] --input-file <input file> [--mixture <dist>]
                [--output-file <output file>] [--num-iterations <NUM>]

    optional arguments:
    -h, --help            show this help message and exit
    -i <input file>, --input-file <input file>
                        Name of input file. Input should be a tab-separated
                        file with two columns consisting of primary TF peak
                        intensities (first column) and partner TF peak
                        intensities (second column). The partner TF is assumed
                        to be TF that cooperatively aids in binding the
                        primary TF to DNA.
    -m <dist>, --mixture <dist>
                        Distribution to be used to model peak intensities. Can
                        be lognormal, gaussian or gamma. Default : lognormal
    -o <output file>, --output-file <output file>
                        File to which output will be written. Default : <input
                        file>.cpi-em
    -n <NUM>, --num-iterations <NUM>
                        Maximum number of EM iterations to execute. Default :
                        10000

##Sample run 

This is a sample run of the cpi-em algorithm, using the `test-input` file in the repository -- 

    ./cpiEm.py -i test-input -m gamma -n 100
    
Output :     

    Fitting a gamma mixture to data
    Maximum number of EM iterations to run : 1000
    EM converged in 134 iterations. Q_2 = -3682.39982174
    Writing output to test-input.cpi-em
    The third column is the probability of a cooperative interaction being present for each entry in the input file
    
The first five lines of `test-input` are :
 
    3.59342	20.07449
    9.98265	6.6389
    7.50183	2.82852
    14.07073	7.64249
    3.34819	3.8665
    
The first five lines of the output file, `test-input.cpi-em`, are : 

    3.59342	20.07449	0.527476792074
    9.98265	6.6389	0.57156824552
    7.50183	2.82852	0.519679135763
    14.07073	7.64249	0.557642531493
    3.34819	3.8665	0.430667839721

