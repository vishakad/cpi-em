# cpi-em
ChIP-seq Peak Intensity - Expectation Maximization (CPI-EM) algorithm for detecting cooperatively bound transcription factor pairs. Paper reference to be updated soon. 

# Manual
## Running cpi-em

The syntax for running cpiEm.py -- 

    ./cpiEm.py <input file name> [<mixturetype> <max EM iterations> <output file name>]')

where, 
`<input file name>` (Mandatory argument) The file name input should point to a file that is a tab-separated file with two columns, which consist of primary TF peak intensities (first column) and partner TF peak intensities (second column). The partner TF is assumed to be TF that cooperatively aids in binding the primary TF to DNA.
`<mixturetype>` (Optional argument) can be `lognormal`, `gaussian` or `gamma`. Default : `lognormal`
`<max EM iterations>` (Optional argument) is the maximum number of EM iterations to execute. Default : `10000`
`<output file name>` (Optional argument) is the file to which output will be written. Default : `<input file name>.cpi-em`

The output of `cpiEm` consists of a three column tab-separated file, with the third column containing the probability that the partner TF helped the primary TF bind DNA at that location. The first two columns are the primary and partner TF peak intensities that were provided as the input.

##Sample run 

This is a sample run of the cpi-em algorithm, using the `test-input` file in the repository -- 

    ./cpiEm.py test-input gamma 1000
    
Output :     

    Fitting a gamma mixture to data
    Maximum number of EM iterations to run : 1000
    EM converged in 134 iterations. Q_2 = -3682.39982174
    Writing output to test-input.cpi-em
    The third column is the probability of a cooperative interaction being present for each entry in the input file
    
The first five lines of `test-input` are :

    4.8744	4.83374
    4.8744	4.44119
    4.8744	14.69983
    4.39167	14.69983
    10.49975	34.43588
    
The first five lines of the output file, `test-input.cpi-em`, are : 

    4.8744	4.83374	0.00156170868692
    4.8744	4.44119	0.00160993743759
    4.8744	14.69983	0.00104735591592
    4.39167	14.69983	0.00117815501185
    10.49975	34.43588	0.000425711313626
