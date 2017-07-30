# TE-NGS
TE-NGS, genomics tools for TE-enriched next generation sequencing and detection of active transposable elements

# Quick start  
TE-NGS is implemented in R and distributed using [packrat](https://github.com/rstudio/packrat) management system for R packages and their dependencies. Currently uses R version 3.1.2.  

Navigate to /scripts where packrat snapshot is built   
 
    $ cd scripts

First time launching R, packrat will build the local library of required packages and dependencies  

Ensure packrat is on  

    > packrat::on()

Check that packrat is up to date 
 
    > packrat::status()

Run the script in test mode: compare output provided as expected in /test, eg. 
    
    $ R --no-restore --no-save --no-readline --quiet < R_get_TE_calls_v0.1.r --args ../test/ ../test/seq_opt_table.txt 201 test

Run on your sample of choice:  

    $ R --no-restore --no-save --no-readline --quiet < R_get_TE_calls_v0.1.r --args ../foo/ ../foo/seq_opt_table.txt 201 sample  

where /foo points to your own bams generated from the TE-enriched [procotol](#)  

See the [project page](https://ekviky.github.io/TE-NGS/) for more information 
