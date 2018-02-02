# TE-NGS
Genomics tools for TE-enriched next generation sequencing and detection of active transposable elements

# Quick start  
TE-NGS is implemented in R and distributed using [packrat](https://github.com/rstudio/packrat) management system for R packages and their dependencies. Currently uses R version 3.1.2 and packrat is built using packrat 0.4.8-1  

Navigate to /distrib and unbundle the packrat tarball:    
 
    > packrat::unbundle(bundle='scripts-2017-07-31.tar.gz', where='/foo/TE-NGS/distrib')  
where /foo is wherever the TE-NGS repo lives locally 

Should see progress message:   

    Untarring 'scripts-2017-07-31.tar.gz' in directory '/foo/TE-NGS/distrib'...  

    ...  

    Done! The project has been unbundled and restored at:  
    - "/foo/TE-NGS/distrib/scripts"  

Navigate to /scripts where packrat snapshot is built:  
    
    $ cd scripts  


Check that packrat is up to date:  

    > packrat::status()

## Testing  

Run the script in test mode to check that output is as expected for a test sample provided in /distrib/test:    

    $  R --no-restore --no-save --no-readline --quiet < R_get_TE_calls_v0.1.r --args ../test/ ../test/seq_opt_table.txt 201 test  


Run the script in sample mode: on your sample of choice  

    $ R --no-restore --no-save --no-readline --quiet < R_get_TE_calls_v0.1.r --args /bar/ /bar/seq_opt_table.txt 201 sample  

where /bar points to your own directory of bams generated from the TE-enriched [protocol](http://rdcu.be/F6w6)  

## Documentation  

See the [project page](https://ekviky.github.io/TE-NGS/) for more information 

# Citing  
Kvikstad et al. *BMC Genomics* (2018) 19:115 
